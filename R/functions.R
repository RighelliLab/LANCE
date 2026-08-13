#' LANCE: tools for lncRNA single-cell analysis
#'
#' Functions for single-cell RNA-seq preprocessing, lncRNA specificity analysis,
#' pseudobulk summarization, and lncRNA-centered gene regulatory network
#' construction extracted from the original `new_pipeline/` analysis workflow.
#'
#' @keywords internal
"_PACKAGE"

#' Load tab-delimited counts into a Seurat object
#'
#' Read a gene-by-cell count matrix and matching cell metadata from tab-delimited files and return a quality-controlled Seurat object.
#'
#' @param counts_file Path to a tab-delimited gene-by-cell count matrix.
#' @param meta_file Path to a tab-delimited metadata table with cell barcodes in the row names.
#' @param mt_pattern Regular expression used to identify mitochondrial genes.
#' @param mt_threshold Maximum allowed mitochondrial percentage retained after filtering.
#' @param min_cells Minimum number of cells in which a feature must be detected.
#' @param min_features Minimum number of detected features required per cell.
#' @param save_path Optional path used to save the resulting Seurat object with `saveRDS()`.
#'
#' @return A Seurat object containing the filtered counts matrix and aligned metadata.
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
load_txt_counts_to_seurat <- function(
    counts_file,
    meta_file,
    mt_pattern = "^MT-",
    mt_threshold = 5,
    min_cells = 1,
    min_features = 1,
    save_path = NULL
) {

  if (!file.exists(counts_file)) {
    stop("Counts file not found: ", counts_file)
  }

  if (!file.exists(meta_file)) {
    stop("Metadata file not found: ", meta_file)
  }

  counts_connection <- if (grepl("\\.gz$", counts_file)) {
    gzfile(counts_file, open = "rt")
  } else {
    file(counts_file, open = "rt")
  }

  meta_connection <- if (grepl("\\.gz$", meta_file)) {
    gzfile(meta_file, open = "rt")
  } else {
    file(meta_file, open = "rt")
  }

  on.exit(close(counts_connection), add = TRUE)
  on.exit(close(meta_connection), add = TRUE)

  counts <- read.table(
    counts_connection,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  meta <- read.table(
    meta_connection,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # The input files are already organized as genes by cells.
  if (anyDuplicated(rownames(counts))) {
    stop("Gene names in the counts matrix are not unique")
  }

  if (anyDuplicated(colnames(counts))) {
    stop("Cell barcodes in the counts matrix are not unique")
  }

  if (anyDuplicated(rownames(meta))) {
    stop("Cell barcodes in the metadata are not unique")
  }

  # Standardize barcodes following the behavior of the original script.
  colnames(counts) <- gsub(
    "-",
    ".",
    colnames(counts),
    fixed = TRUE
  )

  rownames(meta) <- gsub(
    "-",
    ".",
    rownames(meta),
    fixed = TRUE
  )

  if (!setequal(colnames(counts), rownames(meta))) {

    missing_in_meta <- setdiff(
      colnames(counts),
      rownames(meta)
    )

    missing_in_counts <- setdiff(
      rownames(meta),
      colnames(counts)
    )

    stop(
      "Counts and metadata barcodes do not match. ",
      "Missing in metadata: ",
      length(missing_in_meta),
      "; missing in counts: ",
      length(missing_in_counts)
    )
  }

  # Reorder metadata to match the count-matrix columns.
  meta <- meta[colnames(counts), , drop = FALSE]

  dge <- Matrix::Matrix(
    as.matrix(counts),
    sparse = TRUE
  )

  obj <- Seurat::CreateSeuratObject(
    counts = dge,
    min.cells = min_cells,
    min.features = min_features
  )

  # CreateSeuratObject may discard genes or cells based on thresholds:
  # realign metadata to the cells that were actually retained.
  meta <- meta[colnames(obj), , drop = FALSE]

  obj <- Seurat::AddMetaData(
    object = obj,
    metadata = meta
  )

  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    obj,
    pattern = mt_pattern
  )

  obj <- subset(
    obj,
    subset = percent.mt < mt_threshold
  )

  if (!is.null(save_path)) {
    saveRDS(obj, save_path)
  }

  return(obj)
}

# ------------------------------------------------------------------------------
# FUNCTION: load_h5ad_to_seurat
# ------------------------------------------------------------------------------
# Load a normalized AnnData (.h5ad) dataset and convert it into a Seurat object.
#
# h5ad files used in this project contain:
# - normalized expression values in ad$X
# - cell metadata in ad$obs
#
# The expression matrix is stored in both:
# - counts: required by CreateSeuratObject / Seurat object construction
# - data: preserves the already-normalized values for downstream analyses
#
# No additional mitochondrial filtering or normalization is performed here.
# ------------------------------------------------------------------------------

#' Load a normalized AnnData (.h5ad) dataset and convert it into a Seurat object
#'
#' h5ad files used in this project contain:
#' - normalized expression values in ad$X
#' - cell metadata in ad$obs
#' The expression matrix is stored in both:
#' - counts: required by CreateSeuratObject / Seurat object construction
#' - data: preserves the already-normalized values for downstream analyses
#' No additional mitochondrial filtering or normalization is performed here.
#'
#' @param h5ad_file Path to an `.h5ad` file containing normalized expression values and cell metadata.
#' @param save_path Optional path used to save the resulting Seurat object with `saveRDS()`.
#' @param rename_celltype Logical indicating whether a `Celltype` column should be renamed to `Cell_type`.
#'
#' @return A Seurat object created from the AnnData matrix and metadata.
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
load_h5ad_to_seurat <- function(
    h5ad_file,
    save_path = NULL,
    rename_celltype = TRUE
) {
  
  if (!file.exists(h5ad_file)) {
    stop("h5ad file not found: ", h5ad_file)
  }
  
  if (!requireNamespace("anndataR", quietly = TRUE)) {
    stop("Package 'anndataR' is required to read h5ad files")
  }
  
  ad <- anndataR::read_h5ad(h5ad_file)
  
  if (is.null(ad$X)) {
    stop("AnnData object does not contain X")
  }
  
  # AnnData: cells x genes -> Seurat: genes x cells
  expr <- t(ad$X)
  
  # Convert row-oriented sparse matrices (dgRMatrix) to column-oriented
  # sparse matrices, which are handled more reliably by Seurat.
  if (inherits(expr, "dgRMatrix")) {
    expr <- methods::as(expr, "TsparseMatrix")
    expr <- methods::as(expr, "CsparseMatrix")
  } else if (!inherits(expr, "CsparseMatrix")) {
    expr <- methods::as(expr, "CsparseMatrix")
  }
  
  # Recover gene and cell names from AnnData metadata when needed.
  meta <- as.data.frame(ad$obs)
  gene_meta <- as.data.frame(ad$var)
  
  cell_names <- rownames(meta)
  gene_names <- rownames(gene_meta)
  
  if (is.null(cell_names) || length(cell_names) != ncol(expr)) {
    stop("Cell names in ad$obs are missing or inconsistent with ad$X")
  }
  
  if (is.null(gene_names) || length(gene_names) != nrow(expr)) {
    stop("Gene names in ad$var are missing or inconsistent with ad$X")
  }
  
  rownames(expr) <- gene_names
  colnames(expr) <- cell_names
  
  if (anyDuplicated(rownames(expr))) {
    stop("Gene names in the h5ad expression matrix are not unique")
  }
  
  if (anyDuplicated(colnames(expr))) {
    stop("Cell barcodes in the h5ad expression matrix are not unique")
  }
  
  if (anyDuplicated(rownames(meta))) {
    stop("Cell barcodes in ad$obs are not unique")
  }
  
  # Reorder metadata exactly as the expression matrix.
  meta <- meta[colnames(expr), , drop = FALSE]
  
  # Standardize the metadata name used throughout the current pipeline.
  if (rename_celltype &&
      "Celltype" %in% colnames(meta) &&
      !"Cell_type" %in% colnames(meta)) {
    colnames(meta)[colnames(meta) == "Celltype"] <- "Cell_type"
  }
  
  obj <- Seurat::CreateSeuratObject(
    counts = expr,
    assay = "RNA",
    meta.data = meta
  )
  
  # ad$X is already normalized in the PDAC/LU-ADC inputs used in this project.
  obj <- SeuratObject::SetAssayData(
    object = obj,
    assay = "RNA",
    layer = "data",
    new.data = expr
  )
  
  if (!is.null(save_path)) {
    saveRDS(obj, save_path)
  }
  
  return(obj)
}

# ------------------------------------------------------------------------------
# FUNCTION: merge_seurat_objects
# ------------------------------------------------------------------------------
# Add dataset identifiers to Seurat objects and merge them into a single object.
#
# INPUT:
# - obj_list: named list of Seurat objects
# - save_path: optional path for saving merged object
#
# OUTPUT:
# - merged Seurat object
# ------------------------------------------------------------------------------

#' Add dataset identifiers to Seurat objects and merge them into a single object
#'
#' Add dataset identifiers to Seurat objects and merge them into a single object.
#'
#' @param obj_list named list of Seurat objects
#' @param save_path optional path for saving merged object
#'
#' @return merged Seurat object
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
merge_seurat_objects <- function(
    obj_list,
    save_path = NULL
) {

  if (length(obj_list) < 2L) {
    stop("obj_list must contain at least two Seurat objects")
  }
  if (is.null(names(obj_list)) || any(!nzchar(names(obj_list))) ||
      anyDuplicated(names(obj_list))) {
    stop("obj_list must have unique, non-empty names")
  }

  # Add dataset identifier
  for (name in names(obj_list)) {
    obj_list[[name]]$GSE <- name
  }

  # Prefix cell barcodes to preserve dataset provenance.
  merged_obj <- merge(
    x = obj_list[[1]],
    y = obj_list[-1],
    add.cell.ids = names(obj_list),
    project = "merged"
  )

  # Save if requested
  if (!is.null(save_path)) {
    saveRDS(merged_obj, save_path)
  }

  return(merged_obj)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_cell_types
# ------------------------------------------------------------------------------
# Visualize cell type annotations using DimPlot.
#
# INPUT:
# - obj: Seurat object
# - group.by: metadata column for grouping
#
# OUTPUT:
# - ggplot object (DimPlot)
# ------------------------------------------------------------------------------

#' Visualize cell type annotations using DimPlot
#'
#' Visualize cell type annotations using DimPlot.
#'
#' @param obj Seurat object
#' @param group.by metadata column for grouping
#'
#' @return ggplot object (DimPlot)
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_cell_types <- function(obj, group.by) {

  cols <- scales::hue_pal()(length(unique(obj[[group.by]][,1])))

  Seurat::DimPlot(
    obj,
    group.by = group.by,
    cols = cols,
    label = TRUE,
    label.box = TRUE,
    label.color = "white"
  ) +
    ggplot2::theme(legend.position = "none")
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_scrmra_input
# ------------------------------------------------------------------------------
# Prepares normalized expression matrix for scMRMA annotation.
#
# INPUT:
# - obj: Seurat object
#
# OUTPUT:
# - logcounts matrix
# ------------------------------------------------------------------------------

#' Prepares normalized expression matrix for scMRMA annotation
#'
#' Prepares normalized expression matrix for scMRMA annotation.
#'
#' @param obj Seurat object
#'
#' @return logcounts matrix
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
prepare_scrmra_input <- function(obj) {

  SeuratObject::DefaultAssay(obj) <- "RNA"
  if (!"data" %in% SeuratObject::Layers(obj[["RNA"]])) {
    obj <- Seurat::NormalizeData(obj)
  }
  obj <- SeuratObject::JoinLayers(obj)

  logcounts <- Seurat::GetAssayData(obj, layer = "data")

  return(list(obj = obj, logcounts = logcounts))
}

# ------------------------------------------------------------------------------
# FUNCTION: run_scrmra_annotation
# ------------------------------------------------------------------------------
# Runs scMRMA annotation and stores results in Seurat object.
#
# INPUT:
# - obj: Seurat object
# - logcounts: expression matrix
#
# OUTPUT:
# - updated Seurat object
# ------------------------------------------------------------------------------
#' Runs scMRMA annotation and stores results in Seurat object
#'
#' Runs scMRMA annotation and stores results in Seurat object.
#'
#' @param obj Seurat object
#' @param logcounts expression matrix
#'
#' @return updated Seurat object
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
run_scrmra_annotation <- function(obj, logcounts) {

  result <- scMRMA::scMRMA(
    input = logcounts,
    species = "Hs",
    normalizedData = TRUE
  )

  annotation_result <- result$multiR$annotationResult

  if (is.null(rownames(annotation_result)) || !length(rownames(annotation_result))) {
    stop("scMRMA annotationResult has no row names; cannot align annotations to cell barcodes.")
  }

  barcodes <- colnames(obj)
  missing_barcodes <- setdiff(barcodes, rownames(annotation_result))
  if (length(missing_barcodes) > 0L) {
    stop(paste0(
      "scMRMA annotationResult is missing annotations for ",
      length(missing_barcodes),
      " barcode(s), including: ",
      paste(head(missing_barcodes, 5L), collapse = ", ")
    ))
  }

  annotation_values <- annotation_result[barcodes, 1]
  if (anyNA(annotation_values)) {
    stop("scMRMA annotationResult matching produced NA values for one or more barcodes.")
  }

  obj$scMRMA_annotation <- annotation_values
  obj$cell_type1 <- annotation_values

  return(obj)
}

# ------------------------------------------------------------------------------
# FUNCTION: seurat_preprocess
# ------------------------------------------------------------------------------
# Standard preprocessing pipeline for a Seurat object.
#
# This function prepares the object for downstream analysis by:
# - setting the active assay to RNA
# - joining Seurat v5 layers into a unified matrix
# - performing log-normalization
# - defining cell identities based on the "Cell_type" metadata column
#
# IMPORTANT:
# This function does not modify biological content or filtering,
# it only standardizes the object structure for analysis.
#
# INPUT:
# - obj: Seurat object
#
# OUTPUT:
# - preprocessed Seurat object ready for downstream analysis
# ------------------------------------------------------------------------------

#' Standard preprocessing pipeline for a Seurat object
#'
#' This function prepares the object for downstream analysis by:
#' - setting the active assay to RNA
#' - joining Seurat v5 layers into a unified matrix
#' - performing log-normalization
#' - defining cell identities based on the "Cell_type" metadata column
#' IMPORTANT:
#' This function does not modify biological content or filtering,
#' it only standardizes the object structure for analysis.
#'
#' @param obj Seurat object
#'
#' @return preprocessed Seurat object ready for downstream analysis
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
seurat_preprocess <- function(obj) {
  
  SeuratObject::DefaultAssay(obj) <- "RNA"
  
  # Merge Seurat v5 layers first (e.g. counts.1/data.1 after dataset merge).
  obj <- SeuratObject::JoinLayers(
    object = obj,
    assay = "RNA"
  )
  
  available_layers <- SeuratObject::Layers(obj[["RNA"]])
  
  # Raw-count datasets do not yet have a normalized data layer.
  # h5ad datasets already carry normalized values in "data",
  # so they must not be normalized a second time.
  if (!"data" %in% available_layers) {
    obj <- Seurat::NormalizeData(obj)
  }
  
  if (!"Cell_type" %in% colnames(obj@meta.data)) {
    stop("Metadata column 'Cell_type' not found")
  }
  
  Seurat::Idents(obj) <- "Cell_type"
  
  return(obj)
}


# ------------------------------------------------------------------------------
# FUNCTION: gene_biotype_matching
# ------------------------------------------------------------------------------
# Identify lncRNA and protein-coding genes present in a Seurat object
# based on GENCODE annotation.
#
# The function:
# - splits GENCODE annotation into lncRNA and protein-coding genes
# - matches genes with those present in the Seurat object
# - identifies ambiguous genes present in both categories
#
# INPUT:
# - obj: Seurat object
# - gencode: data.frame with gene annotation (SYMBOL + type)
#
# OUTPUT:
# - list containing:
#   - lnc: lncRNA genes in dataset
#   - prot: protein-coding genes in dataset
#   - ambiguous: genes in both categories
#   - lnc_raw: raw lncRNA matches
#   - prot_raw: raw protein-coding matches
# ------------------------------------------------------------------------------

#' Identify lncRNA and protein-coding genes present in a Seurat object
#'
#' based on GENCODE annotation.
#' The function:
#' - splits GENCODE annotation into lncRNA and protein-coding genes
#' - matches genes with those present in the Seurat object
#' - identifies ambiguous genes present in both categories
#'
#' @param obj Seurat object
#' @param gencode data.frame with gene annotation (SYMBOL + type)
#'
#' @return list containing: lnc: lncRNA genes in dataset prot: protein-coding genes in dataset ambiguous: genes in both categories lnc_raw: raw lncRNA matches prot_raw: raw protein-coding matches
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
gene_biotype_matching <- function(obj, gencode) {
  
  lnc <- gencode[gencode$type == "lncRNA", ]
  prot <- gencode[gencode$type == "protein_coding", ]
  
  genes_obj <- rownames(obj)
  
  lnc_raw <- intersect(
    genes_obj,
    lnc$SYMBOL
  )
  
  prot_raw <- intersect(
    genes_obj,
    prot$SYMBOL
  )
  
  ambiguous_genes <- intersect(
    unique(lnc_raw),
    unique(prot_raw)
  )
  
  # In ambiguous cases, prioritize the lncRNA category.
  # Ambiguous genes are therefore retained among lncRNAs.
  lnc_features <- unique(lnc_raw)
  
  # Exclude the same genes from the protein-coding set,
  # so the final feature sets remain disjoint.
  prot_features <- setdiff(
    unique(prot_raw),
    ambiguous_genes
  )
  
  if (length(intersect(lnc_features, prot_features)) > 0L) {
    stop(
      "lncRNA and protein-coding feature sets are not disjoint"
    )
  }
  
  
  return(list(
    lnc = lnc_features,
    prot = prot_features,
    ambiguous = ambiguous_genes,
    lnc_raw = lnc_raw,
    prot_raw = prot_raw
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: compute_gini_dotplot
# ------------------------------------------------------------------------------
# Compute gene-level Gini coefficient using Seurat DotPlot output.
#
# The function:
# - extracts pct.exp from DotPlot
# - computes Gini coefficient per gene across cell types
# - returns tidy table of results
#
# INPUT:
# - obj: Seurat object
# - features_chunk: vector of genes to evaluate
#
# OUTPUT:
# - data.frame with gene and gini score
# ------------------------------------------------------------------------------

#' Compute gene-level Gini coefficient using Seurat DotPlot output
#'
#' The function:
#' - extracts pct.exp from DotPlot
#' - computes Gini coefficient per gene across cell types
#' - returns tidy table of results
#'
#' @param obj Seurat object
#' @param features_chunk vector of genes to evaluate
#' @param min_detection Minimum percentage of detected cells required for a gene to be retained.
#'
#' @return data.frame with gene and gini score
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
compute_gini_dotplot <- function(obj, features_chunk, min_detection = 0) {

  if (!is.numeric(min_detection) || length(min_detection) != 1L ||
      is.na(min_detection) || min_detection < 0 || min_detection > 100) {
    stop("min_detection must be a single number between 0 and 100")
  }
  features_chunk <- intersect(features_chunk, rownames(obj))
  if (!length(features_chunk)) {
    return(data.frame(gene = character(), gini = numeric()))
  }

  dp <- Seurat::DotPlot(obj, features = features_chunk)
  df <- dp$data[, c("features.plot", "pct.exp")]
  detection <- aggregate(pct.exp ~ features.plot, df, max)
  keep <- detection$features.plot[detection$pct.exp >= min_detection]
  df <- df[df$features.plot %in% keep, , drop = FALSE]

  # This is a detection-heterogeneity Gini, not an expression-intensity Gini.
  out <- aggregate(pct.exp ~ features.plot, df, function(x) DescTools::Gini(x))

  colnames(out) <- c("gene", "gini")
  return(out)
}

#' This is a detection-heterogeneity Gini, not an expression-intensity Gini
#'
#' This is a detection-heterogeneity Gini, not an expression-intensity Gini.
#'
#' @param obj Seurat object used to compute DotPlot-based detection summaries.
#' @param features Character vector of gene names to evaluate.
#' @param type Label assigned to the returned feature set.
#' @param chunk_size Number of genes processed in each chunk.
#'
#' @return A data frame with one row per gene and columns for `gene`, `gini`, and `type`.
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
compute_gini_pipeline <- function(obj, features, type, chunk_size = 1000) {

  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || is.na(chunk_size) || chunk_size <= 0) {
    stop("chunk_size must be a single positive number")
  }

  if (length(features) == 0L) {
    return(data.frame(
      gene = character(),
      gini = numeric(),
      type = character()
    ))
  }

  chunks <- split(features,
                  ceiling(seq_along(features) / chunk_size))

  res_list <- lapply(chunks, function(x) {
    compute_gini_dotplot(obj, x)
  })

  res <- do.call(rbind, res_list)
  res$type <- type

  return(res)
}

#' This is a detection-heterogeneity Gini, not an expression-intensity Gini
#'
#' This is a detection-heterogeneity Gini, not an expression-intensity Gini.
#'
#' @param GiniCoef_clean Cleaned Gini summary table without missing values.
#' @param lnc_n Number of lncRNA genes represented in the plot.
#' @param prot_n Number of protein-coding genes represented in the plot.
#'
#' @return A ggplot object comparing Gini distributions across gene classes.
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_gini_results <- function(GiniCoef_clean, lnc_n, prot_n) {
  
  ggplot2::ggplot(GiniCoef_clean, ggplot2::aes(x = type, y = gini, fill = type)) +
    ggplot2::geom_boxplot() +
    ggpubr::stat_compare_means(hjust = 0) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(
      values = c(
        "lncRNA" = "#0072B2",
        "protein_coding" = "#D55E00"
      )
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::annotate(
      x = c("lncRNA", "protein_coding"),
      y = c(0.84, 0.57),
      label = c(
        paste0("(n = ", lnc_n, ")"),
        paste0("(n = ", prot_n, ")")
      ),
      col = "black",
      size = 4,
      geom = "text"
    )
}

# ------------------------------------------------------------------------------
# FUNCTION: check_gene_summary
# ------------------------------------------------------------------------------
# Perform basic validation and expression summary for a single gene.
#
# This function:
# - checks whether a gene is present in the Seurat object
# - computes average expression per cell type (Seurat AverageExpression)
# - extracts DotPlot data for downstream Gini computation
#
# INPUT:
# - obj: Seurat object
# - gene: gene name (character)
#
# OUTPUT:
# - list containing:
#   - present: logical (gene presence in dataset)
#   - avg_expression: average expression per cell type
#   - dotplot_data: Seurat DotPlot output data
# ------------------------------------------------------------------------------

#' Perform basic validation and expression summary for a single gene
#'
#' This function:
#' - checks whether a gene is present in the Seurat object
#' - computes average expression per cell type (Seurat AverageExpression)
#' - extracts DotPlot data for downstream Gini computation
#'
#' @param obj Seurat object
#' @param gene gene name (character)
#'
#' @return list containing: present: logical (gene presence in dataset) avg_expression: average expression per cell type dotplot_data: Seurat DotPlot output data
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
check_gene_summary <- function(obj, gene) {

  present <- gene %in% rownames(obj)

  if (!present) {
    return(list(
      present = FALSE,
      avg_expression = NULL,
      dotplot_data = NULL
    ))
  }

  avg_expression <- Seurat::AverageExpression(obj, features = gene)$RNA

  dotplot_data <- Seurat::DotPlot(obj, features = gene)$data

  return(list(
    present = present,
    avg_expression = avg_expression,
    dotplot_data = dotplot_data
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: compute_gene_gini_single
# ------------------------------------------------------------------------------
# Compute Gini coefficient for a single gene using DotPlot data.
#
# This function:
# - extracts pct.exp values from Seurat DotPlot output
# - computes Gini coefficient across cell types
# - returns tidy table
#
# INPUT:
# - dotplot_data: output from Seurat DotPlot
# - gene: gene name
#
# OUTPUT:
# - data.frame with:
#   - celltype
#   - gini
# ------------------------------------------------------------------------------

#' Compute Gini coefficient for a single gene using DotPlot data
#'
#' This function:
#' - extracts pct.exp values from Seurat DotPlot output
#' - computes Gini coefficient across cell types
#' - returns tidy table
#'
#' @param dotplot_data output from Seurat DotPlot
#' @param gene gene name
#'
#' @return data.frame with: celltype gini
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
compute_gene_gini_single <- function(dotplot_data, gene) {

  df <- dotplot_data[, c("features.plot", "pct.exp")]

  out <- aggregate(pct.exp ~ features.plot, df, DescTools::Gini)

  colnames(out) <- c("gene", "gini")

  return(out)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_gene_density_umap
# ------------------------------------------------------------------------------
# Visualize the expression density of a gene on an existing UMAP.
#
# INPUT:
# - obj: Seurat object containing normalized data and a UMAP reduction
# - gene: gene to visualize
# - reduction: dimensional reduction used for the plot
#
# OUTPUT:
# - Nebulosa density plot
# ------------------------------------------------------------------------------

#' Visualize the expression density of a gene on an existing UMAP
#'
#' Visualize the expression density of a gene on an existing UMAP.
#'
#' @param obj Seurat object containing normalized data and a UMAP reduction
#' @param gene gene to visualize
#' @param reduction dimensional reduction used for the plot
#'
#' @return Nebulosa density plot
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_gene_density_umap <- function(
    obj,
    gene,
    reduction = "umap"
) {
  
  if (!gene %in% rownames(obj)) {
    stop(
      "Gene not found in the Seurat object: ",
      gene
    )
  }
  
  available_reductions <- SeuratObject::Reductions(obj)
  
  if (!reduction %in% available_reductions) {
    stop(
      "Reduction '",
      reduction,
      "' not found. Available reductions: ",
      paste(available_reductions, collapse = ", ")
    )
  }
  
  density_plot <- Nebulosa::plot_density(
    obj,
    features = gene,
    reduction = reduction
  )
  
  return(density_plot)
}

# ------------------------------------------------------------------------------
# FUNCTION: validate_grn_metadata
# ------------------------------------------------------------------------------
# Validate metadata columns required for GRN analysis.
#
# INPUT:
# - obj: Seurat object
# - config: GRN configuration list
#
# OUTPUT:
# - data.frame reporting missing values for each required metadata column
# ------------------------------------------------------------------------------

#' Validate metadata columns required for GRN analysis
#'
#' Validate metadata columns required for GRN analysis.
#'
#' @param obj Seurat object
#' @param config GRN configuration list
#'
#' @return data.frame reporting missing values for each required metadata column
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
validate_grn_metadata <- function(obj, config) {
  
  required_config <- c(
    "cell_type_col",
    "epithelial_label",
    "patient_col",
    "sample_col",
    "dataset_col",
    "condition_col"
  )
  
  missing_config <- setdiff(required_config, names(config))
  
  if (length(missing_config) > 0L) {
    stop(
      "Missing GRN configuration entries: ",
      paste(missing_config, collapse = ", ")
    )
  }
  
  metadata <- obj@meta.data
  
  required_columns <- unlist(
    config[c(
      "cell_type_col",
      "patient_col",
      "sample_col",
      "dataset_col",
      "condition_col"
    )],
    use.names = FALSE
  )
  
  missing_columns <- setdiff(
    required_columns,
    colnames(metadata)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing metadata columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  cell_types <- unique(
    as.character(metadata[[config$cell_type_col]])
  )
  
  if (!config$epithelial_label %in% cell_types) {
    stop(
      "Epithelial label not found in metadata: ",
      config$epithelial_label
    )
  }
  
  missing_summary <- data.frame(
    column = required_columns,
    missing_values = vapply(
      required_columns,
      function(column) {
        sum(
          is.na(metadata[[column]]) |
            trimws(as.character(metadata[[column]])) == ""
        )
      },
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
  
  return(missing_summary)
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_epithelial_samples
# ------------------------------------------------------------------------------
# Summarize epithelial cell counts by dataset, patient, sample, and condition.
#
# INPUT:
# - obj: Seurat object
# - config: GRN configuration list
#
# OUTPUT:
# - data.frame with the number of epithelial cells in each sample
# ------------------------------------------------------------------------------

#' Summarize epithelial cell counts by dataset, patient, sample, and condition
#'
#' Summarize epithelial cell counts by dataset, patient, sample, and condition.
#'
#' @param obj Seurat object
#' @param config GRN configuration list
#'
#' @return data.frame with the number of epithelial cells in each sample
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
summarize_epithelial_samples <- function(obj, config) {
  
  metadata <- obj@meta.data
  
  epithelial_metadata <- metadata[
    as.character(metadata[[config$cell_type_col]]) ==
      config$epithelial_label,
    ,
    drop = FALSE
  ]
  
  if (nrow(epithelial_metadata) == 0L) {
    stop(
      "No cells found for: ",
      config$epithelial_label
    )
  }
  
  grouping_columns <- c(
    config$dataset_col,
    config$patient_col,
    config$sample_col,
    config$condition_col
  )
  
  if (
    "status_col" %in% names(config) &&
    config$status_col %in% colnames(epithelial_metadata)
  ) {
    grouping_columns <- c(
      grouping_columns,
      config$status_col
    )
  }
  
  grouped_metadata <- dplyr::group_by(
    epithelial_metadata,
    dplyr::across(dplyr::all_of(grouping_columns))
  )
  
  sample_summary <- dplyr::summarise(
    grouped_metadata,
    n_epithelial_cells = dplyr::n(),
    .groups = "drop"
  )
  
  sample_summary <- as.data.frame(sample_summary)
  
  sample_summary <- sample_summary[
    do.call(
      order,
      sample_summary[grouping_columns]
    ),
    ,
    drop = FALSE
  ]
  
  rownames(sample_summary) <- NULL
  
  return(sample_summary)
}

# ------------------------------------------------------------------------------
# FUNCTION: select_grn_cells
# ------------------------------------------------------------------------------
# Select cells required for GRN analysis and assign a pseudobulk group.
#
# INPUT:
# - obj: Seurat object
# - config: GRN configuration list
# - group_col: name of the new pseudobulk grouping column
#
# OUTPUT:
# - filtered Seurat object containing only the selected cells
# ------------------------------------------------------------------------------

#' Select cells required for GRN analysis and assign a pseudobulk group
#'
#' Select cells required for GRN analysis and assign a pseudobulk group.
#'
#' @param obj Seurat object
#' @param config GRN configuration list
#' @param group_col name of the new pseudobulk grouping column
#'
#' @return filtered Seurat object containing only the selected cells
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
select_grn_cells <- function(
    obj,
    config,
    group_col = "GRN_pseudobulk_id"
) {
  
  metadata <- obj@meta.data
  
  required_columns <- unique(c(
    config$cell_type_col,
    config$condition_col,
    config$pseudobulk_group_columns
  ))
  
  missing_columns <- setdiff(
    required_columns,
    colnames(metadata)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing metadata columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  cells_to_keep <- rownames(metadata)[
    as.character(metadata[[config$cell_type_col]]) ==
      config$epithelial_label &
      as.character(metadata[[config$condition_col]]) ==
      config$target_condition
  ]
  
  if (length(cells_to_keep) == 0L) {
    stop(
      "No cells found for cell type '",
      config$epithelial_label,
      "' and condition '",
      config$target_condition,
      "'"
    )
  }
  
  grn_obj <- obj[, cells_to_keep]
  
  selected_metadata <- grn_obj@meta.data
  
  invalid_group_values <- vapply(
    config$pseudobulk_group_columns,
    function(column) {
      any(
        is.na(selected_metadata[[column]]) |
          trimws(as.character(selected_metadata[[column]])) == ""
      )
    },
    logical(1)
  )
  
  if (any(invalid_group_values)) {
    stop(
      "Missing values found in pseudobulk grouping columns: ",
      paste(
        config$pseudobulk_group_columns[invalid_group_values],
        collapse = ", "
      )
    )
  }
  
  pseudobulk_id <- do.call(
    paste,
    c(
      selected_metadata[
        config$pseudobulk_group_columns
      ],
      sep = "__"
    )
  )
  
  grn_obj[[group_col]] <- pseudobulk_id
  
  return(grn_obj)
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_grn_counts
# ------------------------------------------------------------------------------
# Join the raw-count layers of the GRN Seurat object.
#
# INPUT:
# - obj: Seurat object containing the selected GRN cells
# - assay: assay containing raw counts
#
# OUTPUT:
# - Seurat object with a single joined counts layer
# ------------------------------------------------------------------------------

#' Join the raw-count layers of the GRN Seurat object
#'
#' Join the raw-count layers of the GRN Seurat object.
#'
#' @param obj Seurat object containing the selected GRN cells
#' @param assay assay containing raw counts
#'
#' @return Seurat object with a single joined counts layer
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
prepare_grn_counts <- function(
    obj,
    assay = "RNA"
) {
  
  if (!assay %in% SeuratObject::Assays(obj)) {
    stop("Assay not found in Seurat object: ", assay)
  }
  
  SeuratObject::DefaultAssay(obj) <- assay
  
  count_layers <- SeuratObject::Layers(
    obj[[assay]],
    search = "^counts"
  )
  
  if (length(count_layers) == 0L) {
    stop("No counts layer found in assay: ", assay)
  }
  
  if (length(count_layers) > 1L) {
    obj <- SeuratObject::JoinLayers(
      object = obj,
      assay = assay
    )
  }
  
  return(obj)
}

# ------------------------------------------------------------------------------
# FUNCTION: aggregate_grn_pseudobulk
# ------------------------------------------------------------------------------
# Aggregate raw gene counts by pseudobulk group.
#
# INPUT:
# - obj: Seurat object containing the selected GRN cells
# - assay: assay containing the raw counts
# - group_col: metadata column defining the pseudobulk groups
#
# OUTPUT:
# - sparse matrix with genes as rows and pseudobulk samples as columns
# ------------------------------------------------------------------------------

#' Aggregate raw gene counts by pseudobulk group
#'
#' Aggregate raw gene counts by pseudobulk group.
#'
#' @param obj Seurat object containing the selected GRN cells
#' @param assay assay containing the raw counts
#' @param group_col metadata column defining the pseudobulk groups
#'
#' @return sparse matrix with genes as rows and pseudobulk samples as columns
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
aggregate_grn_pseudobulk <- function(
    obj,
    assay = "RNA",
    group_col = "GRN_pseudobulk_id"
) {
  
  if (!assay %in% SeuratObject::Assays(obj)) {
    stop("Assay not found in Seurat object: ", assay)
  }
  
  if (!group_col %in% colnames(obj@meta.data)) {
    stop("Pseudobulk grouping column not found: ", group_col)
  }
  
  available_layers <- SeuratObject::Layers(
    obj[[assay]]
  )
  
  if (!"counts" %in% available_layers) {
    stop(
      "A single 'counts' layer was not found. ",
      "Run prepare_grn_counts() first."
    )
  }
  
  counts <- SeuratObject::LayerData(
    object = obj,
    assay = assay,
    layer = "counts"
  )
  
  pseudobulk_groups <- as.character(
    obj@meta.data[colnames(counts), group_col]
  )
  
  if (anyNA(pseudobulk_groups) || any(trimws(pseudobulk_groups) == "")) {
    stop("Missing values found in pseudobulk group assignments")
  }
  
  group_matrix <- Matrix::sparse.model.matrix(
    ~ 0 + factor(pseudobulk_groups)
  )
  
  colnames(group_matrix) <- levels(
    factor(pseudobulk_groups)
  )
  
  pseudobulk_counts <- counts %*% group_matrix
  
  return(pseudobulk_counts)
}

# ------------------------------------------------------------------------------
# FUNCTION: validate_grn_pseudobulk
# ------------------------------------------------------------------------------
# Validate the pseudobulk count matrix after aggregation.
#
# Checks:
# - matrix dimensions
# - pseudobulk sample names
# - duplicated gene/sample names
# - missing or negative values
# - consistency with expected pseudobulk groups
#
# INPUT:
# - pseudobulk_counts: genes x pseudobulk raw-count matrix
# - obj: Seurat object used for pseudobulk aggregation
# - group_col: metadata column defining pseudobulk groups
#
# OUTPUT:
# - data.frame summarizing validation results
# ------------------------------------------------------------------------------

#' Validate the pseudobulk count matrix after aggregation
#'
#' Checks:
#' - matrix dimensions
#' - pseudobulk sample names
#' - duplicated gene/sample names
#' - missing or negative values
#' - consistency with expected pseudobulk groups
#'
#' @param pseudobulk_counts genes x pseudobulk raw-count matrix
#' @param obj Seurat object used for pseudobulk aggregation
#' @param group_col metadata column defining pseudobulk groups
#'
#' @return data.frame summarizing validation results
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
validate_grn_pseudobulk <- function(
    pseudobulk_counts,
    obj,
    group_col = "GRN_pseudobulk_id"
) {
  
  if (!group_col %in% colnames(obj@meta.data)) {
    stop("Pseudobulk grouping column not found: ", group_col)
  }
  
  expected_groups <- sort(
    unique(as.character(obj@meta.data[[group_col]]))
  )
  
  observed_groups <- sort(
    colnames(pseudobulk_counts)
  )
  
  validation <- data.frame(
    check = c(
      "number_of_genes",
      "number_of_pseudobulk_samples",
      "expected_number_of_pseudobulk_samples",
      "duplicated_gene_names",
      "duplicated_pseudobulk_names",
      "missing_values",
      "negative_values",
      "group_names_match"
    ),
    value = c(
      nrow(pseudobulk_counts),
      ncol(pseudobulk_counts),
      length(expected_groups),
      anyDuplicated(rownames(pseudobulk_counts)),
      anyDuplicated(colnames(pseudobulk_counts)),
      sum(is.na(pseudobulk_counts)),
      sum(pseudobulk_counts < 0),
      identical(expected_groups, observed_groups)
    ),
    stringsAsFactors = FALSE
  )
  
  return(validation)
}

# ------------------------------------------------------------------------------
# FUNCTION: normalize_grn_pseudobulk
# ------------------------------------------------------------------------------
# Normalize pseudobulk raw counts using TMM normalization and log-CPM.
#
# INPUT:
# - pseudobulk_counts: genes x pseudobulk raw-count matrix
#
# OUTPUT:
# - list containing:
#   - logcpm: normalized log-CPM expression matrix
#   - norm_factors: TMM normalization factors
#   - library_sizes: raw library sizes
# ------------------------------------------------------------------------------

#' Normalize pseudobulk raw counts using TMM normalization and log-CPM
#'
#' Normalize pseudobulk raw counts using TMM normalization and log-CPM.
#'
#' @param pseudobulk_counts genes x pseudobulk raw-count matrix
#'
#' @return list containing: logcpm: normalized log-CPM expression matrix norm_factors: TMM normalization factors library_sizes: raw library sizes
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
normalize_grn_pseudobulk <- function(pseudobulk_counts) {
  
  if (anyNA(pseudobulk_counts)) {
    stop("Pseudobulk count matrix contains missing values")
  }
  
  if (any(pseudobulk_counts < 0)) {
    stop("Pseudobulk count matrix contains negative values")
  }
  
  y <- edgeR::DGEList(
    counts = pseudobulk_counts
  )
  
  y <- edgeR::normLibSizes(
    y,
    method = "TMM"
  )
  
  logcpm <- edgeR::cpm(
    y,
    log = TRUE,
    prior.count = 1
  )
  
  return(list(
    logcpm = logcpm,
    norm_factors = y$samples$norm.factors,
    library_sizes = y$samples$lib.size
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_grn_gene_sets
# ------------------------------------------------------------------------------
# Prepare lncRNA and protein-coding gene sets for GRN analysis.
#
# INPUT:
# - expression_matrix: genes x pseudobulk normalized expression matrix
# - gencode: GENCODE annotation data.frame containing SYMBOL and type
#
# OUTPUT:
# - list containing:
#   - lnc: lncRNA genes present in the pseudobulk matrix
#   - prot: protein-coding genes present in the pseudobulk matrix
#   - ambiguous: genes assigned to lncRNA and removed from protein-coding
#   - all: combined non-overlapping GRN gene set
# ------------------------------------------------------------------------------

#' Prepare lncRNA and protein-coding gene sets for GRN analysis
#'
#' Prepare lncRNA and protein-coding gene sets for GRN analysis.
#'
#' @param expression_matrix genes x pseudobulk normalized expression matrix
#' @param gencode GENCODE annotation data.frame containing SYMBOL and type
#'
#' @return list containing: lnc: lncRNA genes present in the pseudobulk matrix prot: protein-coding genes present in the pseudobulk matrix ambiguous: genes assigned to lncRNA and removed from protein-coding all: combined non-overlapping GRN gene set
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
prepare_grn_gene_sets <- function(
    expression_matrix,
    gencode
) {
  
  if (is.null(rownames(expression_matrix))) {
    stop("Expression matrix has no gene names")
  }
  
  if (!all(c("SYMBOL", "type") %in% colnames(gencode))) {
    stop("GENCODE annotation must contain SYMBOL and type columns")
  }
  
  genes <- rownames(expression_matrix)
  
  lnc_raw <- intersect(
    genes,
    unique(gencode$SYMBOL[gencode$type == "lncRNA"])
  )
  
  prot_raw <- intersect(
    genes,
    unique(gencode$SYMBOL[gencode$type == "protein_coding"])
  )
  
  ambiguous_genes <- intersect(
    lnc_raw,
    prot_raw
  )
  
  lnc_features <- unique(lnc_raw)
  
  prot_features <- setdiff(
    unique(prot_raw),
    ambiguous_genes
  )
  
  if (length(intersect(lnc_features, prot_features)) > 0L) {
    stop("lncRNA and protein-coding GRN gene sets are not disjoint")
  }
  
  all_features <- c(
    lnc_features,
    prot_features
  )
  
  return(list(
    lnc = lnc_features,
    prot = prot_features,
    ambiguous = ambiguous_genes,
    all = all_features
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: calculate_grn_expression_statistics
# ------------------------------------------------------------------------------
# Calculate gene-level expression statistics for pseudobulk samples.
#
# INPUT:
# - expression_matrix: normalized pseudobulk expression matrix (logCPM)
# - count_matrix: raw pseudobulk count matrix
#
# OUTPUT:
# - data.frame containing:
#   - gene: gene symbol
#   - mean_logcpm: mean normalized expression across pseudobulk samples
#   - median_logcpm: median normalized expression across pseudobulk samples
#   - n_detected_pseudobulk: number of pseudobulk samples with raw counts > 0
# ------------------------------------------------------------------------------

#' Calculate gene-level expression statistics for pseudobulk samples
#'
#' Calculate gene-level expression statistics for pseudobulk samples.
#'
#' @param expression_matrix normalized pseudobulk expression matrix (logCPM)
#' @param count_matrix raw pseudobulk count matrix
#'
#' @return data.frame containing: gene: gene symbol mean_logcpm: mean normalized expression across pseudobulk samples median_logcpm: median normalized expression across pseudobulk samples n_detected_pseudobulk: number of pseudobulk samples with raw counts > 0
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
calculate_grn_expression_statistics <- function(
    expression_matrix,
    count_matrix
) {
  
  if (!all(dim(expression_matrix) == dim(count_matrix))) {
    stop(
      "Expression and count matrices must have identical dimensions"
    )
  }
  
  if (!identical(
    rownames(expression_matrix),
    rownames(count_matrix)
  )) {
    stop(
      "Gene names do not match between expression and count matrices"
    )
  }
  
  expression_statistics <- data.frame(
    gene = rownames(expression_matrix),
    
    mean_logcpm = rowMeans(
      expression_matrix
    ),
    
    median_logcpm = apply(
      expression_matrix,
      1,
      median
    ),
    
    n_detected_pseudobulk = rowSums(
      count_matrix > 0
    ),
    
    stringsAsFactors = FALSE
  )
  
  return(expression_statistics)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_grn_expression_distribution
# ------------------------------------------------------------------------------
# Visualize gene expression statistics across pseudobulk samples.
#
# INPUT:
# - expression_statistics: data.frame generated by
#   calculate_grn_expression_statistics()
#
# OUTPUT:
# - list containing:
#   - mean_logcpm: histogram of mean logCPM
#   - detected_pseudobulk: histogram of detected pseudobulk samples
# ------------------------------------------------------------------------------

#' Visualize gene expression statistics across pseudobulk samples
#'
#' Visualize gene expression statistics across pseudobulk samples.
#'
#' @param expression_statistics data.frame generated by calculate_grn_expression_statistics()
#'
#' @return list containing: mean_logcpm: histogram of mean logCPM detected_pseudobulk: histogram of detected pseudobulk samples
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_grn_expression_distribution <- function(
    expression_statistics
) {
  
  mean_logcpm_plot <- ggplot2::ggplot(
    expression_statistics,
    ggplot2::aes(x = mean_logcpm)
  ) +
    ggplot2::geom_histogram(
      bins = 50,
      color = "black",
      fill = "#4E79A7"
    ) +
    ggplot2::labs(
      x = "Mean logCPM",
      y = "Number of genes",
      title = "Distribution of Mean Gene Expression"
    ) +
    ggplot2::theme_classic()
  
  detected_pseudobulk_plot <- ggplot2::ggplot(
    expression_statistics,
    ggplot2::aes(x = n_detected_pseudobulk)
  ) +
    ggplot2::geom_histogram(
      binwidth = 1,
      color = "black",
      fill = "#E15759"
    ) +
    ggplot2::labs(
      x = "Detected pseudobulk samples",
      y = "Number of genes",
      title = "Distribution of Gene Detection"
    ) +
    ggplot2::theme_classic()
  
  return(list(
    mean_logcpm = mean_logcpm_plot,
    detected_pseudobulk = detected_pseudobulk_plot
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_grn_correlation_matrix
# ------------------------------------------------------------------------------
# Prepare the normalized expression matrix for GRN correlation analysis.
#
# This function:
# - restricts the normalized pseudobulk matrix to the selected GRN gene set
# - preserves the order of the requested genes
# - identifies genes with zero variance across pseudobulk samples
# - does not remove any genes
#
# INPUT:
# - expression_matrix: normalized pseudobulk expression matrix (logCPM)
# - gene_sets: list generated by prepare_grn_gene_sets()
#
# OUTPUT:
# - list containing:
#   - expression: normalized GRN expression matrix
#   - zero_variance_genes: genes with no expression variability across samples
#   - n_genes: total number of GRN genes
#   - n_samples: total number of pseudobulk samples
# ------------------------------------------------------------------------------

#' Prepare the normalized expression matrix for GRN correlation analysis
#'
#' This function:
#' - restricts the normalized pseudobulk matrix to the selected GRN gene set
#' - preserves the order of the requested genes
#' - identifies genes with zero variance across pseudobulk samples
#' - does not remove any genes
#'
#' @param expression_matrix normalized pseudobulk expression matrix (logCPM)
#' @param gene_sets list generated by prepare_grn_gene_sets()
#'
#' @return list containing: expression: normalized GRN expression matrix zero_variance_genes: genes with no expression variability across samples n_genes: total number of GRN genes n_samples: total number of pseudobulk samples
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
prepare_grn_correlation_matrix <- function(
    expression_matrix,
    gene_sets
) {
  
  required_sets <- c("lnc", "prot", "all")
  
  if (!all(required_sets %in% names(gene_sets))) {
    stop(
      "gene_sets must contain: ",
      paste(required_sets, collapse = ", ")
    )
  }
  
  missing_genes <- setdiff(
    gene_sets$all,
    rownames(expression_matrix)
  )
  
  if (length(missing_genes) > 0L) {
    stop(
      "Some GRN genes are missing from the expression matrix: ",
      length(missing_genes)
    )
  }
  
  grn_expression <- expression_matrix[
    gene_sets$all,
    ,
    drop = FALSE
  ]
  
  gene_variance <- apply(
    grn_expression,
    1,
    stats::var
  )
  
  zero_variance_genes <- names(
    gene_variance[
      is.na(gene_variance) | gene_variance == 0
    ]
  )
  
  return(list(
    expression = grn_expression,
    zero_variance_genes = zero_variance_genes,
    n_genes = nrow(grn_expression),
    n_samples = ncol(grn_expression)
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: compute_grn_spearman_correlations
# ------------------------------------------------------------------------------
# Compute the complete gene-gene Spearman correlation matrix for GRN analysis.
#
# This function:
# - uses normalized pseudobulk expression values
# - computes Spearman correlations across pseudobulk samples
# - preserves all genes and all pairwise correlations
# - does not apply correlation thresholds or network filtering
#
# INPUT:
# - expression_matrix: normalized GRN expression matrix
#   (genes x pseudobulk samples)
#
# OUTPUT:
# - square numeric matrix containing:
#   - rows: genes
#   - columns: genes
#   - values: Spearman correlation coefficients
# ------------------------------------------------------------------------------

#' Compute the complete gene-gene Spearman correlation matrix for GRN analysis
#'
#' This function:
#' - uses normalized pseudobulk expression values
#' - computes Spearman correlations across pseudobulk samples
#' - preserves all genes and all pairwise correlations
#' - does not apply correlation thresholds or network filtering
#'
#' @param expression_matrix normalized GRN expression matrix (genes x pseudobulk samples)
#'
#' @return square numeric matrix containing: rows: genes columns: genes values: Spearman correlation coefficients
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
compute_grn_spearman_correlations <- function(
    expression_matrix
) {
  
  if (is.null(rownames(expression_matrix))) {
    stop("Expression matrix must contain gene names")
  }
  
  if (anyDuplicated(rownames(expression_matrix))) {
    stop("Expression matrix contains duplicated gene names")
  }
  
  if (anyNA(expression_matrix)) {
    stop("Expression matrix contains missing values")
  }
  
  gene_variance <- apply(
    expression_matrix,
    1,
    stats::var
  )
  
  if (any(is.na(gene_variance) | gene_variance == 0)) {
    stop(
      "Expression matrix contains genes with zero or undefined variance"
    )
  }
  
  correlation_matrix <- stats::cor(
    t(expression_matrix),
    method = "spearman",
    use = "everything"
  )
  
  return(correlation_matrix)
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_grn_correlations
# ------------------------------------------------------------------------------
# Summarize the distribution of pairwise Spearman correlations for GRN analysis.
#
# This function:
# - considers only unique gene-gene pairs
# - excludes self-correlations on the matrix diagonal
# - evaluates positive and negative correlations separately
# - reports how many pairs exceed a set of absolute correlation thresholds
# - does not modify or filter the original correlation matrix
#
# INPUT:
# - correlation_matrix: square gene-gene Spearman correlation matrix
# - thresholds: numeric vector of absolute correlation thresholds to evaluate
#
# OUTPUT:
# - list containing:
#   - summary: general statistics for unique pairwise correlations
#   - threshold_counts: number and proportion of gene pairs exceeding each
#     absolute correlation threshold
# ------------------------------------------------------------------------------

#' Summarize the distribution of pairwise Spearman correlations for GRN analysis
#'
#' This function:
#' - considers only unique gene-gene pairs
#' - excludes self-correlations on the matrix diagonal
#' - evaluates positive and negative correlations separately
#' - reports how many pairs exceed a set of absolute correlation thresholds
#' - does not modify or filter the original correlation matrix
#'
#' @param correlation_matrix square gene-gene Spearman correlation matrix
#' @param thresholds numeric vector of absolute correlation thresholds to evaluate
#'
#' @return list containing: summary: general statistics for unique pairwise correlations threshold_counts: number and proportion of gene pairs exceeding each absolute correlation threshold
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
summarize_grn_correlations <- function(
    correlation_matrix,
    thresholds = seq(0.1, 0.9, by = 0.1)
) {
  
  if (nrow(correlation_matrix) != ncol(correlation_matrix)) {
    stop("Correlation matrix must be square")
  }
  
  if (anyNA(correlation_matrix)) {
    stop("Correlation matrix contains missing values")
  }
  
  upper_correlations <- correlation_matrix[
    upper.tri(correlation_matrix, diag = FALSE)
  ]
  
  total_pairs <- length(upper_correlations)
  
  correlation_summary <- data.frame(
    statistic = c(
      "total_unique_pairs",
      "minimum_correlation",
      "maximum_correlation",
      "mean_correlation",
      "median_correlation",
      "positive_pairs",
      "negative_pairs",
      "zero_pairs"
    ),
    value = c(
      total_pairs,
      min(upper_correlations),
      max(upper_correlations),
      mean(upper_correlations),
      median(upper_correlations),
      sum(upper_correlations > 0),
      sum(upper_correlations < 0),
      sum(upper_correlations == 0)
    ),
    stringsAsFactors = FALSE
  )
  
  threshold_counts <- do.call(
    rbind,
    lapply(
      thresholds,
      function(threshold) {
        
        n_absolute <- sum(
          abs(upper_correlations) >= threshold
        )
        
        n_positive <- sum(
          upper_correlations >= threshold
        )
        
        n_negative <- sum(
          upper_correlations <= -threshold
        )
        
        data.frame(
          threshold = threshold,
          n_pairs = n_absolute,
          n_positive = n_positive,
          n_negative = n_negative,
          proportion = n_absolute / total_pairs,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  return(list(
    summary = correlation_summary,
    threshold_counts = threshold_counts
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_grn_correlations_by_gene_type
# ------------------------------------------------------------------------------
# Summarize Spearman correlations according to GRN gene-pair type.
#
# This function:
# - distinguishes lncRNA-lncRNA, lncRNA-protein-coding, and
#   protein-coding-protein-coding gene pairs
# - considers each unique gene-gene pair only once
# - excludes self-correlations
# - evaluates positive and negative correlations separately
# - reports the number and proportion of pairs exceeding each absolute
#   correlation threshold
# - does not modify or filter the original correlation matrix
#
# INPUT:
# - correlation_matrix: square gene-gene Spearman correlation matrix
# - gene_sets: list generated by prepare_grn_gene_sets()
# - thresholds: numeric vector of absolute correlation thresholds to evaluate
#
# OUTPUT:
# - data.frame containing:
#   - pair_type: gene-pair category
#   - threshold: absolute Spearman correlation threshold
#   - total_pairs: total number of possible pairs in the category
#   - n_pairs: number of pairs exceeding the absolute threshold
#   - n_positive: number of positively correlated pairs above the threshold
#   - n_negative: number of negatively correlated pairs below the threshold
#   - proportion: proportion of category-specific pairs exceeding the threshold
# ------------------------------------------------------------------------------

#' Summarize Spearman correlations according to GRN gene-pair type
#'
#' This function:
#' - distinguishes lncRNA-lncRNA, lncRNA-protein-coding, and
#' protein-coding-protein-coding gene pairs
#' - considers each unique gene-gene pair only once
#' - excludes self-correlations
#' - evaluates positive and negative correlations separately
#' - reports the number and proportion of pairs exceeding each absolute
#' correlation threshold
#' - does not modify or filter the original correlation matrix
#'
#' @param correlation_matrix square gene-gene Spearman correlation matrix
#' @param gene_sets list generated by prepare_grn_gene_sets()
#' @param thresholds numeric vector of absolute correlation thresholds to evaluate
#'
#' @return data.frame containing: pair_type: gene-pair category threshold: absolute Spearman correlation threshold total_pairs: total number of possible pairs in the category n_pairs: number of pairs exceeding the absolute threshold n_positive: number of positively correlated pairs above the threshold n_negative: number of negatively correlated pairs below the threshold proportion: proportion of category-specific pairs exceeding the threshold
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
summarize_grn_correlations_by_gene_type <- function(
    correlation_matrix,
    gene_sets,
    thresholds = seq(0.1, 0.9, by = 0.1)
) {
  
  if (nrow(correlation_matrix) != ncol(correlation_matrix)) {
    stop("Correlation matrix must be square")
  }
  
  if (anyNA(correlation_matrix)) {
    stop("Correlation matrix contains missing values")
  }
  
  if (is.null(rownames(correlation_matrix))) {
    stop("Correlation matrix must contain gene names")
  }
  
  lnc_genes <- intersect(
    gene_sets$lnc,
    rownames(correlation_matrix)
  )
  
  prot_genes <- intersect(
    gene_sets$prot,
    rownames(correlation_matrix)
  )
  
  lnc_lnc <- correlation_matrix[
    lnc_genes,
    lnc_genes,
    drop = FALSE
  ]
  
  lnc_lnc <- lnc_lnc[
    upper.tri(lnc_lnc, diag = FALSE)
  ]
  
  lnc_prot <- as.vector(
    correlation_matrix[
      lnc_genes,
      prot_genes,
      drop = FALSE
    ]
  )
  
  prot_prot <- correlation_matrix[
    prot_genes,
    prot_genes,
    drop = FALSE
  ]
  
  prot_prot <- prot_prot[
    upper.tri(prot_prot, diag = FALSE)
  ]
  
  correlation_groups <- list(
    lncRNA_lncRNA = lnc_lnc,
    lncRNA_protein_coding = lnc_prot,
    protein_coding_protein_coding = prot_prot
  )
  
  results <- do.call(
    rbind,
    lapply(
      names(correlation_groups),
      function(pair_type) {
        
        correlations <- correlation_groups[[pair_type]]
        total_pairs <- length(correlations)
        
        do.call(
          rbind,
          lapply(
            thresholds,
            function(threshold) {
              
              n_pairs <- sum(
                abs(correlations) >= threshold
              )
              
              n_positive <- sum(
                correlations >= threshold
              )
              
              n_negative <- sum(
                correlations <= -threshold
              )
              
              data.frame(
                pair_type = pair_type,
                threshold = threshold,
                total_pairs = total_pairs,
                n_pairs = n_pairs,
                n_positive = n_positive,
                n_negative = n_negative,
                proportion = n_pairs / total_pairs,
                stringsAsFactors = FALSE
              )
            }
          )
        )
      }
    )
  )
  
  rownames(results) <- NULL
  
  return(results)
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_grn_correlations_by_detection
# ------------------------------------------------------------------------------
# Summarize strong GRN correlations according to gene detection frequency.
#
# This function:
# - evaluates correlations at selected absolute Spearman thresholds
# - considers each unique gene-gene pair only once
# - associates each gene with its pseudobulk detection frequency
# - summarizes the minimum detection frequency within each correlated gene pair
# - distinguishes pairs involving at least one lncRNA
# - does not modify or filter the original correlation matrix
#
# INPUT:
# - correlation_matrix: square gene-gene Spearman correlation matrix
# - expression_statistics: data.frame generated by
#   calculate_grn_expression_statistics()
# - gene_sets: list generated by prepare_grn_gene_sets()
# - thresholds: numeric vector of absolute correlation thresholds to evaluate
#
# OUTPUT:
# - data.frame containing:
#   - threshold: absolute Spearman correlation threshold
#   - pair_group: all pairs or pairs involving at least one lncRNA
#   - n_pairs: number of pairs exceeding the threshold
#   - median_min_detection: median of the minimum detection frequency per pair
#   - q1_min_detection: first quartile of minimum detection frequency per pair
#   - q3_min_detection: third quartile of minimum detection frequency per pair
#   - pairs_min_detection_le_1: pairs with at least one gene detected in <= 1
#     pseudobulk sample
#   - pairs_min_detection_le_5: pairs with at least one gene detected in <= 5
#     pseudobulk samples
# ------------------------------------------------------------------------------

#' Summarize strong GRN correlations according to gene detection frequency
#'
#' This function:
#' - evaluates correlations at selected absolute Spearman thresholds
#' - considers each unique gene-gene pair only once
#' - associates each gene with its pseudobulk detection frequency
#' - summarizes the minimum detection frequency within each correlated gene pair
#' - distinguishes pairs involving at least one lncRNA
#' - does not modify or filter the original correlation matrix
#'
#' @param correlation_matrix square gene-gene Spearman correlation matrix
#' @param expression_statistics data.frame generated by calculate_grn_expression_statistics()
#' @param gene_sets list generated by prepare_grn_gene_sets()
#' @param thresholds numeric vector of absolute correlation thresholds to evaluate
#'
#' @return data.frame containing: threshold: absolute Spearman correlation threshold pair_group: all pairs or pairs involving at least one lncRNA n_pairs: number of pairs exceeding the threshold median_min_detection: median of the minimum detection frequency per pair q1_min_detection: first quartile of minimum detection frequency per pair q3_min_detection: third quartile of minimum detection frequency per pair pairs_min_detection_le_1: pairs with at least one gene detected in <= 1 pseudobulk sample pairs_min_detection_le_5: pairs with at least one gene detected in <= 5 pseudobulk samples
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
summarize_grn_correlations_by_detection <- function(
    correlation_matrix,
    expression_statistics,
    gene_sets,
    thresholds = c(0.7, 0.8, 0.9)
) {
  
  genes <- rownames(correlation_matrix)
  
  detection <- expression_statistics$n_detected_pseudobulk
  names(detection) <- expression_statistics$gene
  
  if (!all(genes %in% names(detection))) {
    stop(
      "Detection statistics are missing for some genes in the correlation matrix"
    )
  }
  
  detection <- detection[genes]
  
  lnc_status <- genes %in% gene_sets$lnc
  
  results <- list()
  result_index <- 1L
  
  for (threshold in thresholds) {
    
    pair_indices <- which(
      abs(correlation_matrix) >= threshold &
        upper.tri(correlation_matrix, diag = FALSE),
      arr.ind = TRUE
    )
    
    pair_min_detection <- pmin(
      detection[pair_indices[, 1]],
      detection[pair_indices[, 2]]
    )
    
    pair_has_lnc <- (
      lnc_status[pair_indices[, 1]] |
        lnc_status[pair_indices[, 2]]
    )
    
    groups <- list(
      all_pairs = rep(TRUE, length(pair_min_detection)),
      lncRNA_involving = pair_has_lnc
    )
    
    for (group_name in names(groups)) {
      
      keep <- groups[[group_name]]
      current_detection <- pair_min_detection[keep]
      
      results[[result_index]] <- data.frame(
        threshold = threshold,
        pair_group = group_name,
        n_pairs = length(current_detection),
        median_min_detection = median(current_detection),
        q1_min_detection = as.numeric(
          quantile(current_detection, 0.25)
        ),
        q3_min_detection = as.numeric(
          quantile(current_detection, 0.75)
        ),
        pairs_min_detection_le_1 = sum(
          current_detection <= 1
        ),
        pairs_min_detection_le_5 = sum(
          current_detection <= 5
        ),
        stringsAsFactors = FALSE
      )
      
      result_index <- result_index + 1L
    }
  }
  
  results <- do.call(
    rbind,
    results
  )
  
  rownames(results) <- NULL
  
  return(results)
}

# ------------------------------------------------------------------------------
# FUNCTION: evaluate_grn_detection_thresholds
# ------------------------------------------------------------------------------
# Evaluate candidate pseudobulk detection thresholds for GRN expression
# filtering.
#
# This function:
# - evaluates multiple minimum pseudobulk detection thresholds
# - counts the number of retained lncRNA and protein-coding genes
# - calculates retention proportions relative to the original GRN gene sets
# - supports threshold selection before expression filtering
# - does not modify or filter the expression matrix
#
# INPUT:
# - expression_statistics: data.frame generated by
#   calculate_grn_expression_statistics()
# - gene_sets: list generated by prepare_grn_gene_sets()
# - thresholds: integer vector defining the minimum number of pseudobulk samples
#   in which a gene must be detected
#
# OUTPUT:
# - data.frame containing:
#   - threshold: minimum number of pseudobulk samples required for detection
#   - retained_lnc: number of retained lncRNA genes
#   - retained_prot: number of retained protein-coding genes
#   - retained_total: total number of retained GRN genes
#   - lnc_retention: proportion of original lncRNA genes retained
#   - prot_retention: proportion of original protein-coding genes retained
#   - total_retention: proportion of all original GRN genes retained
# ------------------------------------------------------------------------------

#' Evaluate candidate pseudobulk detection thresholds for GRN expression
#'
#' filtering.
#' This function:
#' - evaluates multiple minimum pseudobulk detection thresholds
#' - counts the number of retained lncRNA and protein-coding genes
#' - calculates retention proportions relative to the original GRN gene sets
#' - supports threshold selection before expression filtering
#' - does not modify or filter the expression matrix
#'
#' @param expression_statistics data.frame generated by calculate_grn_expression_statistics()
#' @param gene_sets list generated by prepare_grn_gene_sets()
#' @param thresholds integer vector defining the minimum number of pseudobulk samples in which a gene must be detected
#'
#' @return data.frame containing: threshold: minimum number of pseudobulk samples required for detection retained_lnc: number of retained lncRNA genes retained_prot: number of retained protein-coding genes retained_total: total number of retained GRN genes lnc_retention: proportion of original lncRNA genes retained prot_retention: proportion of original protein-coding genes retained total_retention: proportion of all original GRN genes retained
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
evaluate_grn_detection_thresholds <- function(
    expression_statistics,
    gene_sets,
    thresholds = c(1, 2, 3, 5, 10, 15)
) {
  
  required_columns <- c(
    "gene",
    "n_detected_pseudobulk"
  )
  
  if (!all(required_columns %in% colnames(expression_statistics))) {
    stop(
      "expression_statistics must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (!all(c("lnc", "prot", "all") %in% names(gene_sets))) {
    stop(
      "gene_sets must contain lnc, prot, and all gene sets"
    )
  }
  
  if (
    anyNA(thresholds) ||
    any(thresholds < 1) ||
    any(thresholds %% 1 != 0)
  ) {
    stop(
      "thresholds must contain positive integer values"
    )
  }
  
  results <- lapply(
    thresholds,
    function(threshold) {
      
      retained_genes <- expression_statistics$gene[
        expression_statistics$n_detected_pseudobulk >= threshold
      ]
      
      retained_lnc <- intersect(
        retained_genes,
        gene_sets$lnc
      )
      
      retained_prot <- intersect(
        retained_genes,
        gene_sets$prot
      )
      
      retained_total <- length(retained_lnc) +
        length(retained_prot)
      
      data.frame(
        threshold = threshold,
        retained_lnc = length(retained_lnc),
        retained_prot = length(retained_prot),
        retained_total = retained_total,
        lnc_retention =
          length(retained_lnc) / length(gene_sets$lnc),
        prot_retention =
          length(retained_prot) / length(gene_sets$prot),
        total_retention =
          retained_total / length(gene_sets$all),
        stringsAsFactors = FALSE
      )
    }
  )
  
  results <- do.call(
    rbind,
    results
  )
  
  rownames(results) <- NULL
  
  return(results)
}

# ------------------------------------------------------------------------------
# FUNCTION: apply_grn_expression_filter
# ------------------------------------------------------------------------------
# Apply pseudobulk detection-based expression filtering for GRN analysis.
#
# This function:
# - retains genes detected in at least the selected number of pseudobulk samples
# - preserves lncRNA and protein-coding gene annotations
# - filters both normalized expression and raw pseudobulk count matrices
# - returns filtered gene sets for downstream GRN inference
#
# INPUT:
# - expression_matrix: normalized pseudobulk expression matrix (logCPM)
# - count_matrix: raw pseudobulk count matrix
# - expression_statistics: data.frame generated by
#   calculate_grn_expression_statistics()
# - gene_sets: list generated by prepare_grn_gene_sets()
# - min_detected_pseudobulk: minimum number of pseudobulk samples in which
#   a gene must be detected
#
# OUTPUT:
# - list containing:
#   - expression: filtered normalized pseudobulk expression matrix
#   - counts: filtered raw pseudobulk count matrix
#   - gene_sets: filtered GRN gene sets
#   - retained_genes: genes retained after filtering
#   - removed_genes: genes removed by filtering
# ------------------------------------------------------------------------------

#' Apply pseudobulk detection-based expression filtering for GRN analysis
#'
#' This function:
#' - retains genes detected in at least the selected number of pseudobulk samples
#' - preserves lncRNA and protein-coding gene annotations
#' - filters both normalized expression and raw pseudobulk count matrices
#' - returns filtered gene sets for downstream GRN inference
#'
#' @param expression_matrix normalized pseudobulk expression matrix (logCPM)
#' @param count_matrix raw pseudobulk count matrix
#' @param expression_statistics data.frame generated by calculate_grn_expression_statistics()
#' @param gene_sets list generated by prepare_grn_gene_sets()
#' @param min_detected_pseudobulk minimum number of pseudobulk samples in which a gene must be detected
#'
#' @return list containing: expression: filtered normalized pseudobulk expression matrix counts: filtered raw pseudobulk count matrix gene_sets: filtered GRN gene sets retained_genes: genes retained after filtering removed_genes: genes removed by filtering
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
apply_grn_expression_filter <- function(
    expression_matrix,
    count_matrix,
    expression_statistics,
    gene_sets,
    min_detected_pseudobulk
) {
  
  if (!identical(
    rownames(expression_matrix),
    rownames(count_matrix)
  )) {
    stop(
      "Gene names do not match between expression and count matrices"
    )
  }
  
  if (!all(c(
    "gene",
    "n_detected_pseudobulk"
  ) %in% colnames(expression_statistics))) {
    stop(
      "expression_statistics must contain gene and n_detected_pseudobulk"
    )
  }
  
  retained_genes <- expression_statistics$gene[
    expression_statistics$n_detected_pseudobulk >=
      min_detected_pseudobulk
  ]
  
  retained_genes <- intersect(
    gene_sets$all,
    retained_genes
  )
  
  removed_genes <- setdiff(
    gene_sets$all,
    retained_genes
  )
  
  filtered_expression <- expression_matrix[
    retained_genes,
    ,
    drop = FALSE
  ]
  
  filtered_counts <- count_matrix[
    retained_genes,
    ,
    drop = FALSE
  ]
  
  filtered_lnc <- intersect(
    gene_sets$lnc,
    retained_genes
  )
  
  filtered_prot <- intersect(
    gene_sets$prot,
    retained_genes
  )
  
  filtered_gene_sets <- list(
    lnc = filtered_lnc,
    prot = filtered_prot,
    all = c(
      filtered_lnc,
      filtered_prot
    )
  )
  
  return(list(
    expression = filtered_expression,
    counts = filtered_counts,
    gene_sets = filtered_gene_sets,
    retained_genes = retained_genes,
    removed_genes = removed_genes
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_genie3_regulators
# ------------------------------------------------------------------------------
# Prepare candidate transcription factor regulators for GENIE3 inference.
#
# This function:
# - selects genes classified as transcription factors in the Human TF database
# - extracts HGNC gene symbols for confirmed TFs
# - intersects TFs with filtered protein-coding GRN genes
# - ensures that candidate regulators are present in the filtered expression
#   matrix
#
# INPUT:
# - human_tf: Human TF database data.frame
# - filtered_gene_sets: filtered GRN gene sets generated by
#   apply_grn_expression_filter()
# - expression_matrix: filtered normalized GRN expression matrix
#
# OUTPUT:
# - character vector containing candidate TF regulators for GENIE3
# ------------------------------------------------------------------------------

#' Prepare candidate transcription factor regulators for GENIE3 inference
#'
#' This function:
#' - selects genes classified as transcription factors in the Human TF database
#' - extracts HGNC gene symbols for confirmed TFs
#' - intersects TFs with filtered protein-coding GRN genes
#' - ensures that candidate regulators are present in the filtered expression
#' matrix
#'
#' @param human_tf Human TF database data.frame
#' @param filtered_gene_sets filtered GRN gene sets generated by apply_grn_expression_filter()
#' @param expression_matrix filtered normalized GRN expression matrix
#'
#' @return character vector containing candidate TF regulators for GENIE3
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
prepare_genie3_regulators <- function(
    human_tf,
    filtered_gene_sets,
    expression_matrix
) {
  
  required_columns <- c(
    "HGNC.symbol",
    "Is.TF."
  )
  
  if (!all(required_columns %in% colnames(human_tf))) {
    stop(
      "Human TF database must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  human_tf_genes <- unique(
    human_tf$HGNC.symbol[
      human_tf$Is.TF. == "Yes"
    ]
  )
  
  candidate_regulators <- intersect(
    human_tf_genes,
    filtered_gene_sets$prot
  )
  
  candidate_regulators <- intersect(
    candidate_regulators,
    rownames(expression_matrix)
  )
  
  if (length(candidate_regulators) == 0L) {
    stop("No candidate TF regulators matched the filtered GRN data")
  }
  
  return(candidate_regulators)
}

# ------------------------------------------------------------------------------
# FUNCTION: integrate_genie3_spearman
# ------------------------------------------------------------------------------
# Integrate GENIE3 regulatory importance scores with Spearman correlations.
#
# This function:
# - extracts TF-target relationships evaluated by GENIE3
# - retrieves the corresponding Spearman correlation for each TF-target pair
# - preserves GENIE3 directionality from regulator to target
# - removes self-regulatory interactions
# - returns an integrated edge table for downstream GRN construction
#
# INPUT:
# - genie3_weights: GENIE3 importance matrix with regulators as rows and
#   target genes as columns
# - spearman_matrix: gene-by-gene Spearman correlation matrix
#
# OUTPUT:
# - data.frame containing:
#   regulator: candidate transcription factor
#   target: target gene
#   genie3_weight: GENIE3 importance score
#   spearman_rho: corresponding Spearman correlation
# ------------------------------------------------------------------------------

#' Integrate GENIE3 regulatory importance scores with Spearman correlations
#'
#' This function:
#' - extracts TF-target relationships evaluated by GENIE3
#' - retrieves the corresponding Spearman correlation for each TF-target pair
#' - preserves GENIE3 directionality from regulator to target
#' - removes self-regulatory interactions
#' - returns an integrated edge table for downstream GRN construction
#'
#' @param genie3_weights GENIE3 importance matrix with regulators as rows and target genes as columns
#' @param spearman_matrix gene-by-gene Spearman correlation matrix
#'
#' @return data.frame containing: regulator: candidate transcription factor target: target gene genie3_weight: GENIE3 importance score spearman_rho: corresponding Spearman correlation
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
integrate_genie3_spearman <- function(
    genie3_weights,
    spearman_matrix
) {
  
  regulators <- rownames(genie3_weights)
  targets <- colnames(genie3_weights)
  
  if (is.null(regulators) || is.null(targets)) {
    stop("GENIE3 weight matrix must contain row and column names")
  }
  
  if (!all(regulators %in% rownames(spearman_matrix)) ||
      !all(targets %in% colnames(spearman_matrix))) {
    stop("GENIE3 genes are not fully represented in the Spearman matrix")
  }
  
  spearman_subset <- spearman_matrix[
    regulators,
    targets,
    drop = FALSE
  ]
  
  regulator_index <- rep(
    seq_along(regulators),
    times = length(targets)
  )
  
  target_index <- rep(
    seq_along(targets),
    each = length(regulators)
  )
  
  keep <- regulators[regulator_index] != targets[target_index]
  
  integrated_edges <- data.frame(
    regulator = regulators[regulator_index[keep]],
    target = targets[target_index[keep]],
    genie3_weight = as.vector(genie3_weights)[keep],
    spearman_rho = as.vector(spearman_subset)[keep],
    stringsAsFactors = FALSE
  )
  
  return(integrated_edges)
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_grn_integrated_edges
# ------------------------------------------------------------------------------
# Summarize GENIE3 importance scores and Spearman correlations for integrated
# candidate regulatory interactions.
#
# This function:
# - summarizes the distribution of GENIE3 importance scores
# - summarizes the distribution of absolute Spearman correlations
# - calculates selected quantiles for both measures
# - evaluates the number of candidate interactions retained at selected
#   GENIE3 and Spearman thresholds
#
# No interactions are removed at this stage.
#
# INPUT:
# - integrated_edges: data.frame generated by integrate_genie3_spearman()
# - genie3_quantiles: quantiles used to summarize GENIE3 importance scores
# - spearman_thresholds: absolute Spearman correlation thresholds to evaluate
#
# OUTPUT:
# - list containing:
#   genie3_summary: summary statistics for GENIE3 importance scores
#   spearman_summary: summary statistics for absolute Spearman correlations
#   genie3_quantiles: GENIE3 importance values at selected quantiles
#   threshold_counts: number of interactions retained for combinations of
#     GENIE3 quantiles and absolute Spearman thresholds
# ------------------------------------------------------------------------------

#' Summarize GENIE3 importance scores and Spearman correlations for integrated
#'
#' candidate regulatory interactions.
#' This function:
#' - summarizes the distribution of GENIE3 importance scores
#' - summarizes the distribution of absolute Spearman correlations
#' - calculates selected quantiles for both measures
#' - evaluates the number of candidate interactions retained at selected
#' GENIE3 and Spearman thresholds
#' No interactions are removed at this stage.
#'
#' @param integrated_edges data.frame generated by integrate_genie3_spearman()
#' @param genie3_quantiles quantiles used to summarize GENIE3 importance scores
#' @param spearman_thresholds absolute Spearman correlation thresholds to evaluate
#'
#' @return list containing: genie3_summary: summary statistics for GENIE3 importance scores spearman_summary: summary statistics for absolute Spearman correlations genie3_quantiles: GENIE3 importance values at selected quantiles threshold_counts: number of interactions retained for combinations of GENIE3 quantiles and absolute Spearman thresholds
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
summarize_grn_integrated_edges <- function(
    integrated_edges,
    genie3_quantiles = c(0.90, 0.95, 0.99, 0.995, 0.999),
    spearman_thresholds = c(0.3, 0.4, 0.5, 0.6, 0.7)
) {
  
  genie3_weight <- integrated_edges$genie3_weight
  abs_spearman <- abs(integrated_edges$spearman_rho)
  
  genie3_summary <- summary(genie3_weight)
  spearman_summary <- summary(abs_spearman)
  
  genie3_cutoffs <- quantile(
    genie3_weight,
    probs = genie3_quantiles,
    na.rm = TRUE,
    names = FALSE
  )
  
  genie3_quantile_table <- data.frame(
    quantile = genie3_quantiles,
    genie3_cutoff = genie3_cutoffs
  )
  
  threshold_counts <- do.call(
    rbind,
    lapply(seq_along(genie3_quantiles), function(i) {
      
      do.call(
        rbind,
        lapply(spearman_thresholds, function(rho_threshold) {
          
          keep <- genie3_weight >= genie3_cutoffs[i] &
            abs_spearman >= rho_threshold
          
          data.frame(
            genie3_quantile = genie3_quantiles[i],
            genie3_cutoff = genie3_cutoffs[i],
            spearman_threshold = rho_threshold,
            n_edges = sum(keep),
            proportion = mean(keep)
          )
        })
      )
    })
  )
  
  return(
    list(
      genie3_summary = genie3_summary,
      spearman_summary = spearman_summary,
      genie3_quantiles = genie3_quantile_table,
      threshold_counts = threshold_counts
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: evaluate_grn_candidate_networks
# ------------------------------------------------------------------------------
# Evaluate candidate GRN compositions across GENIE3 and Spearman thresholds.
#
# This function:
# - applies candidate GENIE3 quantile cutoffs and absolute Spearman thresholds
# - counts retained regulatory interactions
# - counts unique regulators and targets
# - quantifies lncRNA and protein-coding target representation
# - reports the proportion of filtered lncRNA genes represented as targets
# - does not define or return the final GRN
#
# INPUT:
# - integrated_edges: data.frame generated by integrate_genie3_spearman()
# - filtered_gene_sets: filtered GRN gene sets generated by
#   apply_grn_expression_filter()
# - genie3_quantile_table: data.frame generated by
#   summarize_grn_integrated_edges()$genie3_quantiles
# - spearman_thresholds: absolute Spearman thresholds to evaluate
#
# OUTPUT:
# - data.frame containing:
#   - genie3_quantile
#   - genie3_cutoff
#   - spearman_threshold
#   - n_edges
#   - n_regulators
#   - n_targets
#   - n_lnc_targets
#   - n_prot_targets
#   - lnc_target_retention
# ------------------------------------------------------------------------------

#' Evaluate candidate GRN compositions across GENIE3 and Spearman thresholds
#'
#' This function:
#' - applies candidate GENIE3 quantile cutoffs and absolute Spearman thresholds
#' - counts retained regulatory interactions
#' - counts unique regulators and targets
#' - quantifies lncRNA and protein-coding target representation
#' - reports the proportion of filtered lncRNA genes represented as targets
#' - does not define or return the final GRN
#'
#' @param integrated_edges data.frame generated by integrate_genie3_spearman()
#' @param filtered_gene_sets filtered GRN gene sets generated by apply_grn_expression_filter()
#' @param genie3_quantile_table data.frame generated by summarize_grn_integrated_edges()$genie3_quantiles
#' @param spearman_thresholds absolute Spearman thresholds to evaluate
#'
#' @return data.frame containing: genie3_quantile genie3_cutoff spearman_threshold n_edges n_regulators n_targets n_lnc_targets n_prot_targets lnc_target_retention
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
evaluate_grn_candidate_networks <- function(
    integrated_edges,
    filtered_gene_sets,
    genie3_quantile_table,
    spearman_thresholds = c(0.3, 0.4, 0.5, 0.6, 0.7)
) {
  
  results <- list()
  result_index <- 1L
  
  for (i in seq_len(nrow(genie3_quantile_table))) {
    
    genie3_quantile <- genie3_quantile_table$quantile[i]
    genie3_cutoff <- genie3_quantile_table$genie3_cutoff[i]
    
    for (rho_threshold in spearman_thresholds) {
      
      keep <- integrated_edges$genie3_weight >= genie3_cutoff &
        abs(integrated_edges$spearman_rho) >= rho_threshold
      
      candidate_edges <- integrated_edges[keep, , drop = FALSE]
      
      candidate_targets <- unique(candidate_edges$target)
      
      lnc_targets <- intersect(
        candidate_targets,
        filtered_gene_sets$lnc
      )
      
      prot_targets <- intersect(
        candidate_targets,
        filtered_gene_sets$prot
      )
      
      results[[result_index]] <- data.frame(
        genie3_quantile = genie3_quantile,
        genie3_cutoff = genie3_cutoff,
        spearman_threshold = rho_threshold,
        n_edges = nrow(candidate_edges),
        n_regulators = length(unique(candidate_edges$regulator)),
        n_targets = length(candidate_targets),
        n_lnc_targets = length(lnc_targets),
        n_prot_targets = length(prot_targets),
        lnc_target_retention =
          length(lnc_targets) / length(filtered_gene_sets$lnc),
        stringsAsFactors = FALSE
      )
      
      result_index <- result_index + 1L
    }
  }
  
  results <- do.call(rbind, results)
  
  rownames(results) <- NULL
  
  return(results)
}

# ------------------------------------------------------------------------------
# FUNCTION: evaluate_grn_lnc_edge_representation
# ------------------------------------------------------------------------------
# Evaluate lncRNA edge representation across selected candidate GRN thresholds.
#
# This function:
# - applies selected GENIE3 and Spearman thresholds
# - counts TF-to-lncRNA and TF-to-protein-coding regulatory interactions
# - counts unique lncRNA and protein-coding targets
# - reports the proportion of retained edges targeting lncRNA genes
# - does not construct the final GRN
#
# INPUT:
# - integrated_edges: data.frame generated by integrate_genie3_spearman()
# - filtered_gene_sets: filtered GRN gene sets generated by
#   apply_grn_expression_filter()
# - candidate_settings: data.frame containing selected combinations of:
#   - genie3_quantile
#   - genie3_cutoff
#   - spearman_threshold
#
# OUTPUT:
# - data.frame containing:
#   - genie3_quantile
#   - genie3_cutoff
#   - spearman_threshold
#   - n_edges
#   - n_lnc_edges
#   - n_prot_edges
#   - lnc_edge_proportion
#   - n_lnc_targets
#   - n_prot_targets
# ------------------------------------------------------------------------------

#' Evaluate lncRNA edge representation across selected candidate GRN thresholds
#'
#' This function:
#' - applies selected GENIE3 and Spearman thresholds
#' - counts TF-to-lncRNA and TF-to-protein-coding regulatory interactions
#' - counts unique lncRNA and protein-coding targets
#' - reports the proportion of retained edges targeting lncRNA genes
#' - does not construct the final GRN
#'
#' @param integrated_edges data.frame generated by integrate_genie3_spearman()
#' @param filtered_gene_sets filtered GRN gene sets generated by apply_grn_expression_filter()
#' @param candidate_settings data.frame containing selected combinations of: - genie3_quantile - genie3_cutoff - spearman_threshold
#'
#' @return data.frame containing: genie3_quantile genie3_cutoff spearman_threshold n_edges n_lnc_edges n_prot_edges lnc_edge_proportion n_lnc_targets n_prot_targets
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
evaluate_grn_lnc_edge_representation <- function(
    integrated_edges,
    filtered_gene_sets,
    candidate_settings
) {
  
  results <- lapply(
    seq_len(nrow(candidate_settings)),
    function(i) {
      
      genie3_quantile <- candidate_settings$genie3_quantile[i]
      genie3_cutoff <- candidate_settings$genie3_cutoff[i]
      rho_threshold <- candidate_settings$spearman_threshold[i]
      
      keep <- integrated_edges$genie3_weight >= genie3_cutoff &
        abs(integrated_edges$spearman_rho) >= rho_threshold
      
      candidate_edges <- integrated_edges[keep, , drop = FALSE]
      
      lnc_status <- candidate_edges$target %in% filtered_gene_sets$lnc
      prot_status <- candidate_edges$target %in% filtered_gene_sets$prot
      
      data.frame(
        genie3_quantile = genie3_quantile,
        genie3_cutoff = genie3_cutoff,
        spearman_threshold = rho_threshold,
        n_edges = nrow(candidate_edges),
        n_lnc_edges = sum(lnc_status),
        n_prot_edges = sum(prot_status),
        lnc_edge_proportion = sum(lnc_status) / nrow(candidate_edges),
        n_lnc_targets = length(unique(candidate_edges$target[lnc_status])),
        n_prot_targets = length(unique(candidate_edges$target[prot_status])),
        stringsAsFactors = FALSE
      )
    }
  )
  
  results <- do.call(rbind, results)
  
  rownames(results) <- NULL
  
  return(results)
}

# ------------------------------------------------------------------------------
# FUNCTION: construct_final_grn
# ------------------------------------------------------------------------------
# Construct the final GRN edge list using selected GENIE3 and Spearman criteria.
#
# This function:
# - filters integrated TF-target interactions using a GENIE3 importance cutoff
# - applies an absolute Spearman correlation threshold
# - preserves GENIE3 importance and Spearman correlation values
# - classifies the sign of the Spearman association
# - returns the final GRN edge list for downstream network analysis
#
# INPUT:
# - integrated_edges: data.frame generated by integrate_genie3_spearman()
# - genie3_cutoff: minimum GENIE3 importance score
# - spearman_threshold: minimum absolute Spearman correlation
#
# OUTPUT:
# - data.frame containing:
#   - regulator
#   - target
#   - genie3_weight
#   - spearman_rho
#   - correlation_sign
# ------------------------------------------------------------------------------

#' Construct the final GRN edge list using selected GENIE3 and Spearman criteria
#'
#' This function:
#' - filters integrated TF-target interactions using a GENIE3 importance cutoff
#' - applies an absolute Spearman correlation threshold
#' - preserves GENIE3 importance and Spearman correlation values
#' - classifies the sign of the Spearman association
#' - returns the final GRN edge list for downstream network analysis
#'
#' @param integrated_edges data.frame generated by integrate_genie3_spearman()
#' @param genie3_cutoff minimum GENIE3 importance score
#' @param spearman_threshold minimum absolute Spearman correlation
#'
#' @return data.frame containing: regulator target genie3_weight spearman_rho correlation_sign
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
construct_final_grn <- function(
    integrated_edges,
    genie3_cutoff,
    spearman_threshold
) {
  
  keep <- integrated_edges$genie3_weight >= genie3_cutoff &
    abs(integrated_edges$spearman_rho) >= spearman_threshold
  
  final_edges <- integrated_edges[
    keep,
    ,
    drop = FALSE
  ]
  
  final_edges$correlation_sign <- ifelse(
    final_edges$spearman_rho >= 0,
    "positive",
    "negative"
  )
  
  rownames(final_edges) <- NULL
  
  return(final_edges)
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_grn_node_table
# ------------------------------------------------------------------------------
# Prepare an annotated node table for the final GRN.
#
# This function:
# - identifies all unique genes represented in the final GRN
# - classifies genes as lncRNA or protein-coding
# - identifies genes acting as transcription factor regulators
# - records whether each gene appears as a regulator and/or target
# - prepares node annotations for downstream hub analysis and visualization
#
# INPUT:
# - grn_edges: final GRN edge list generated by construct_final_grn()
# - gene_sets: filtered GRN gene sets generated by
#   apply_grn_expression_filter()
# - tf_regulators: candidate transcription factors used for GENIE3 inference
#
# OUTPUT:
# - data.frame containing:
#   - gene
#   - gene_type
#   - is_tf
#   - is_regulator
#   - is_target
# ------------------------------------------------------------------------------

#' Prepare an annotated node table for the final GRN
#'
#' This function:
#' - identifies all unique genes represented in the final GRN
#' - classifies genes as lncRNA or protein-coding
#' - identifies genes acting as transcription factor regulators
#' - records whether each gene appears as a regulator and/or target
#' - prepares node annotations for downstream hub analysis and visualization
#'
#' @param grn_edges final GRN edge list generated by construct_final_grn()
#' @param gene_sets filtered GRN gene sets generated by apply_grn_expression_filter()
#' @param tf_regulators candidate transcription factors used for GENIE3 inference
#'
#' @return data.frame containing: gene gene_type is_tf is_regulator is_target
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
prepare_grn_node_table <- function(
    grn_edges,
    gene_sets,
    tf_regulators
) {
  
  genes <- unique(c(
    grn_edges$regulator,
    grn_edges$target
  ))
  
  gene_type <- ifelse(
    genes %in% gene_sets$lnc,
    "lncRNA",
    ifelse(
      genes %in% gene_sets$prot,
      "protein_coding",
      "unknown"
    )
  )
  
  nodes <- data.frame(
    gene = genes,
    gene_type = gene_type,
    is_tf = genes %in% tf_regulators,
    is_regulator = genes %in% grn_edges$regulator,
    is_target = genes %in% grn_edges$target,
    stringsAsFactors = FALSE
  )
  
  return(nodes)
}

# ------------------------------------------------------------------------------
# FUNCTION: calculate_grn_node_connectivity
# ------------------------------------------------------------------------------
# Calculate node-level connectivity metrics for the final GRN.
#
# This function:
# - calculates out-degree as the number of outgoing regulatory interactions
# - calculates in-degree as the number of incoming regulatory interactions
# - calculates total degree as the sum of in-degree and out-degree
# - integrates connectivity metrics with GRN node annotations
#
# INPUT:
# - grn_edges: final GRN edge list generated by construct_final_grn()
# - grn_nodes: annotated node table generated by prepare_grn_node_table()
#
# OUTPUT:
# - data.frame containing:
#   - gene
#   - gene_type
#   - is_tf
#   - is_regulator
#   - is_target
#   - out_degree
#   - in_degree
#   - total_degree
# ------------------------------------------------------------------------------

#' Calculate node-level connectivity metrics for the final GRN
#'
#' This function:
#' - calculates out-degree as the number of outgoing regulatory interactions
#' - calculates in-degree as the number of incoming regulatory interactions
#' - calculates total degree as the sum of in-degree and out-degree
#' - integrates connectivity metrics with GRN node annotations
#'
#' @param grn_edges final GRN edge list generated by construct_final_grn()
#' @param grn_nodes annotated node table generated by prepare_grn_node_table()
#'
#' @return data.frame containing: gene gene_type is_tf is_regulator is_target out_degree in_degree total_degree
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
calculate_grn_node_connectivity <- function(
    grn_edges,
    grn_nodes
) {
  
  out_degree <- table(grn_edges$regulator)
  in_degree <- table(grn_edges$target)
  
  grn_connectivity <- grn_nodes
  
  grn_connectivity$out_degree <- as.integer(
    out_degree[grn_connectivity$gene]
  )
  
  grn_connectivity$in_degree <- as.integer(
    in_degree[grn_connectivity$gene]
  )
  
  grn_connectivity$out_degree[
    is.na(grn_connectivity$out_degree)
  ] <- 0L
  
  grn_connectivity$in_degree[
    is.na(grn_connectivity$in_degree)
  ] <- 0L
  
  grn_connectivity$total_degree <-
    grn_connectivity$out_degree +
    grn_connectivity$in_degree
  
  return(grn_connectivity)
}

# ------------------------------------------------------------------------------
# FUNCTION: evaluate_grn_hub_thresholds
# ------------------------------------------------------------------------------
# Evaluate candidate degree thresholds for GRN hub identification.
#
# This function:
# - evaluates high quantiles of out-degree, in-degree, and total degree
# - calculates the corresponding degree cutoff
# - counts the number of nodes satisfying each cutoff
# - evaluates regulatory hubs only among genes acting as regulators
# - evaluates target hubs only among genes acting as targets
# - supports data-driven hub threshold selection
#
# INPUT:
# - node_connectivity: node table generated by
#   calculate_grn_node_connectivity()
# - quantiles: degree quantiles to evaluate
#
# OUTPUT:
# - data.frame containing:
#   - hub_type
#   - quantile
#   - degree_cutoff
#   - n_nodes
# ------------------------------------------------------------------------------

#' Evaluate candidate degree thresholds for GRN hub identification
#'
#' This function:
#' - evaluates high quantiles of out-degree, in-degree, and total degree
#' - calculates the corresponding degree cutoff
#' - counts the number of nodes satisfying each cutoff
#' - evaluates regulatory hubs only among genes acting as regulators
#' - evaluates target hubs only among genes acting as targets
#' - supports data-driven hub threshold selection
#'
#' @param node_connectivity node table generated by calculate_grn_node_connectivity()
#' @param quantiles degree quantiles to evaluate
#'
#' @return data.frame containing: hub_type quantile degree_cutoff n_nodes
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
evaluate_grn_hub_thresholds <- function(
    node_connectivity,
    quantiles = c(0.90, 0.95, 0.99, 0.995)
) {
  
  regulatory_nodes <- node_connectivity[
    node_connectivity$is_regulator,
    ,
    drop = FALSE
  ]
  
  target_nodes <- node_connectivity[
    node_connectivity$is_target,
    ,
    drop = FALSE
  ]
  
  evaluate_metric <- function(values, hub_type) {
    
    cutoffs <- as.numeric(
      quantile(
        values,
        probs = quantiles,
        na.rm = TRUE,
        names = FALSE
      )
    )
    
    do.call(
      rbind,
      lapply(seq_along(quantiles), function(i) {
        
        data.frame(
          hub_type = hub_type,
          quantile = quantiles[i],
          degree_cutoff = cutoffs[i],
          n_nodes = sum(values >= cutoffs[i]),
          stringsAsFactors = FALSE
        )
      })
    )
  }
  
  regulatory_results <- evaluate_metric(
    regulatory_nodes$out_degree,
    "regulatory"
  )
  
  target_results <- evaluate_metric(
    target_nodes$in_degree,
    "target"
  )
  
  global_results <- evaluate_metric(
    node_connectivity$total_degree,
    "global"
  )
  
  results <- rbind(
    regulatory_results,
    target_results,
    global_results
  )
  
  rownames(results) <- NULL
  
  return(results)
}

# ------------------------------------------------------------------------------
# FUNCTION: identify_grn_hubs
# ------------------------------------------------------------------------------
# Identify regulatory, target, and global hubs using degree quantiles.
#
# This function:
# - identifies regulatory hubs from the upper tail of out-degree among regulators
# - identifies target hubs from the upper tail of in-degree among targets
# - identifies global hubs from the upper tail of total degree
# - preserves node annotations and connectivity metrics
# - returns hub-specific degree cutoffs for reproducibility
#
# INPUT:
# - node_connectivity: node table generated by
#   calculate_grn_node_connectivity()
# - hub_quantile: quantile used to define high-connectivity nodes
#
# OUTPUT:
# - list containing:
#   - regulatory_hubs
#   - target_hubs
#   - global_hubs
#   - thresholds
# ------------------------------------------------------------------------------

#' Identify regulatory, target, and global hubs using degree quantiles
#'
#' This function:
#' - identifies regulatory hubs from the upper tail of out-degree among regulators
#' - identifies target hubs from the upper tail of in-degree among targets
#' - identifies global hubs from the upper tail of total degree
#' - preserves node annotations and connectivity metrics
#' - returns hub-specific degree cutoffs for reproducibility
#'
#' @param node_connectivity node table generated by calculate_grn_node_connectivity()
#' @param hub_quantile quantile used to define high-connectivity nodes
#'
#' @return list containing: regulatory_hubs target_hubs global_hubs thresholds
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
identify_grn_hubs <- function(
    node_connectivity,
    hub_quantile = 0.99
) {
  
  regulatory_nodes <- node_connectivity[
    node_connectivity$is_regulator,
    ,
    drop = FALSE
  ]
  
  target_nodes <- node_connectivity[
    node_connectivity$is_target,
    ,
    drop = FALSE
  ]
  
  regulatory_cutoff <- as.numeric(
    quantile(
      regulatory_nodes$out_degree,
      probs = hub_quantile,
      na.rm = TRUE
    )
  )
  
  target_cutoff <- as.numeric(
    quantile(
      target_nodes$in_degree,
      probs = hub_quantile,
      na.rm = TRUE
    )
  )
  
  global_cutoff <- as.numeric(
    quantile(
      node_connectivity$total_degree,
      probs = hub_quantile,
      na.rm = TRUE
    )
  )
  
  regulatory_hubs <- regulatory_nodes[
    regulatory_nodes$out_degree >= regulatory_cutoff,
    ,
    drop = FALSE
  ]
  
  target_hubs <- target_nodes[
    target_nodes$in_degree >= target_cutoff,
    ,
    drop = FALSE
  ]
  
  global_hubs <- node_connectivity[
    node_connectivity$total_degree >= global_cutoff,
    ,
    drop = FALSE
  ]
  
  regulatory_hubs <- regulatory_hubs[
    order(-regulatory_hubs$out_degree),
    ,
    drop = FALSE
  ]
  
  target_hubs <- target_hubs[
    order(-target_hubs$in_degree),
    ,
    drop = FALSE
  ]
  
  global_hubs <- global_hubs[
    order(-global_hubs$total_degree),
    ,
    drop = FALSE
  ]
  
  thresholds <- data.frame(
    hub_type = c("regulatory", "target", "global"),
    quantile = hub_quantile,
    degree_cutoff = c(
      regulatory_cutoff,
      target_cutoff,
      global_cutoff
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      regulatory_hubs = regulatory_hubs,
      target_hubs = target_hubs,
      global_hubs = global_hubs,
      thresholds = thresholds
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: analyze_grn_hub_lnc_connectivity
# ------------------------------------------------------------------------------
# Analyze lncRNA connectivity of regulatory hubs in the final GRN.
#
# This function:
# - selects interactions originating from regulatory hubs
# - counts lncRNA and protein-coding targets for each hub
# - calculates the proportion of hub targets represented by lncRNA genes
# - summarizes positive and negative lncRNA-associated interactions
# - ranks regulatory hubs by lncRNA connectivity
#
# INPUT:
# - grn_edges: final GRN edge list generated by construct_final_grn()
# - regulatory_hubs: regulatory hub table generated by identify_grn_hubs()
# - gene_sets: filtered GRN gene sets
#
# OUTPUT:
# - data.frame containing:
#   - regulator
#   - total_edges
#   - lnc_edges
#   - prot_edges
#   - lnc_edge_proportion
#   - positive_lnc_edges
#   - negative_lnc_edges
# ------------------------------------------------------------------------------

#' Analyze lncRNA connectivity of regulatory hubs in the final GRN
#'
#' This function:
#' - selects interactions originating from regulatory hubs
#' - counts lncRNA and protein-coding targets for each hub
#' - calculates the proportion of hub targets represented by lncRNA genes
#' - summarizes positive and negative lncRNA-associated interactions
#' - ranks regulatory hubs by lncRNA connectivity
#'
#' @param grn_edges final GRN edge list generated by construct_final_grn()
#' @param regulatory_hubs regulatory hub table generated by identify_grn_hubs()
#' @param gene_sets filtered GRN gene sets
#'
#' @return data.frame containing: regulator total_edges lnc_edges prot_edges lnc_edge_proportion positive_lnc_edges negative_lnc_edges
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
analyze_grn_hub_lnc_connectivity <- function(
    grn_edges,
    regulatory_hubs,
    gene_sets
) {
  
  hub_edges <- grn_edges[
    grn_edges$regulator %in% regulatory_hubs$gene,
    ,
    drop = FALSE
  ]
  
  results <- lapply(
    regulatory_hubs$gene,
    function(regulator) {
      
      regulator_edges <- hub_edges[
        hub_edges$regulator == regulator,
        ,
        drop = FALSE
      ]
      
      lnc_status <- regulator_edges$target %in% gene_sets$lnc
      prot_status <- regulator_edges$target %in% gene_sets$prot
      
      data.frame(
        regulator = regulator,
        total_edges = nrow(regulator_edges),
        lnc_edges = sum(lnc_status),
        prot_edges = sum(prot_status),
        lnc_edge_proportion = sum(lnc_status) / nrow(regulator_edges),
        positive_lnc_edges = sum(
          lnc_status &
            regulator_edges$correlation_sign == "positive"
        ),
        negative_lnc_edges = sum(
          lnc_status &
            regulator_edges$correlation_sign == "negative"
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  
  results <- do.call(rbind, results)
  
  results <- results[
    order(-results$lnc_edges),
    ,
    drop = FALSE
  ]
  
  rownames(results) <- NULL
  
  return(results)
}

# ------------------------------------------------------------------------------
# FUNCTION: analyze_grn_lnc_target_hubs
# ------------------------------------------------------------------------------
# Analyze lncRNA genes among high-connectivity GRN target hubs.
#
# This function:
# - selects lncRNA genes among identified target hubs
# - retrieves all regulatory interactions targeting each lncRNA hub
# - summarizes the number of incoming TF interactions
# - calculates average and maximum GENIE3 importance scores
# - calculates average and maximum absolute Spearman correlations
# - ranks lncRNA target hubs by in-degree
#
# INPUT:
# - grn_edges: final GRN edge list generated by construct_final_grn()
# - target_hubs: target hub table generated by identify_grn_hubs()
# - gene_sets: filtered GRN gene sets
#
# OUTPUT:
# - data.frame containing:
#   - gene
#   - in_degree
#   - n_regulators
#   - mean_genie3_weight
#   - max_genie3_weight
#   - mean_abs_spearman
#   - max_abs_spearman
# ------------------------------------------------------------------------------

#' Analyze lncRNA genes among high-connectivity GRN target hubs
#'
#' This function:
#' - selects lncRNA genes among identified target hubs
#' - retrieves all regulatory interactions targeting each lncRNA hub
#' - summarizes the number of incoming TF interactions
#' - calculates average and maximum GENIE3 importance scores
#' - calculates average and maximum absolute Spearman correlations
#' - ranks lncRNA target hubs by in-degree
#'
#' @param grn_edges final GRN edge list generated by construct_final_grn()
#' @param target_hubs target hub table generated by identify_grn_hubs()
#' @param gene_sets filtered GRN gene sets
#'
#' @return data.frame containing: gene in_degree n_regulators mean_genie3_weight max_genie3_weight mean_abs_spearman max_abs_spearman
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
analyze_grn_lnc_target_hubs <- function(
    grn_edges,
    target_hubs,
    gene_sets
) {
  
  lnc_target_hubs <- target_hubs[
    target_hubs$gene %in% gene_sets$lnc,
    ,
    drop = FALSE
  ]
  
  if (nrow(lnc_target_hubs) == 0L) {
    return(data.frame())
  }
  
  results <- lapply(
    lnc_target_hubs$gene,
    function(gene) {
      
      gene_edges <- grn_edges[
        grn_edges$target == gene,
        ,
        drop = FALSE
      ]
      
      data.frame(
        gene = gene,
        in_degree = nrow(gene_edges),
        n_regulators = length(unique(gene_edges$regulator)),
        mean_genie3_weight = mean(gene_edges$genie3_weight),
        max_genie3_weight = max(gene_edges$genie3_weight),
        mean_abs_spearman = mean(abs(gene_edges$spearman_rho)),
        max_abs_spearman = max(abs(gene_edges$spearman_rho)),
        stringsAsFactors = FALSE
      )
    }
  )
  
  results <- do.call(rbind, results)
  
  results <- results[
    order(
      -results$in_degree,
      -results$max_genie3_weight
    ),
    ,
    drop = FALSE
  ]
  
  rownames(results) <- NULL
  
  return(results)
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_grn_hub_lnc_subnetwork
# ------------------------------------------------------------------------------
# Prepare a focused regulatory hub-lncRNA subnetwork for visualization.
#
# This function:
# - selects interactions originating from identified regulatory hubs
# - retains only interactions targeting lncRNA genes
# - preserves GENIE3 importance and Spearman correlation attributes
# - annotates nodes as regulatory hubs or lncRNA targets
# - returns edge and node tables for downstream visualization
#
# INPUT:
# - grn_edges: final GRN edge list generated by construct_final_grn()
# - regulatory_hubs: regulatory hub table generated by identify_grn_hubs()
# - gene_sets: filtered GRN gene sets
#
# OUTPUT:
# - list containing:
#   - edges: regulatory hub-to-lncRNA interactions
#   - nodes: annotated nodes represented in the subnetwork
# ------------------------------------------------------------------------------

#' Prepare a focused regulatory hub-lncRNA subnetwork for visualization
#'
#' This function:
#' - selects interactions originating from identified regulatory hubs
#' - retains only interactions targeting lncRNA genes
#' - preserves GENIE3 importance and Spearman correlation attributes
#' - annotates nodes as regulatory hubs or lncRNA targets
#' - returns edge and node tables for downstream visualization
#'
#' @param grn_edges final GRN edge list generated by construct_final_grn()
#' @param regulatory_hubs regulatory hub table generated by identify_grn_hubs()
#' @param gene_sets filtered GRN gene sets
#'
#' @return list containing: edges: regulatory hub-to-lncRNA interactions nodes: annotated nodes represented in the subnetwork
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
prepare_grn_hub_lnc_subnetwork <- function(
    grn_edges,
    regulatory_hubs,
    gene_sets
) {
  
  edges <- grn_edges[
    grn_edges$regulator %in% regulatory_hubs$gene &
      grn_edges$target %in% gene_sets$lnc,
    ,
    drop = FALSE
  ]
  
  node_genes <- unique(c(
    edges$regulator,
    edges$target
  ))
  
  nodes <- data.frame(
    gene = node_genes,
    node_type = ifelse(
      node_genes %in% regulatory_hubs$gene,
      "regulatory_hub",
      "lncRNA"
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      edges = edges,
      nodes = nodes
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_focused_grn_hub_lnc_network
# ------------------------------------------------------------------------------
# Prepare a focused hub-lncRNA network for visualization.
#
# This function:
# - calculates the number of regulatory hubs connected to each lncRNA target
# - retains lncRNA targets connected to at least a specified number of hubs
# - preserves all hub-to-lncRNA interactions for retained targets
# - creates an annotated node table for visualization
#
# This filtering is intended only to reduce visual complexity and does not
# modify the final GRN used for downstream analyses.
#
# INPUT:
# - hub_lnc_subnetwork: output generated by prepare_grn_hub_lnc_subnetwork()
# - min_hub_connections: minimum number of regulatory hubs connected to an
#   lncRNA target
#
# OUTPUT:
# - list containing:
#   - edges
#   - nodes
# ------------------------------------------------------------------------------

#' Prepare a focused hub-lncRNA network for visualization
#'
#' This function:
#' - calculates the number of regulatory hubs connected to each lncRNA target
#' - retains lncRNA targets connected to at least a specified number of hubs
#' - preserves all hub-to-lncRNA interactions for retained targets
#' - creates an annotated node table for visualization
#' This filtering is intended only to reduce visual complexity and does not
#' modify the final GRN used for downstream analyses.
#'
#' @param hub_lnc_subnetwork output generated by prepare_grn_hub_lnc_subnetwork()
#' @param min_hub_connections minimum number of regulatory hubs connected to an lncRNA target
#'
#' @return list containing: edges nodes
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
prepare_focused_grn_hub_lnc_network <- function(
    hub_lnc_subnetwork,
    min_hub_connections = 2
) {
  
  edges <- hub_lnc_subnetwork$edges
  
  lnc_hub_counts <- table(edges$target)
  
  retained_lnc <- names(
    lnc_hub_counts[lnc_hub_counts >= min_hub_connections]
  )
  
  focused_edges <- edges[
    edges$target %in% retained_lnc,
    ,
    drop = FALSE
  ]
  
  node_genes <- unique(c(
    focused_edges$regulator,
    focused_edges$target
  ))
  
  nodes <- data.frame(
    gene = node_genes,
    node_type = ifelse(
      node_genes %in% focused_edges$regulator,
      "regulatory_hub",
      "lncRNA"
    ),
    stringsAsFactors = FALSE
  )
  
  nodes$degree_in_subnetwork <- vapply(
    nodes$gene,
    function(gene) {
      sum(focused_edges$regulator == gene) +
        sum(focused_edges$target == gene)
    },
    integer(1)
  )
  
  return(
    list(
      edges = focused_edges,
      nodes = nodes
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_focused_grn_hub_lnc_network
# ------------------------------------------------------------------------------
# Visualize the focused regulatory hub-lncRNA network.
#
# This function:
# - converts focused edge and node tables into a directed igraph object
# - initializes regulatory hubs across the plotting area
# - places lncRNA targets close to their connected regulatory hubs
# - refines node positions using a force-directed layout
# - distinguishes regulatory hubs from lncRNA targets by shape and color
# - labels all regulatory hubs and all lncRNA targets
# - uses automatic label repulsion to avoid text overlap
# - adds a legend describing node types
# - preserves all nodes and edges in the focused network
#
# Layout and labeling affect visualization only and do not modify the network.
#
# INPUT:
# - focused_network: output generated by
#   prepare_focused_grn_hub_lnc_network()
# - seed: random seed used for reproducible network layout
#
# OUTPUT:
# - igraph object used for network visualization
# ------------------------------------------------------------------------------

#' Visualize the focused regulatory hub-lncRNA network
#'
#' This function:
#' - converts focused edge and node tables into a directed igraph object
#' - initializes regulatory hubs across the plotting area
#' - places lncRNA targets close to their connected regulatory hubs
#' - refines node positions using a force-directed layout
#' - distinguishes regulatory hubs from lncRNA targets by shape and color
#' - labels all regulatory hubs and all lncRNA targets
#' - uses automatic label repulsion to avoid text overlap
#' - adds a legend describing node types
#' - preserves all nodes and edges in the focused network
#' Layout and labeling affect visualization only and do not modify the network.
#'
#' @param focused_network output generated by prepare_focused_grn_hub_lnc_network()
#' @param seed random seed used for reproducible network layout
#'
#' @return igraph object used for network visualization
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_focused_grn_hub_lnc_network <- function(
    focused_network,
    seed = 1234
) {
  
  edges <- focused_network$edges
  nodes <- focused_network$nodes
  
  graph <- igraph::graph_from_data_frame(
    d = edges,
    directed = TRUE,
    vertices = nodes
  )
  
  node_is_hub <- igraph::V(graph)$node_type == "regulatory_hub"
  
  hub_names <- igraph::V(graph)$name[node_is_hub]
  lnc_names <- igraph::V(graph)$name[!node_is_hub]
  
  # --------------------------------------------------------------------------
  # Initial hub-and-spoke layout
  # --------------------------------------------------------------------------
  
  set.seed(seed)
  
  initial_layout <- matrix(
    0,
    nrow = igraph::vcount(graph),
    ncol = 2,
    dimnames = list(
      igraph::V(graph)$name,
      c("x", "y")
    )
  )
  
  n_hubs <- length(hub_names)
  
  hub_angles <- seq(
    0,
    2 * pi,
    length.out = n_hubs + 1
  )[-(n_hubs + 1)]
  
  initial_layout[hub_names, "x"] <- 2.7 * cos(hub_angles)
  initial_layout[hub_names, "y"] <- 1.5 * sin(hub_angles)
  
  for (lnc in lnc_names) {
    
    connected_hubs <- unique(
      edges$regulator[
        edges$target == lnc &
          edges$regulator %in% hub_names
      ]
    )
    
    if (length(connected_hubs) > 0L) {
      
      hub_center <- colMeans(
        initial_layout[
          connected_hubs,
          ,
          drop = FALSE
        ]
      )
      
      angle <- stats::runif(
        1,
        0,
        2 * pi
      )
      
      radius <- stats::runif(
        1,
        0.35,
        0.75
      )
      
      initial_layout[lnc, "x"] <-
        hub_center["x"] +
        radius * cos(angle)
      
      initial_layout[lnc, "y"] <-
        hub_center["y"] +
        radius * sin(angle)
      
    } else {
      
      initial_layout[lnc, ] <- c(
        stats::runif(1, -2.5, 2.5),
        stats::runif(1, -1.5, 1.5)
      )
    }
  }
  
  initial_layout <- initial_layout[
    igraph::V(graph)$name,
    ,
    drop = FALSE
  ]
  
  # --------------------------------------------------------------------------
  # Force-directed refinement
  # --------------------------------------------------------------------------
  
  layout <- igraph::layout_with_fr(
    graph,
    coords = initial_layout,
    niter = 5000,
    grid = "nogrid"
  )
  
  layout <- igraph::norm_coords(
    layout,
    xmin = -3.0,
    xmax = 3.0,
    ymin = -1.8,
    ymax = 1.8
  )
  
  # --------------------------------------------------------------------------
  # Prepare node coordinates
  # --------------------------------------------------------------------------
  
  plot_nodes <- data.frame(
    gene = igraph::V(graph)$name,
    node_type = igraph::V(graph)$node_type,
    degree_in_subnetwork = igraph::V(graph)$degree_in_subnetwork,
    x = layout[, 1],
    y = layout[, 2],
    stringsAsFactors = FALSE
  )
  
  plot_nodes$is_hub <-
    plot_nodes$node_type == "regulatory_hub"
  
  # --------------------------------------------------------------------------
  # Prepare edge coordinates
  # --------------------------------------------------------------------------
  
  plot_edges <- edges
  
  plot_edges$x <- plot_nodes$x[
    match(
      plot_edges$regulator,
      plot_nodes$gene
    )
  ]
  
  plot_edges$y <- plot_nodes$y[
    match(
      plot_edges$regulator,
      plot_nodes$gene
    )
  ]
  
  plot_edges$xend <- plot_nodes$x[
    match(
      plot_edges$target,
      plot_nodes$gene
    )
  ]
  
  plot_edges$yend <- plot_nodes$y[
    match(
      plot_edges$target,
      plot_nodes$gene
    )
  ]
  
  # --------------------------------------------------------------------------
  # Separate node types for labeling
  # --------------------------------------------------------------------------
  
  hub_data <- plot_nodes[
    plot_nodes$is_hub,
    ,
    drop = FALSE
  ]
  
  lnc_data <- plot_nodes[
    !plot_nodes$is_hub,
    ,
    drop = FALSE
  ]
  
  # --------------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------------
  
  plot <- ggplot2::ggplot() +
    
    # Directed GRN edges
    ggplot2::geom_segment(
      data = plot_edges,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      linewidth = 0.35,
      color = "grey55",
      alpha = 0.55,
      arrow = grid::arrow(
        length = grid::unit(0.08, "inches"),
        type = "closed"
      )
    ) +
    
    # lncRNA target nodes
    ggplot2::geom_point(
      data = lnc_data,
      ggplot2::aes(
        x = x,
        y = y,
        shape = "lncRNA target",
        fill = "lncRNA target"
      ),
      size = 2.2,
      stroke = 0.45,
      color = "grey25"
    ) +
    
    # Regulatory hub nodes
    ggplot2::geom_point(
      data = hub_data,
      ggplot2::aes(
        x = x,
        y = y,
        shape = "Regulatory hub",
        fill = "Regulatory hub"
      ),
      size = 3.2,
      stroke = 0.55,
      color = "grey20"
    ) +
    
    # lncRNA labels
    ggrepel::geom_text_repel(
      data = lnc_data,
      ggplot2::aes(
        x = x,
        y = y,
        label = gene
      ),
      size = 2.4,
      color = "grey20",
      seed = seed,
      box.padding = 0.35,
      point.padding = 0.30,
      min.segment.length = 0,
      segment.color = "grey70",
      segment.size = 0.25,
      force = 2.5,
      force_pull = 0.4,
      max.overlaps = Inf,
      max.time = 10
    ) +
    
    # Regulatory hub labels
    ggrepel::geom_text_repel(
      data = hub_data,
      ggplot2::aes(
        x = x,
        y = y,
        label = gene
      ),
      size = 3.1,
      fontface = "bold",
      color = "grey10",
      seed = seed,
      box.padding = 0.8,
      point.padding = 0.7,
      min.segment.length = 0,
      segment.color = "grey60",
      segment.size = 0.30,
      force = 5,
      force_pull = 0.2,
      max.overlaps = Inf,
      max.time = 10
    ) +
    
    ggplot2::scale_shape_manual(
      name = NULL,
      values = c(
        "Regulatory hub" = 22,
        "lncRNA target" = 21
      ),
      breaks = c(
        "Regulatory hub",
        "lncRNA target"
      )
    ) +
    
    ggplot2::scale_fill_manual(
      name = NULL,
      values = c(
        "Regulatory hub" = "#E68600",
        "lncRNA target" = "#9CCBEA"
      ),
      breaks = c(
        "Regulatory hub",
        "lncRNA target"
      )
    ) +
    
    ggplot2::coord_cartesian(
      xlim = c(-3.7, 3.7),
      ylim = c(-2.35, 2.35),
      clip = "off"
    ) +
    
    ggplot2::labs(
      title = "Regulatory Hub-lncRNA Network"
    ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 15,
        margin = ggplot2::margin(
          b = 12
        )
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(
        size = 9
      ),
      plot.margin = ggplot2::margin(
        t = 15,
        r = 30,
        b = 15,
        l = 30
      )
    )
  
  print(plot)
  
  invisible(graph)
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_grn_lnc_target_hub_network
# ------------------------------------------------------------------------------
# Prepare a focused lncRNA target hub network for visualization.
#
# This function:
# - ranks lncRNA target hubs according to in-degree
# - selects the most connected lncRNA target hubs for visualization
# - retrieves all incoming TF interactions associated with selected lncRNA hubs
# - creates an annotated node table distinguishing TF regulators and lncRNA hubs
#
# Selection of the displayed lncRNA hubs is performed only for graphical
# readability and does not modify the final GRN or hub classification.
#
# INPUT:
# - grn_edges: final GRN edge list generated by construct_final_grn()
# - lnc_target_hubs: lncRNA target hub table generated by
#   analyze_grn_lnc_target_hubs()
# - top_n: number of highest in-degree lncRNA target hubs to display
#
# OUTPUT:
# - list containing:
#   - edges
#   - nodes
# ------------------------------------------------------------------------------

#' Prepare a focused lncRNA target hub network for visualization
#'
#' This function:
#' - ranks lncRNA target hubs according to in-degree
#' - selects the most connected lncRNA target hubs for visualization
#' - retrieves all incoming TF interactions associated with selected lncRNA hubs
#' - creates an annotated node table distinguishing TF regulators and lncRNA hubs
#' Selection of the displayed lncRNA hubs is performed only for graphical
#' readability and does not modify the final GRN or hub classification.
#'
#' @param grn_edges final GRN edge list generated by construct_final_grn()
#' @param lnc_target_hubs lncRNA target hub table generated by analyze_grn_lnc_target_hubs()
#' @param top_n number of highest in-degree lncRNA target hubs to display
#'
#' @return list containing: edges nodes
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
prepare_grn_lnc_target_hub_network <- function(
    grn_edges,
    lnc_target_hubs,
    top_n = 10
) {
  
  selected_hubs <- head(
    lnc_target_hubs$gene[
      order(-lnc_target_hubs$in_degree)
    ],
    top_n
  )
  
  edges <- grn_edges[
    grn_edges$target %in% selected_hubs,
    ,
    drop = FALSE
  ]
  
  node_genes <- unique(c(
    edges$regulator,
    edges$target
  ))
  
  nodes <- data.frame(
    gene = node_genes,
    node_type = ifelse(
      node_genes %in% selected_hubs,
      "lncRNA_target_hub",
      "TF_regulator"
    ),
    stringsAsFactors = FALSE
  )
  
  nodes$degree_in_subnetwork <- vapply(
    nodes$gene,
    function(gene) {
      sum(edges$regulator == gene) +
        sum(edges$target == gene)
    },
    integer(1)
  )
  
  return(
    list(
      edges = edges,
      nodes = nodes
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_grn_lnc_target_hub_network
# ------------------------------------------------------------------------------
# Visualize the lncRNA target hub network.
#
# This function:
# - converts target-centered edge and node tables into a directed igraph object
# - initializes lncRNA target hubs across the plotting area
# - places TF regulators close to their connected lncRNA target hubs
# - refines node positions using a force-directed layout
# - distinguishes lncRNA target hubs from TF regulators by shape and color
# - labels all lncRNA target hubs and all TF regulators
# - uses automatic label repulsion to reduce text overlap
# - adds a legend describing node types
# - preserves all nodes and edges in the prepared network
#
# Layout and labeling affect visualization only and do not modify the network.
#
# INPUT:
# - target_hub_network: output generated by
#   prepare_grn_lnc_target_hub_network()
# - seed: random seed used for reproducible network layout
#
# OUTPUT:
# - igraph object used for network visualization
# ------------------------------------------------------------------------------

#' Visualize the lncRNA target hub network
#'
#' This function:
#' - converts target-centered edge and node tables into a directed igraph object
#' - initializes lncRNA target hubs across the plotting area
#' - places TF regulators close to their connected lncRNA target hubs
#' - refines node positions using a force-directed layout
#' - distinguishes lncRNA target hubs from TF regulators by shape and color
#' - labels all lncRNA target hubs and all TF regulators
#' - uses automatic label repulsion to reduce text overlap
#' - adds a legend describing node types
#' - preserves all nodes and edges in the prepared network
#' Layout and labeling affect visualization only and do not modify the network.
#'
#' @param target_hub_network output generated by prepare_grn_lnc_target_hub_network()
#' @param seed random seed used for reproducible network layout
#'
#' @return igraph object used for network visualization
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_grn_lnc_target_hub_network <- function(
    target_hub_network,
    seed = 1234
) {
  
  edges <- target_hub_network$edges
  nodes <- target_hub_network$nodes
  
  graph <- igraph::graph_from_data_frame(
    d = edges,
    directed = TRUE,
    vertices = nodes
  )
  
  node_is_hub <-
    igraph::V(graph)$node_type == "lncRNA_target_hub"
  
  hub_names <- igraph::V(graph)$name[node_is_hub]
  
  tf_names <- igraph::V(graph)$name[!node_is_hub]
  
  # --------------------------------------------------------------------------
  # Initial target-centered layout
  # --------------------------------------------------------------------------
  
  set.seed(seed)
  
  initial_layout <- matrix(
    0,
    nrow = igraph::vcount(graph),
    ncol = 2,
    dimnames = list(
      igraph::V(graph)$name,
      c("x", "y")
    )
  )
  
  n_hubs <- length(hub_names)
  
  hub_angles <- seq(
    0,
    2 * pi,
    length.out = n_hubs + 1
  )[-(n_hubs + 1)]
  
  initial_layout[hub_names, "x"] <-
    2.7 * cos(hub_angles)
  
  initial_layout[hub_names, "y"] <-
    1.5 * sin(hub_angles)
  
  # Place each TF close to the lncRNA target hubs it regulates.
  
  for (tf in tf_names) {
    
    connected_hubs <- unique(
      edges$target[
        edges$regulator == tf &
          edges$target %in% hub_names
      ]
    )
    
    if (length(connected_hubs) > 0L) {
      
      hub_center <- colMeans(
        initial_layout[
          connected_hubs,
          ,
          drop = FALSE
        ]
      )
      
      angle <- stats::runif(
        1,
        0,
        2 * pi
      )
      
      radius <- stats::runif(
        1,
        0.55,
        1.10
      )
      
      initial_layout[tf, "x"] <-
        hub_center["x"] +
        radius * cos(angle)
      
      initial_layout[tf, "y"] <-
        hub_center["y"] +
        radius * sin(angle)
      
    } else {
      
      initial_layout[tf, ] <- c(
        stats::runif(1, -2.5, 2.5),
        stats::runif(1, -1.5, 1.5)
      )
    }
  }
  
  initial_layout <- initial_layout[
    igraph::V(graph)$name,
    ,
    drop = FALSE
  ]
  
  # --------------------------------------------------------------------------
  # Force-directed refinement
  # --------------------------------------------------------------------------
  
  layout <- igraph::layout_with_fr(
    graph,
    coords = initial_layout,
    niter = 5000,
    grid = "nogrid"
  )
  
  layout <- igraph::norm_coords(
    layout,
    xmin = -4.2,
    xmax = 4.2,
    ymin = -2.8,
    ymax = 2.8
  )
  
  # --------------------------------------------------------------------------
  # Prepare node coordinates
  # --------------------------------------------------------------------------
  
  plot_nodes <- data.frame(
    gene = igraph::V(graph)$name,
    node_type = igraph::V(graph)$node_type,
    degree_in_subnetwork =
      igraph::V(graph)$degree_in_subnetwork,
    x = layout[, 1],
    y = layout[, 2],
    stringsAsFactors = FALSE
  )
  
  plot_nodes$is_hub <-
    plot_nodes$node_type == "lncRNA_target_hub"
  
  # --------------------------------------------------------------------------
  # Prepare edge coordinates
  # --------------------------------------------------------------------------
  
  plot_edges <- edges
  
  plot_edges$x <- plot_nodes$x[
    match(
      plot_edges$regulator,
      plot_nodes$gene
    )
  ]
  
  plot_edges$y <- plot_nodes$y[
    match(
      plot_edges$regulator,
      plot_nodes$gene
    )
  ]
  
  plot_edges$xend <- plot_nodes$x[
    match(
      plot_edges$target,
      plot_nodes$gene
    )
  ]
  
  plot_edges$yend <- plot_nodes$y[
    match(
      plot_edges$target,
      plot_nodes$gene
    )
  ]
  
  # --------------------------------------------------------------------------
  # Separate node types
  # --------------------------------------------------------------------------
  
  hub_data <- plot_nodes[
    plot_nodes$is_hub,
    ,
    drop = FALSE
  ]
  
  tf_data <- plot_nodes[
    !plot_nodes$is_hub,
    ,
    drop = FALSE
  ]
  
  # --------------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------------
  
  plot <- ggplot2::ggplot() +
    
    # Directed TF -> lncRNA edges
    ggplot2::geom_segment(
      data = plot_edges,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      linewidth = 0.35,
      color = "grey55",
      alpha = 0.55,
      arrow = grid::arrow(
        length = grid::unit(0.08, "inches"),
        type = "closed"
      )
    ) +
    
    # TF regulator nodes
    ggplot2::geom_point(
      data = tf_data,
      ggplot2::aes(
        x = x,
        y = y,
        shape = "TF regulator",
        fill = "TF regulator"
      ),
      size = 2.2,
      stroke = 0.45,
      color = "grey25"
    ) +
    
    # lncRNA target hub nodes
    ggplot2::geom_point(
      data = hub_data,
      ggplot2::aes(
        x = x,
        y = y,
        shape = "lncRNA target hub",
        fill = "lncRNA target hub"
      ),
      size = 3.2,
      stroke = 0.55,
      color = "grey20"
    ) +
    
    # ------------------------------------------------------------------------
  # All TF regulator labels
  # ------------------------------------------------------------------------
  
  ggrepel::geom_text_repel(
    data = tf_data,
    ggplot2::aes(
      x = x,
      y = y,
      label = gene
    ),
    size = 2.1,
    color = "grey25",
    seed = seed,
    box.padding = 0.30,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.color = "grey75",
    segment.size = 0.20,
    force = 4,
    force_pull = 0.15,
    max.overlaps = Inf,
    max.time = 30,
    max.iter = 30000
  ) +
    
    # ------------------------------------------------------------------------
  # All lncRNA target hub labels
  # ------------------------------------------------------------------------
  
  ggrepel::geom_text_repel(
    data = hub_data,
    ggplot2::aes(
      x = x,
      y = y,
      label = gene
    ),
    size = 3.1,
    fontface = "bold",
    color = "grey10",
    seed = seed,
    box.padding = 0.9,
    point.padding = 0.8,
    min.segment.length = 0,
    segment.color = "grey60",
    segment.size = 0.30,
    force = 6,
    force_pull = 0.15,
    max.overlaps = Inf,
    max.time = 30,
    max.iter = 30000
  ) +
    
    # ------------------------------------------------------------------------
  # Node legend
  # ------------------------------------------------------------------------
  
  ggplot2::scale_shape_manual(
    name = NULL,
    values = c(
      "lncRNA target hub" = 22,
      "TF regulator" = 21
    ),
    breaks = c(
      "lncRNA target hub",
      "TF regulator"
    )
  ) +
    
    ggplot2::scale_fill_manual(
      name = NULL,
      values = c(
        "lncRNA target hub" = "#E68600",
        "TF regulator" = "#9CCBEA"
      ),
      breaks = c(
        "lncRNA target hub",
        "TF regulator"
      )
    ) +
    
    ggplot2::coord_cartesian(
      xlim = c(-5.0, 5.0),
      ylim = c(-3.5, 3.5),
      clip = "off"
    ) +
    
    ggplot2::labs(
      title = "lncRNA Target Hub Network"
    ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = 15,
        margin = ggplot2::margin(
          b = 12
        )
      ),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(
        size = 9
      ),
      plot.margin = ggplot2::margin(
        t = 15,
        r = 30,
        b = 15,
        l = 30
      )
    )
  
  print(plot)
  
  invisible(graph)
}

# ------------------------------------------------------------------------------
# FUNCTION: export_grn_core_results
# ------------------------------------------------------------------------------
# Export principal GRN result tables in tab-separated text format.
#
# This function:
# - creates a dedicated results table directory
# - exports the final GRN edge list
# - exports annotated GRN nodes and node connectivity
# - exports regulatory, target, and global hub tables
# - exports regulatory hub-lncRNA connectivity results
# - exports lncRNA target hub results
#
# INPUT:
# - grn_edges: final GRN edge list
# - grn_nodes: annotated GRN node table
# - node_connectivity: GRN node connectivity table
# - hubs: hub-analysis result generated by identify_grn_hubs()
# - hub_lnc_connectivity: regulatory hub-lncRNA connectivity table
# - lnc_target_hubs: lncRNA target hub table
# - output_dir: directory in which the exported tables are saved
#
# OUTPUT:
# - named character vector containing the paths of exported files
# ------------------------------------------------------------------------------

#' Export principal GRN result tables in tab-separated text format
#'
#' This function:
#' - creates a dedicated results table directory
#' - exports the final GRN edge list
#' - exports annotated GRN nodes and node connectivity
#' - exports regulatory, target, and global hub tables
#' - exports regulatory hub-lncRNA connectivity results
#' - exports lncRNA target hub results
#'
#' @param grn_edges final GRN edge list
#' @param grn_nodes annotated GRN node table
#' @param node_connectivity GRN node connectivity table
#' @param hubs hub-analysis result generated by identify_grn_hubs()
#' @param hub_lnc_connectivity regulatory hub-lncRNA connectivity table
#' @param lnc_target_hubs lncRNA target hub table
#' @param output_dir directory in which the exported tables are saved
#'
#' @return named character vector containing the paths of exported files
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
export_grn_core_results <- function(
    grn_edges,
    grn_nodes,
    node_connectivity,
    hubs,
    hub_lnc_connectivity,
    lnc_target_hubs,
    output_dir
) {
  
  table_dir <- file.path(
    output_dir,
    "tables"
  )
  
  dir.create(
    table_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  tables <- list(
    grn_final_edges = grn_edges,
    grn_nodes = grn_nodes,
    grn_node_connectivity = node_connectivity,
    grn_regulatory_hubs = hubs$regulatory_hubs,
    grn_target_hubs = hubs$target_hubs,
    grn_global_hubs = hubs$global_hubs,
    grn_hub_lnc_connectivity = hub_lnc_connectivity,
    grn_lnc_target_hubs = lnc_target_hubs
  )
  
  file_paths <- vapply(
    names(tables),
    function(table_name) {
      
      file_path <- file.path(
        table_dir,
        paste0(table_name, ".tsv")
      )
      
      utils::write.table(
        tables[[table_name]],
        file = file_path,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE
      )
      
      file_path
    },
    character(1)
  )
  
  return(file_paths)
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_final_grn
# ------------------------------------------------------------------------------
# Summarize the principal structural properties of the final GRN.
#
# This function:
# - counts nodes and edges
# - counts unique regulators and targets
# - summarizes positive and negative regulatory interactions
# - summarizes lncRNA and protein-coding node composition
# - reports the number of regulatory, target, and global hubs
# - reports the number of lncRNA target hubs
#
# INPUT:
# - grn_edges: final GRN edge list
# - grn_nodes: annotated GRN node table
# - hubs: hub-analysis result generated by identify_grn_hubs()
# - lnc_target_hubs: lncRNA target hub table
#
# OUTPUT:
# - data.frame containing GRN summary statistics
# ------------------------------------------------------------------------------

#' Summarize the principal structural properties of the final GRN
#'
#' This function:
#' - counts nodes and edges
#' - counts unique regulators and targets
#' - summarizes positive and negative regulatory interactions
#' - summarizes lncRNA and protein-coding node composition
#' - reports the number of regulatory, target, and global hubs
#' - reports the number of lncRNA target hubs
#'
#' @param grn_edges final GRN edge list
#' @param grn_nodes annotated GRN node table
#' @param hubs hub-analysis result generated by identify_grn_hubs()
#' @param lnc_target_hubs lncRNA target hub table
#'
#' @return data.frame containing GRN summary statistics
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
summarize_final_grn <- function(
    grn_edges,
    grn_nodes,
    hubs,
    lnc_target_hubs
) {
  
  summary_table <- data.frame(
    metric = c(
      "n_edges",
      "n_nodes",
      "n_regulators",
      "n_targets",
      "n_positive_edges",
      "n_negative_edges",
      "n_lncRNA_nodes",
      "n_protein_coding_nodes",
      "n_TF_nodes",
      "n_regulatory_hubs",
      "n_target_hubs",
      "n_global_hubs",
      "n_lncRNA_target_hubs"
    ),
    value = c(
      nrow(grn_edges),
      nrow(grn_nodes),
      length(unique(grn_edges$regulator)),
      length(unique(grn_edges$target)),
      sum(grn_edges$correlation_sign == "positive"),
      sum(grn_edges$correlation_sign == "negative"),
      sum(grn_nodes$gene_type == "lncRNA"),
      sum(grn_nodes$gene_type == "protein_coding"),
      sum(grn_nodes$is_tf),
      nrow(hubs$regulatory_hubs),
      nrow(hubs$target_hubs),
      nrow(hubs$global_hubs),
      nrow(lnc_target_hubs)
    ),
    stringsAsFactors = FALSE
  )
  
  return(summary_table)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_grn_detection_thresholds
# ------------------------------------------------------------------------------
# Visualize the effect of pseudobulk detection thresholds on GRN gene retention.
#
# This function:
# - compares retained lncRNA and protein-coding gene counts across thresholds
# - compares retention proportions relative to the original gene sets
# - highlights the selected detection threshold
# - supports visual justification of expression-filtering threshold selection
#
# INPUT:
# - detection_thresholds: data.frame generated by
#   evaluate_grn_detection_thresholds()
# - selected_threshold: detection threshold selected for final GRN filtering
#
# OUTPUT:
# - list containing:
#   - retained_counts: retained gene count plot
#   - retention_proportion: retained gene proportion plot
# ------------------------------------------------------------------------------

#' Visualize the effect of pseudobulk detection thresholds on GRN gene retention
#'
#' This function:
#' - compares retained lncRNA and protein-coding gene counts across thresholds
#' - compares retention proportions relative to the original gene sets
#' - highlights the selected detection threshold
#' - supports visual justification of expression-filtering threshold selection
#'
#' @param detection_thresholds data.frame generated by evaluate_grn_detection_thresholds()
#' @param selected_threshold detection threshold selected for final GRN filtering
#'
#' @return list containing: retained_counts: retained gene count plot retention_proportion: retained gene proportion plot
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_grn_detection_thresholds <- function(
    detection_thresholds,
    selected_threshold = 2
) {
  
  required_columns <- c(
    "threshold",
    "retained_lnc",
    "retained_prot",
    "lnc_retention",
    "prot_retention"
  )
  
  if (!all(required_columns %in% colnames(detection_thresholds))) {
    stop(
      "detection_thresholds must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  count_data <- rbind(
    data.frame(
      threshold = detection_thresholds$threshold,
      gene_type = "lncRNA",
      retained_genes = detection_thresholds$retained_lnc
    ),
    data.frame(
      threshold = detection_thresholds$threshold,
      gene_type = "protein_coding",
      retained_genes = detection_thresholds$retained_prot
    )
  )
  
  retention_data <- rbind(
    data.frame(
      threshold = detection_thresholds$threshold,
      gene_type = "lncRNA",
      retention = detection_thresholds$lnc_retention
    ),
    data.frame(
      threshold = detection_thresholds$threshold,
      gene_type = "protein_coding",
      retention = detection_thresholds$prot_retention
    )
  )
  
  retained_counts_plot <- ggplot2::ggplot(
    count_data,
    ggplot2::aes(
      x = threshold,
      y = retained_genes,
      color = gene_type,
      group = gene_type
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(
      xintercept = selected_threshold,
      linetype = "dashed"
    ) +
    ggplot2::scale_x_continuous(
      breaks = detection_thresholds$threshold
    ) +
    ggplot2::labs(
      x = "Minimum number of pseudobulk samples with detected expression",
      y = "Retained genes",
      color = "Gene type",
      title = "Gene Retention Across Expression Detection Thresholds"
    ) +
    ggplot2::theme_classic()
  
  retention_proportion_plot <- ggplot2::ggplot(
    retention_data,
    ggplot2::aes(
      x = threshold,
      y = retention * 100,
      color = gene_type,
      group = gene_type
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(
      xintercept = selected_threshold,
      linetype = "dashed"
    ) +
    ggplot2::scale_x_continuous(
      breaks = detection_thresholds$threshold
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 100)
    ) +
    ggplot2::labs(
      x = "Minimum number of pseudobulk samples with detected expression",
      y = "Genes retained (%)",
      color = "Gene type",
      title = "Gene Retention Proportion Across Expression Detection Thresholds"
    ) +
    ggplot2::theme_classic()
  
  return(
    list(
      retained_counts = retained_counts_plot,
      retention_proportion = retention_proportion_plot
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_grn_spearman_distribution
# ------------------------------------------------------------------------------
# Visualize the distribution of unfiltered pairwise Spearman correlations.
#
# This function:
# - extracts unique gene-gene correlations from the upper matrix triangle
# - excludes self-correlations from the matrix diagonal
# - visualizes the complete unfiltered correlation distribution
# - marks commonly evaluated absolute correlation thresholds
#
# INPUT:
# - correlation_matrix: unfiltered Spearman correlation matrix
# - thresholds: absolute Spearman correlation thresholds to display
#
# OUTPUT:
# - ggplot object showing the Spearman correlation distribution
# ------------------------------------------------------------------------------

#' Visualize the distribution of unfiltered pairwise Spearman correlations
#'
#' This function:
#' - extracts unique gene-gene correlations from the upper matrix triangle
#' - excludes self-correlations from the matrix diagonal
#' - visualizes the complete unfiltered correlation distribution
#' - marks commonly evaluated absolute correlation thresholds
#'
#' @param correlation_matrix unfiltered Spearman correlation matrix
#' @param thresholds absolute Spearman correlation thresholds to display
#'
#' @return ggplot object showing the Spearman correlation distribution
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_grn_spearman_distribution <- function(
    correlation_matrix,
    thresholds = c(0.3, 0.4, 0.5, 0.6, 0.7)
) {
  
  correlation_values <- correlation_matrix[
    upper.tri(correlation_matrix)
  ]
  
  correlation_values <- correlation_values[
    is.finite(correlation_values)
  ]
  
  plot_data <- data.frame(
    spearman_rho = correlation_values
  )
  
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = spearman_rho)
  ) +
    ggplot2::geom_histogram(
      bins = 100,
      boundary = 0
    ) +
    ggplot2::geom_vline(
      xintercept = c(-thresholds, thresholds),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(-1, 1, by = 0.2)
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-1, 1)
    ) +
    ggplot2::labs(
      x = "Spearman correlation (rho)",
      y = "Number of gene pairs",
      title = "Distribution of Unfiltered Spearman Correlations"
    ) +
    ggplot2::theme_classic()
  
  return(plot)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_grn_genie3_importance_distribution
# ------------------------------------------------------------------------------
# Visualize the distribution of GENIE3 importance scores.
#
# This function:
# - extracts GENIE3 importance scores from the integrated edge table
# - visualizes the complete importance-score distribution
# - calculates selected upper quantiles
# - marks quantile-based thresholds used for candidate network evaluation
# - labels the displayed quantile thresholds
#
# INPUT:
# - integrated_edges: integrated GENIE3-Spearman edge table
# - quantiles: GENIE3 importance quantiles to display
#
# OUTPUT:
# - ggplot object showing the GENIE3 importance distribution
# ------------------------------------------------------------------------------

#' Visualize the distribution of GENIE3 importance scores
#'
#' This function:
#' - extracts GENIE3 importance scores from the integrated edge table
#' - visualizes the complete importance-score distribution
#' - calculates selected upper quantiles
#' - marks quantile-based thresholds used for candidate network evaluation
#' - labels the displayed quantile thresholds
#'
#' @param integrated_edges integrated GENIE3-Spearman edge table
#' @param quantiles GENIE3 importance quantiles to display
#'
#' @return ggplot object showing the GENIE3 importance distribution
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_grn_genie3_importance_distribution <- function(
    integrated_edges,
    quantiles = c(0.90, 0.95, 0.99, 0.995, 0.999)
) {
  
  if (!"genie3_weight" %in% colnames(integrated_edges)) {
    stop("integrated_edges must contain a 'genie3_weight' column")
  }
  
  importance_values <- integrated_edges$genie3_weight
  
  importance_values <- importance_values[
    is.finite(importance_values)
  ]
  
  importance_values <- importance_values[
    importance_values >= 0
  ]
  
  quantile_values <- stats::quantile(
    importance_values,
    probs = quantiles,
    names = FALSE
  )
  
  quantile_data <- data.frame(
    quantile = quantiles,
    cutoff = quantile_values
  )
  
  quantile_data$label <- paste0(
    quantile_data$quantile * 100,
    "%"
  )
  
  plot_data <- data.frame(
    genie3_weight = importance_values
  )
  
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = genie3_weight)
  ) +
    ggplot2::geom_histogram(
      bins = 100
    ) +
    ggplot2::geom_vline(
      data = quantile_data,
      ggplot2::aes(xintercept = cutoff),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_text(
      data = quantile_data,
      ggplot2::aes(
        x = cutoff,
        y = Inf,
        label = label
      ),
      angle = 90,
      vjust = 1.5,
      hjust = 1.1,
      size = 3,
      inherit.aes = FALSE
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, 0.025)
    ) +
    ggplot2::labs(
      x = "GENIE3 importance score",
      y = "Number of candidate TF-target interactions",
      title = "Distribution of GENIE3 Importance Scores"
    ) +
    ggplot2::theme_classic()
  
  return(plot)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_grn_threshold_selection
# ------------------------------------------------------------------------------
# Visualize the effect of combined GENIE3 and Spearman thresholds on network size.
#
# This function:
# - represents GENIE3 quantiles and absolute Spearman thresholds as a grid
# - shows the number of retained edges for each threshold combination
# - uses a logarithmic fill scale to accommodate large differences in network size
# - highlights the threshold combination selected for final GRN construction
#
# INPUT:
# - threshold_counts: threshold comparison table generated by
#   summarize_grn_integrated_edges()
# - selected_genie3_quantile: GENIE3 quantile selected for final GRN construction
# - selected_spearman_threshold: absolute Spearman threshold selected for the
#   final GRN
#
# OUTPUT:
# - ggplot object showing network size across threshold combinations
# ------------------------------------------------------------------------------

#' Visualize the effect of combined GENIE3 and Spearman thresholds on network size
#'
#' This function:
#' - represents GENIE3 quantiles and absolute Spearman thresholds as a grid
#' - shows the number of retained edges for each threshold combination
#' - uses a logarithmic fill scale to accommodate large differences in network size
#' - highlights the threshold combination selected for final GRN construction
#'
#' @param threshold_counts threshold comparison table generated by summarize_grn_integrated_edges()
#' @param selected_genie3_quantile GENIE3 quantile selected for final GRN construction
#' @param selected_spearman_threshold absolute Spearman threshold selected for the final GRN
#'
#' @return ggplot object showing network size across threshold combinations
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_grn_threshold_selection <- function(
    threshold_counts,
    selected_genie3_quantile = 0.995,
    selected_spearman_threshold = 0.5
) {
  
  required_columns <- c(
    "genie3_quantile",
    "spearman_threshold",
    "n_edges"
  )
  
  if (!all(required_columns %in% colnames(threshold_counts))) {
    stop(
      "threshold_counts must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  plot_data <- threshold_counts
  
  plot_data$genie3_label <- paste0(
    plot_data$genie3_quantile * 100,
    "%"
  )
  
  plot_data$genie3_label <- factor(
    plot_data$genie3_label,
    levels = paste0(sort(unique(plot_data$genie3_quantile)) * 100, "%")
  )
  
  plot_data$spearman_label <- as.character(
    plot_data$spearman_threshold
  )
  
  plot_data$selected <-
    plot_data$genie3_quantile == selected_genie3_quantile &
    plot_data$spearman_threshold == selected_spearman_threshold
  
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = spearman_label,
      y = genie3_label,
      fill = n_edges
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = format(n_edges, big.mark = ",")),
      size = 3
    ) +
    ggplot2::geom_tile(
      data = plot_data[plot_data$selected, , drop = FALSE],
      fill = NA,
      linewidth = 1.2
    ) +
    ggplot2::scale_fill_viridis_c(
      trans = "log10",
      name = "Retained edges"
    ) +
    ggplot2::labs(
      x = "Absolute Spearman correlation threshold",
      y = "GENIE3 importance quantile",
      title = "GRN Network Size Across GENIE3 and Spearman Thresholds"
    ) +
    ggplot2::theme_classic()
  
  return(plot)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_grn_lnc_retention_vs_network_size
# ------------------------------------------------------------------------------
# Visualize the trade-off between network size and lncRNA target retention.
#
# This function:
# - compares retained network size across candidate threshold combinations
# - represents the proportion of retained lncRNA targets
# - distinguishes candidate networks according to GENIE3 quantile
# - highlights the threshold combination selected for final GRN construction
#
# INPUT:
# - candidate_networks: output generated by evaluate_grn_candidate_networks()
# - selected_genie3_quantile: GENIE3 quantile selected for final GRN construction
# - selected_spearman_threshold: absolute Spearman threshold selected for the
#   final GRN
#
# OUTPUT:
# - ggplot object showing the trade-off between network size and lncRNA retention
# ------------------------------------------------------------------------------

#' Visualize the trade-off between network size and lncRNA target retention
#'
#' This function:
#' - compares retained network size across candidate threshold combinations
#' - represents the proportion of retained lncRNA targets
#' - distinguishes candidate networks according to GENIE3 quantile
#' - highlights the threshold combination selected for final GRN construction
#'
#' @param candidate_networks output generated by evaluate_grn_candidate_networks()
#' @param selected_genie3_quantile GENIE3 quantile selected for final GRN construction
#' @param selected_spearman_threshold absolute Spearman threshold selected for the final GRN
#'
#' @return ggplot object showing the trade-off between network size and lncRNA retention
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_grn_lnc_retention_vs_network_size <- function(
    candidate_networks,
    selected_genie3_quantile = 0.995,
    selected_spearman_threshold = 0.5
) {
  
  required_columns <- c(
    "genie3_quantile",
    "spearman_threshold",
    "n_edges",
    "lnc_target_retention"
  )
  
  if (!all(required_columns %in% colnames(candidate_networks))) {
    stop(
      "candidate_networks must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  plot_data <- candidate_networks
  
  plot_data$genie3_label <- paste0(
    plot_data$genie3_quantile * 100,
    "%"
  )
  
  plot_data$genie3_label <- factor(
    plot_data$genie3_label,
    levels = paste0(sort(unique(plot_data$genie3_quantile)) * 100, "%")
  )
  
  plot_data$selected <-
    plot_data$genie3_quantile == selected_genie3_quantile &
    plot_data$spearman_threshold == selected_spearman_threshold
  
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = n_edges,
      y = lnc_target_retention * 100,
      color = genie3_label,
      group = genie3_label
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      size = 2.5
    ) +
    ggplot2::geom_point(
      data = plot_data[plot_data$selected, , drop = FALSE],
      shape = 21,
      size = 5,
      stroke = 1.2,
      fill = NA,
      color = "black"
    ) +
    ggplot2::scale_x_log10(
      labels = scales::label_comma()
    ) +
    ggplot2::labs(
      x = "Number of retained edges",
      y = "lncRNA targets retained (%)",
      color = "GENIE3 quantile",
      title = "lncRNA Retention vs GRN Network Size"
    ) +
    ggplot2::theme_classic()
  
  return(plot)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_final_grn_gene_composition
# ------------------------------------------------------------------------------
# Visualize gene-biotype composition of the final GRN.
#
# This function:
# - extracts lncRNA and protein-coding node counts from the final GRN summary
# - calculates the corresponding proportions
# - visualizes final GRN gene composition
# - annotates bars with node counts and percentages
#
# INPUT:
# - final_summary: summary table generated by summarize_final_grn()
#
# OUTPUT:
# - ggplot object showing final GRN gene-biotype composition
# ------------------------------------------------------------------------------

#' Visualize gene-biotype composition of the final GRN
#'
#' This function:
#' - extracts lncRNA and protein-coding node counts from the final GRN summary
#' - calculates the corresponding proportions
#' - visualizes final GRN gene composition
#' - annotates bars with node counts and percentages
#'
#' @param final_summary summary table generated by summarize_final_grn()
#'
#' @return ggplot object showing final GRN gene-biotype composition
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_final_grn_gene_composition <- function(
    final_summary
) {
  
  required_metrics <- c(
    "n_lncRNA_nodes",
    "n_protein_coding_nodes"
  )
  
  if (!all(required_metrics %in% final_summary$metric)) {
    stop(
      "final_summary must contain: ",
      paste(required_metrics, collapse = ", ")
    )
  }
  
  lnc_count <- final_summary$value[
    final_summary$metric == "n_lncRNA_nodes"
  ]
  
  prot_count <- final_summary$value[
    final_summary$metric == "n_protein_coding_nodes"
  ]
  
  plot_data <- data.frame(
    gene_type = c(
      "lncRNA",
      "Protein-coding"
    ),
    n_nodes = c(
      lnc_count,
      prot_count
    ),
    stringsAsFactors = FALSE
  )
  
  plot_data$proportion <- plot_data$n_nodes / sum(plot_data$n_nodes)
  
  plot_data$label <- paste0(
    format(plot_data$n_nodes, big.mark = ","),
    "\n(",
    round(plot_data$proportion * 100, 1),
    "%)"
  )
  
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = gene_type,
      y = n_nodes,
      fill = gene_type
    )
  ) +
    ggplot2::geom_col(
      width = 0.65
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      vjust = -0.4,
      size = 4
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0.18))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Number of nodes",
      title = "Gene Biotype Composition of the Final GRN"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none"
    )
  
  return(plot)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_grn_regulatory_hubs
# ------------------------------------------------------------------------------
# Visualize regulatory hubs identified in the final GRN.
#
# This function:
# - ranks regulatory hubs according to out-degree
# - displays the number of outgoing regulatory interactions for each hub
# - preserves the hub definition generated by identify_grn_hubs()
# - provides a compact visualization of the major GRN regulators
#
# INPUT:
# - regulatory_hubs: regulatory hub table generated by identify_grn_hubs()
#
# OUTPUT:
# - ggplot object showing regulatory hub out-degree
# ------------------------------------------------------------------------------

#' Visualize regulatory hubs identified in the final GRN
#'
#' This function:
#' - ranks regulatory hubs according to out-degree
#' - displays the number of outgoing regulatory interactions for each hub
#' - preserves the hub definition generated by identify_grn_hubs()
#' - provides a compact visualization of the major GRN regulators
#'
#' @param regulatory_hubs regulatory hub table generated by identify_grn_hubs()
#'
#' @return ggplot object showing regulatory hub out-degree
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_grn_regulatory_hubs <- function(
    regulatory_hubs
) {
  
  required_columns <- c(
    "gene",
    "out_degree"
  )
  
  if (!all(required_columns %in% colnames(regulatory_hubs))) {
    stop(
      "regulatory_hubs must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  plot_data <- regulatory_hubs
  
  plot_data <- plot_data[
    order(plot_data$out_degree, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  
  plot_data$gene <- factor(
    plot_data$gene,
    levels = rev(plot_data$gene)
  )
  
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = gene,
      y = out_degree
    )
  ) +
    ggplot2::geom_col(
      width = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = out_degree),
      hjust = -0.2,
      size = 3.5
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Out-degree",
      title = "Regulatory Hubs in the Final GRN"
    ) +
    ggplot2::theme_classic()
  
  return(plot)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_grn_hub_lnc_connectivity
# ------------------------------------------------------------------------------
# Visualize lncRNA connectivity among regulatory hubs.
#
# This function:
# - compares lncRNA-targeting connectivity across regulatory hubs
# - displays the number of outgoing edges targeting lncRNA genes
# - annotates each hub with the proportion of its regulatory interactions
#   directed toward lncRNA targets
# - ranks hubs according to lncRNA connectivity
#
# INPUT:
# - hub_lnc_connectivity: output generated by
#   analyze_grn_hub_lnc_connectivity()
#
# OUTPUT:
# - ggplot object showing regulatory hub connectivity toward lncRNA targets
# ------------------------------------------------------------------------------

#' Visualize lncRNA connectivity among regulatory hubs
#'
#' This function:
#' - compares lncRNA-targeting connectivity across regulatory hubs
#' - displays the number of outgoing edges targeting lncRNA genes
#' - annotates each hub with the proportion of its regulatory interactions
#' directed toward lncRNA targets
#' - ranks hubs according to lncRNA connectivity
#'
#' @param hub_lnc_connectivity output generated by analyze_grn_hub_lnc_connectivity()
#'
#' @return ggplot object showing regulatory hub connectivity toward lncRNA targets
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_grn_hub_lnc_connectivity <- function(
    hub_lnc_connectivity
) {
  
  required_columns <- c(
    "regulator",
    "lnc_edges",
    "lnc_edge_proportion"
  )
  
  if (!all(required_columns %in% colnames(hub_lnc_connectivity))) {
    stop(
      "hub_lnc_connectivity must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  plot_data <- hub_lnc_connectivity
  
  plot_data <- plot_data[
    order(plot_data$lnc_edges, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  
  plot_data$regulator <- factor(
    plot_data$regulator,
    levels = rev(plot_data$regulator)
  )
  
  plot_data$label <- paste0(
    plot_data$lnc_edges,
    " (",
    round(plot_data$lnc_edge_proportion * 100, 1),
    "%)"
  )
  
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = regulator,
      y = lnc_edges
    )
  ) +
    ggplot2::geom_col(
      width = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = label),
      hjust = -0.15,
      size = 3.3
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.22))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Number of lncRNA-targeting edges",
      title = "lncRNA Connectivity of Regulatory Hubs in the Final GRN"
    ) +
    ggplot2::theme_classic()
  
  return(plot)
}

# ==============================================================================
# INTEGRATED lncRNA-TF-PCG GRN
# ==============================================================================

# ------------------------------------------------------------------------------
# FUNCTION: infer_integrated_grn_relationships
# ------------------------------------------------------------------------------

# Infer regulator-target relationships for one interaction class.
#
# This function:
# - restricts regulators and targets to genes present in the expression matrix
# - runs GENIE3 using the specified regulator and target gene sets
# - integrates GENIE3 importance scores with the corresponding Spearman
#   correlations
# - annotates regulator type, target type, and interaction class
# - does not apply final GENIE3 or Spearman thresholds
#
# INPUT:
# - expression_matrix: filtered normalized pseudobulk expression matrix
# - spearman_matrix: filtered gene-gene Spearman correlation matrix
# - regulators: genes used as candidate regulators
# - targets: genes used as candidate targets
# - regulator_type: label describing the regulator class
# - target_type: label describing the target class
# - edge_type: label describing the interaction class
# - seed: random seed used for reproducible GENIE3 inference
#
# OUTPUT:
# - list containing:
#   - weights: GENIE3 importance matrix
#   - edges: integrated and annotated regulator-target interactions
# ------------------------------------------------------------------------------

#' Infer regulator-target relationships for one interaction class
#'
#' This function:
#' - restricts regulators and targets to genes present in the expression matrix
#' - runs GENIE3 using the specified regulator and target gene sets
#' - integrates GENIE3 importance scores with the corresponding Spearman
#' correlations
#' - annotates regulator type, target type, and interaction class
#' - does not apply final GENIE3 or Spearman thresholds
#'
#' @param expression_matrix filtered normalized pseudobulk expression matrix
#' @param spearman_matrix filtered gene-gene Spearman correlation matrix
#' @param regulators genes used as candidate regulators
#' @param targets genes used as candidate targets
#' @param regulator_type label describing the regulator class
#' @param target_type label describing the target class
#' @param edge_type label describing the interaction class
#' @param seed random seed used for reproducible GENIE3 inference
#' @param n_cores Number of CPU cores passed to `GENIE3::GENIE3()`.
#'
#' @return list containing: weights: GENIE3 importance matrix edges: integrated and annotated regulator-target interactions
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
infer_integrated_grn_relationships <- function(
    expression_matrix,
    spearman_matrix,
    regulators,
    targets,
    regulator_type,
    target_type,
    edge_type,
    seed = 1234,
    n_cores = 1
) {
  
  regulators <- intersect(regulators, rownames(expression_matrix))
  targets <- intersect(targets, rownames(expression_matrix))
  
  if (length(regulators) == 0L) {
    stop("No candidate regulators are present in the expression matrix")
  }
  
  if (length(targets) == 0L) {
    stop("No candidate targets are present in the expression matrix")
  }
  
  if (!all(regulators %in% rownames(spearman_matrix))) {
    stop("Some regulators are missing from the Spearman matrix")
  }
  
  if (!all(targets %in% colnames(spearman_matrix))) {
    stop("Some targets are missing from the Spearman matrix")
  }
  
  set.seed(seed)
  
  weights <- GENIE3::GENIE3(
    exprMatrix = as.matrix(expression_matrix),
    regulators = regulators,
    targets = targets,
    nCores = n_cores
  )
  
  edges <- integrate_genie3_spearman(
    genie3_weights = weights,
    spearman_matrix = spearman_matrix
  )
  
  edges$regulator_type <- regulator_type
  edges$target_type <- target_type
  edges$edge_type <- edge_type
  
  return(list(
    weights = weights,
    edges = edges
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: inspect_integrated_grn_score_distribution
# ------------------------------------------------------------------------------

# Inspect GENIE3 and Spearman score distributions for one interaction class.
#
# This function:
# - validates the required edge-table columns
# - summarizes the number of candidate interactions
# - calculates GENIE3 importance-score quantiles
# - calculates signed Spearman correlation quantiles
# - calculates absolute Spearman correlation quantiles
# - prepares histogram data for lightweight distribution visualization
# - does not remove or modify any network edges
#
# INPUT:
# - edges: integrated GRN interaction table
# - edge_type: interaction class label
# - genie3_breaks: number of histogram bins for GENIE3 scores
# - spearman_breaks: number of histogram bins for Spearman correlations
#
# OUTPUT:
# - list containing:
#   - summary: general interaction-class summary
#   - genie3_quantiles: GENIE3 importance-score quantiles
#   - spearman_quantiles: signed Spearman quantiles
#   - abs_spearman_quantiles: absolute Spearman quantiles
#   - genie3_histogram: histogram information for GENIE3 scores
#   - spearman_histogram: histogram information for Spearman rho
# ------------------------------------------------------------------------------

#' Inspect GENIE3 and Spearman score distributions for one interaction class
#'
#' This function:
#' - validates the required edge-table columns
#' - summarizes the number of candidate interactions
#' - calculates GENIE3 importance-score quantiles
#' - calculates signed Spearman correlation quantiles
#' - calculates absolute Spearman correlation quantiles
#' - prepares histogram data for lightweight distribution visualization
#' - does not remove or modify any network edges
#'
#' @param edges integrated GRN interaction table
#' @param edge_type interaction class label
#' @param genie3_breaks number of histogram bins for GENIE3 scores
#' @param spearman_breaks number of histogram bins for Spearman correlations
#'
#' @return list containing: summary: general interaction-class summary genie3_quantiles: GENIE3 importance-score quantiles spearman_quantiles: signed Spearman quantiles abs_spearman_quantiles: absolute Spearman quantiles genie3_histogram: histogram information for GENIE3 scores spearman_histogram: histogram information for Spearman rho
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
inspect_integrated_grn_score_distribution <- function(
    edges,
    edge_type,
    genie3_breaks = 100,
    spearman_breaks = 100
) {
  
  required_columns <- c("genie3_weight", "spearman_rho")
  
  if (!all(required_columns %in% colnames(edges))) {
    stop("Required GENIE3 or Spearman columns are missing from the edge table")
  }
  
  probabilities <- c(
    0,
    0.50,
    0.75,
    0.90,
    0.95,
    0.975,
    0.99,
    0.995,
    0.999,
    1
  )
  
  genie3_quantiles <- stats::quantile(
    edges$genie3_weight,
    probs = probabilities,
    na.rm = TRUE,
    names = TRUE
  )
  
  spearman_quantiles <- stats::quantile(
    edges$spearman_rho,
    probs = probabilities,
    na.rm = TRUE,
    names = TRUE
  )
  
  abs_spearman_quantiles <- stats::quantile(
    abs(edges$spearman_rho),
    probs = probabilities,
    na.rm = TRUE,
    names = TRUE
  )
  
  summary <- data.frame(
    edge_type = edge_type,
    number_of_edges = nrow(edges),
    zero_genie3_weights = sum(edges$genie3_weight == 0, na.rm = TRUE),
    missing_genie3_weights = sum(is.na(edges$genie3_weight)),
    missing_spearman_values = sum(is.na(edges$spearman_rho)),
    stringsAsFactors = FALSE
  )
  
  genie3_histogram <- graphics::hist(
    edges$genie3_weight,
    breaks = genie3_breaks,
    plot = FALSE
  )
  
  spearman_histogram <- graphics::hist(
    edges$spearman_rho,
    breaks = seq(-1, 1, length.out = spearman_breaks + 1),
    plot = FALSE
  )
  
  return(list(
    summary = summary,
    genie3_quantiles = genie3_quantiles,
    spearman_quantiles = spearman_quantiles,
    abs_spearman_quantiles = abs_spearman_quantiles,
    genie3_histogram = genie3_histogram,
    spearman_histogram = spearman_histogram
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: evaluate_integrated_grn_filter_thresholds
# ------------------------------------------------------------------------------

# Evaluate candidate filtering thresholds for one integrated GRN interaction
# class.
#
# This function:
# - validates the required GENIE3 and Spearman columns
# - evaluates multiple class-specific GENIE3 quantile thresholds
# - evaluates multiple absolute Spearman correlation thresholds
# - counts the number of interactions retained by each threshold combination
# - calculates the fraction of retained interactions relative to the original
#   interaction class
# - does not modify or filter the input edge table
#
# INPUT:
# - edges: interaction table generated by infer_integrated_grn_relationships()
# - edge_type: label identifying the interaction class
# - genie3_quantiles: GENIE3 quantiles to evaluate as class-specific thresholds
# - spearman_thresholds: absolute Spearman correlation thresholds to evaluate
#
# OUTPUT:
# - data frame containing one row for each GENIE3/Spearman threshold combination,
#   with:
#   - edge_type
#   - genie3_quantile
#   - genie3_threshold
#   - spearman_threshold
#   - retained_edges
#   - retained_fraction
# ------------------------------------------------------------------------------

#' Evaluate candidate filtering thresholds for one integrated GRN interaction
#'
#' class.
#' This function:
#' - validates the required GENIE3 and Spearman columns
#' - evaluates multiple class-specific GENIE3 quantile thresholds
#' - evaluates multiple absolute Spearman correlation thresholds
#' - counts the number of interactions retained by each threshold combination
#' - calculates the fraction of retained interactions relative to the original
#' interaction class
#' - does not modify or filter the input edge table
#'
#' @param edges interaction table generated by infer_integrated_grn_relationships()
#' @param edge_type label identifying the interaction class
#' @param genie3_quantiles GENIE3 quantiles to evaluate as class-specific thresholds
#' @param spearman_thresholds absolute Spearman correlation thresholds to evaluate
#'
#' @return data frame containing one row for each GENIE3/Spearman threshold combination, with: edge_type genie3_quantile genie3_threshold spearman_threshold retained_edges retained_fraction
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
evaluate_integrated_grn_filter_thresholds <- function(
    edges,
    edge_type,
    genie3_quantiles = c(0.90, 0.95, 0.975, 0.99),
    spearman_thresholds = c(0.3, 0.4, 0.5)
) {
  
  required_columns <- c("genie3_weight", "spearman_rho")
  
  if (!all(required_columns %in% colnames(edges))) {
    stop("Required GENIE3 or Spearman columns are missing from the edge table")
  }
  
  results <- list()
  result_index <- 1L
  
  for (genie3_quantile in genie3_quantiles) {
    
    genie3_threshold <- stats::quantile(
      edges$genie3_weight,
      probs = genie3_quantile,
      na.rm = TRUE,
      names = FALSE
    )
    
    for (spearman_threshold in spearman_thresholds) {
      
      keep <- edges$genie3_weight >= genie3_threshold &
        abs(edges$spearman_rho) >= spearman_threshold
      
      retained_edges <- sum(keep, na.rm = TRUE)
      
      results[[result_index]] <- data.frame(
        edge_type = edge_type,
        genie3_quantile = genie3_quantile,
        genie3_threshold = genie3_threshold,
        spearman_threshold = spearman_threshold,
        retained_edges = retained_edges,
        retained_fraction = retained_edges / nrow(edges),
        stringsAsFactors = FALSE
      )
      
      result_index <- result_index + 1L
    }
  }
  
  do.call(rbind, results)
}

# ------------------------------------------------------------------------------
# FUNCTION: filter_integrated_grn_relationships
# ------------------------------------------------------------------------------

# Filter one integrated GRN interaction class using a class-specific GENIE3
# quantile threshold and an absolute Spearman correlation threshold.
#
# This function:
# - validates the required GENIE3 and Spearman columns
# - calculates the GENIE3 threshold from the specified quantile
# - retains interactions with GENIE3 importance at or above the threshold
# - retains interactions with absolute Spearman correlation at or above the
#   specified threshold
# - preserves all original edge annotations
# - reports the applied thresholds and the number of retained interactions
#
# INPUT:
# - edges: interaction table generated by infer_integrated_grn_relationships()
# - genie3_quantile: class-specific GENIE3 quantile used for filtering
# - spearman_threshold: minimum absolute Spearman correlation
#
# OUTPUT:
# - list containing:
#   - edges: filtered interaction table
#   - summary: filtering summary and applied thresholds
# ------------------------------------------------------------------------------

#' Filter one integrated GRN interaction class using a class-specific GENIE3
#'
#' quantile threshold and an absolute Spearman correlation threshold.
#' This function:
#' - validates the required GENIE3 and Spearman columns
#' - calculates the GENIE3 threshold from the specified quantile
#' - retains interactions with GENIE3 importance at or above the threshold
#' - retains interactions with absolute Spearman correlation at or above the
#' specified threshold
#' - preserves all original edge annotations
#' - reports the applied thresholds and the number of retained interactions
#'
#' @param edges interaction table generated by infer_integrated_grn_relationships()
#' @param genie3_quantile class-specific GENIE3 quantile used for filtering
#' @param spearman_threshold minimum absolute Spearman correlation
#'
#' @return list containing: edges: filtered interaction table summary: filtering summary and applied thresholds
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
filter_integrated_grn_relationships <- function(
    edges,
    genie3_quantile = 0.95,
    spearman_threshold = 0.5
) {
  
  required_columns <- c("genie3_weight", "spearman_rho", "edge_type")
  
  if (!all(required_columns %in% colnames(edges))) {
    stop("Required GRN columns are missing from the edge table")
  }
  
  genie3_threshold <- stats::quantile(
    edges$genie3_weight,
    probs = genie3_quantile,
    na.rm = TRUE,
    names = FALSE
  )
  
  keep <- edges$genie3_weight >= genie3_threshold &
    abs(edges$spearman_rho) >= spearman_threshold
  
  filtered_edges <- edges[keep & !is.na(keep), , drop = FALSE]
  
  summary <- data.frame(
    edge_type = unique(edges$edge_type),
    original_edges = nrow(edges),
    genie3_quantile = genie3_quantile,
    genie3_threshold = genie3_threshold,
    spearman_threshold = spearman_threshold,
    retained_edges = nrow(filtered_edges),
    retained_fraction = nrow(filtered_edges) / nrow(edges),
    stringsAsFactors = FALSE
  )
  
  return(list(
    edges = filtered_edges,
    summary = summary
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: integrate_filtered_grn_interaction_classes
# ------------------------------------------------------------------------------

# Integrate independently filtered GRN interaction classes into a single
# regulatory network.
#
# This function:
# - validates the required columns in each filtered interaction table
# - combines all filtered interaction classes
# - checks for missing regulator, target, score, and annotation values
# - checks for duplicated regulator-target-edge_type combinations
# - constructs a node table from all regulators and targets
# - assigns each node to the lncRNA, TF, or PCG gene class
# - calculates in-degree, out-degree, and total degree in the integrated network
# - preserves all filtered edges without applying additional selection
#
# INPUT:
# - edge_tables: named list of filtered interaction tables
# - gene_sets: integrated GRN gene sets containing lncRNA, TF, and PCG
#
# OUTPUT:
# - list containing:
#   - edges: integrated filtered interaction table
#   - nodes: integrated node annotation table
#   - summary: integration and validation summary
# ------------------------------------------------------------------------------

#' Integrate independently filtered GRN interaction classes into a single
#'
#' regulatory network.
#' This function:
#' - validates the required columns in each filtered interaction table
#' - combines all filtered interaction classes
#' - checks for missing regulator, target, score, and annotation values
#' - checks for duplicated regulator-target-edge_type combinations
#' - constructs a node table from all regulators and targets
#' - assigns each node to the lncRNA, TF, or PCG gene class
#' - calculates in-degree, out-degree, and total degree in the integrated network
#' - preserves all filtered edges without applying additional selection
#'
#' @param edge_tables named list of filtered interaction tables
#' @param gene_sets integrated GRN gene sets containing lncRNA, TF, and PCG
#'
#' @return list containing: edges: integrated filtered interaction table nodes: integrated node annotation table summary: integration and validation summary
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
integrate_filtered_grn_interaction_classes <- function(
    edge_tables,
    gene_sets
) {
  
  required_columns <- c(
    "regulator",
    "target",
    "genie3_weight",
    "spearman_rho",
    "regulator_type",
    "target_type",
    "edge_type"
  )
  
  if (!is.list(edge_tables) || length(edge_tables) == 0L) {
    stop("edge_tables must be a non-empty list")
  }
  
  for (i in seq_along(edge_tables)) {
    
    if (!all(required_columns %in% colnames(edge_tables[[i]]))) {
      stop(
        paste0(
          "Required columns are missing from interaction table: ",
          names(edge_tables)[i]
        )
      )
    }
  }
  
  integrated_edges <- do.call(
    rbind,
    lapply(
      edge_tables,
      function(x) x[, required_columns, drop = FALSE]
    )
  )
  
  rownames(integrated_edges) <- NULL
  
  duplicated_edges <- duplicated(
    integrated_edges[, c(
      "regulator",
      "target",
      "edge_type"
    )]
  )
  
  node_names <- unique(c(
    integrated_edges$regulator,
    integrated_edges$target
  ))
  
  node_type <- ifelse(
    node_names %in% gene_sets$lncRNA,
    "lncRNA",
    ifelse(
      node_names %in% gene_sets$TF,
      "TF",
      ifelse(
        node_names %in% gene_sets$PCG,
        "PCG",
        NA_character_
      )
    )
  )
  
  in_degree <- table(integrated_edges$target)
  out_degree <- table(integrated_edges$regulator)
  
  nodes <- data.frame(
    gene = node_names,
    node_type = node_type,
    in_degree = as.integer(in_degree[node_names]),
    out_degree = as.integer(out_degree[node_names]),
    stringsAsFactors = FALSE
  )
  
  nodes$in_degree[is.na(nodes$in_degree)] <- 0L
  nodes$out_degree[is.na(nodes$out_degree)] <- 0L
  nodes$total_degree <- nodes$in_degree + nodes$out_degree
  
  summary <- data.frame(
    metric = c(
      "number_of_edges",
      "number_of_nodes",
      "duplicated_edges",
      "missing_regulators",
      "missing_targets",
      "missing_genie3_weights",
      "missing_spearman_values",
      "missing_regulator_types",
      "missing_target_types",
      "missing_edge_types",
      "unclassified_nodes"
    ),
    value = c(
      nrow(integrated_edges),
      nrow(nodes),
      sum(duplicated_edges),
      sum(is.na(integrated_edges$regulator)),
      sum(is.na(integrated_edges$target)),
      sum(is.na(integrated_edges$genie3_weight)),
      sum(is.na(integrated_edges$spearman_rho)),
      sum(is.na(integrated_edges$regulator_type)),
      sum(is.na(integrated_edges$target_type)),
      sum(is.na(integrated_edges$edge_type)),
      sum(is.na(nodes$node_type))
    ),
    stringsAsFactors = FALSE
  )
  
  return(list(
    edges = integrated_edges,
    nodes = nodes,
    summary = summary
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: annotate_integrated_grn_with_lnctard
# ------------------------------------------------------------------------------

# Annotate integrated GRN lncRNA-target interactions using LncTarD.
#
# This function:
# - validates the required integrated GRN and LncTarD columns
# - selects LncTarD interactions with a lncRNA regulator
# - matches directed regulator-target pairs between the GRN and LncTarD
# - annotates lncRNA-to-PCG and lncRNA-to-TF GRN edges
# - counts the number of supporting LncTarD records for each interaction
# - stores the corresponding LncTarD record identifiers
# - preserves all integrated GRN edges regardless of database support
# - does not modify GENIE3 or Spearman scores
#
# INPUT:
# - integrated_edges: filtered and integrated GRN edge table
# - lnctard: LncTarD interaction table
#
# OUTPUT:
# - list containing:
#   - edges: complete integrated GRN edge table with LncTarD annotations
#   - lnc_edges: lncRNA-regulated GRN interactions with LncTarD annotations
#   - supported_edges: lncRNA-regulated interactions supported by LncTarD
#   - summary: LncTarD annotation summary
# ------------------------------------------------------------------------------

#' Annotate integrated GRN lncRNA-target interactions using LncTarD
#'
#' This function:
#' - validates the required integrated GRN and LncTarD columns
#' - selects LncTarD interactions with a lncRNA regulator
#' - matches directed regulator-target pairs between the GRN and LncTarD
#' - annotates lncRNA-to-PCG and lncRNA-to-TF GRN edges
#' - counts the number of supporting LncTarD records for each interaction
#' - stores the corresponding LncTarD record identifiers
#' - preserves all integrated GRN edges regardless of database support
#' - does not modify GENIE3 or Spearman scores
#'
#' @param integrated_edges filtered and integrated GRN edge table
#' @param lnctard LncTarD interaction table
#'
#' @return list containing: edges: complete integrated GRN edge table with LncTarD annotations lnc_edges: lncRNA-regulated GRN interactions with LncTarD annotations supported_edges: lncRNA-regulated interactions supported by LncTarD summary: LncTarD annotation summary
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
annotate_integrated_grn_with_lnctard <- function(
    integrated_edges,
    lnctard
) {
  
  required_grn_columns <- c(
    "regulator",
    "target",
    "regulator_type",
    "target_type",
    "edge_type"
  )
  
  required_lnctard_columns <- c(
    "Regulator",
    "Target",
    "RegulatorType",
    "TargetType",
    "RID"
  )
  
  if (!all(required_grn_columns %in% colnames(integrated_edges))) {
    stop("Required integrated GRN columns are missing")
  }
  
  if (!all(required_lnctard_columns %in% colnames(lnctard))) {
    stop("Required LncTarD columns are missing")
  }
  
  lnctard_lnc <- lnctard[
    tolower(trimws(lnctard$RegulatorType)) == "lncrna",
    ,
    drop = FALSE
  ]
  
  lnctard_lnc$regulator_key <- toupper(trimws(lnctard_lnc$Regulator))
  lnctard_lnc$target_key <- toupper(trimws(lnctard_lnc$Target))
  
  lnctard_lnc$pair_key <- paste(
    lnctard_lnc$regulator_key,
    lnctard_lnc$target_key,
    sep = "||"
  )
  
  lnctard_pair_records <- split(
    lnctard_lnc,
    lnctard_lnc$pair_key
  )
  
  lnctard_pair_summary <- do.call(
    rbind,
    lapply(
      names(lnctard_pair_records),
      function(pair_key) {
        
        x <- lnctard_pair_records[[pair_key]]
        
        data.frame(
          pair_key = pair_key,
          lnctard_records = nrow(x),
          lnctard_rids = paste(unique(x$RID), collapse = ";"),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  annotated_edges <- integrated_edges
  
  annotated_edges$pair_key <- paste(
    toupper(trimws(annotated_edges$regulator)),
    toupper(trimws(annotated_edges$target)),
    sep = "||"
  )
  
  annotated_edges$lnctard_supported <- FALSE
  annotated_edges$lnctard_records <- 0L
  annotated_edges$lnctard_rids <- NA_character_
  
  lnc_edge_index <- annotated_edges$regulator_type == "lncRNA"
  
  match_index <- match(
    annotated_edges$pair_key[lnc_edge_index],
    lnctard_pair_summary$pair_key
  )
  
  supported_index <- !is.na(match_index)
  
  lnc_rows <- which(lnc_edge_index)
  supported_rows <- lnc_rows[supported_index]
  
  annotated_edges$lnctard_supported[supported_rows] <- TRUE
  
  annotated_edges$lnctard_records[supported_rows] <-
    lnctard_pair_summary$lnctard_records[
      match_index[supported_index]
    ]
  
  annotated_edges$lnctard_rids[supported_rows] <-
    lnctard_pair_summary$lnctard_rids[
      match_index[supported_index]
    ]
  
  annotated_edges$pair_key <- NULL
  
  lnc_edges <- annotated_edges[
    annotated_edges$regulator_type == "lncRNA",
    ,
    drop = FALSE
  ]
  
  supported_edges <- lnc_edges[
    lnc_edges$lnctard_supported,
    ,
    drop = FALSE
  ]
  
  summary <- data.frame(
    metric = c(
      "integrated_grn_edges",
      "lncRNA_regulated_edges",
      "lncRNA_to_PCG_edges",
      "lncRNA_to_TF_edges",
      "lnctard_supported_lncRNA_edges",
      "lnctard_unsupported_lncRNA_edges",
      "lnctard_supported_lncRNA_to_PCG",
      "lnctard_supported_lncRNA_to_TF"
    ),
    value = c(
      nrow(annotated_edges),
      nrow(lnc_edges),
      sum(lnc_edges$edge_type == "lncRNA_to_PCG"),
      sum(lnc_edges$edge_type == "lncRNA_to_TF"),
      sum(lnc_edges$lnctard_supported),
      sum(!lnc_edges$lnctard_supported),
      sum(
        lnc_edges$edge_type == "lncRNA_to_PCG" &
          lnc_edges$lnctard_supported
      ),
      sum(
        lnc_edges$edge_type == "lncRNA_to_TF" &
          lnc_edges$lnctard_supported
      )
    ),
    stringsAsFactors = FALSE
  )
  
  return(list(
    edges = annotated_edges,
    lnc_edges = lnc_edges,
    supported_edges = supported_edges,
    summary = summary
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: analyze_integrated_grn_structure
# ------------------------------------------------------------------------------

# Analyze the global structure of the integrated GRN.
#
# This function:
# - calculates the total number of nodes and edges
# - summarizes nodes by gene class
# - summarizes edges by interaction class
# - distinguishes positive and negative Spearman associations
# - calculates in-degree, out-degree, and total degree for every node
# - identifies the most highly connected nodes
# - summarizes LncTarD-supported interactions when available
#
# INPUT:
# - edges: integrated GRN edge table
# - nodes: integrated GRN node table
# - top_n: number of highest-degree nodes to return
#
# OUTPUT:
# - list containing:
#   - network_summary: global network statistics
#   - node_type_summary: number of nodes by gene class
#   - edge_type_summary: number of edges by interaction class
#   - spearman_sign_summary: positive/negative edge counts
#   - degree_table: node-level degree statistics
#   - top_degree_nodes: highest total-degree nodes
#   - top_in_degree_nodes: highest in-degree nodes
#   - top_out_degree_nodes: highest out-degree nodes
#   - lnctard_summary: LncTarD support statistics
# ------------------------------------------------------------------------------

#' Analyze the global structure of the integrated GRN
#'
#' This function:
#' - calculates the total number of nodes and edges
#' - summarizes nodes by gene class
#' - summarizes edges by interaction class
#' - distinguishes positive and negative Spearman associations
#' - calculates in-degree, out-degree, and total degree for every node
#' - identifies the most highly connected nodes
#' - summarizes LncTarD-supported interactions when available
#'
#' @param edges integrated GRN edge table
#' @param nodes integrated GRN node table
#' @param top_n number of highest-degree nodes to return
#'
#' @return list containing: network_summary: global network statistics node_type_summary: number of nodes by gene class edge_type_summary: number of edges by interaction class spearman_sign_summary: positive/negative edge counts degree_table: node-level degree statistics top_degree_nodes: highest total-degree nodes top_in_degree_nodes: highest in-degree nodes top_out_degree_nodes: highest out-degree nodes lnctard_summary: LncTarD support statistics
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
analyze_integrated_grn_structure <- function(
    edges,
    nodes,
    top_n = 20
) {
  
  required_edge_columns <- c(
    "regulator",
    "target",
    "spearman_rho",
    "regulator_type",
    "target_type",
    "edge_type"
  )
  
  if (!all(required_edge_columns %in% colnames(edges))) {
    stop("Required integrated GRN edge columns are missing")
  }
  
  required_node_columns <- c(
    "gene",
    "node_type"
  )
  
  if (!all(required_node_columns %in% colnames(nodes))) {
    stop("Required integrated GRN node columns are missing")
  }
  
  in_degree <- table(edges$target)
  out_degree <- table(edges$regulator)
  
  degree_table <- data.frame(
    gene = nodes$gene,
    node_type = nodes$node_type,
    stringsAsFactors = FALSE
  )
  
  degree_table$in_degree <- as.integer(
    in_degree[degree_table$gene]
  )
  
  degree_table$out_degree <- as.integer(
    out_degree[degree_table$gene]
  )
  
  degree_table$in_degree[is.na(degree_table$in_degree)] <- 0L
  degree_table$out_degree[is.na(degree_table$out_degree)] <- 0L
  
  degree_table$total_degree <-
    degree_table$in_degree +
    degree_table$out_degree
  
  network_summary <- data.frame(
    metric = c(
      "number_of_nodes",
      "number_of_edges",
      "positive_spearman_edges",
      "negative_spearman_edges",
      "zero_spearman_edges"
    ),
    value = c(
      nrow(nodes),
      nrow(edges),
      sum(edges$spearman_rho > 0, na.rm = TRUE),
      sum(edges$spearman_rho < 0, na.rm = TRUE),
      sum(edges$spearman_rho == 0, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  node_type_summary <- as.data.frame(
    table(nodes$node_type),
    stringsAsFactors = FALSE
  )
  
  colnames(node_type_summary) <- c(
    "node_type",
    "number_of_nodes"
  )
  
  edge_type_summary <- as.data.frame(
    table(edges$edge_type),
    stringsAsFactors = FALSE
  )
  
  colnames(edge_type_summary) <- c(
    "edge_type",
    "number_of_edges"
  )
  
  spearman_sign <- ifelse(
    edges$spearman_rho > 0,
    "positive",
    ifelse(
      edges$spearman_rho < 0,
      "negative",
      "zero"
    )
  )
  
  spearman_sign_summary <- as.data.frame(
    table(spearman_sign),
    stringsAsFactors = FALSE
  )
  
  colnames(spearman_sign_summary) <- c(
    "spearman_sign",
    "number_of_edges"
  )
  
  top_degree_nodes <- head(
    degree_table[
      order(
        degree_table$total_degree,
        decreasing = TRUE
      ),
    ],
    top_n
  )
  
  top_in_degree_nodes <- head(
    degree_table[
      order(
        degree_table$in_degree,
        decreasing = TRUE
      ),
    ],
    top_n
  )
  
  top_out_degree_nodes <- head(
    degree_table[
      order(
        degree_table$out_degree,
        decreasing = TRUE
      ),
    ],
    top_n
  )
  
  if ("lnctard_supported" %in% colnames(edges)) {
    
    lnctard_summary <- data.frame(
      metric = c(
        "lnctard_supported_lncRNA_edges",
        "lnctard_unsupported_lncRNA_edges"
      ),
      value = c(
        sum(edges$lnctard_supported, na.rm = TRUE),
        sum(
          edges$regulator_type == "lncRNA" &
            !edges$lnctard_supported,
          na.rm = TRUE
        )
      ),
      stringsAsFactors = FALSE
    )
    
  } else {
    
    lnctard_summary <- NULL
  }
  
  return(list(
    network_summary = network_summary,
    node_type_summary = node_type_summary,
    edge_type_summary = edge_type_summary,
    spearman_sign_summary = spearman_sign_summary,
    degree_table = degree_table,
    top_degree_nodes = top_degree_nodes,
    top_in_degree_nodes = top_in_degree_nodes,
    top_out_degree_nodes = top_out_degree_nodes,
    lnctard_summary = lnctard_summary
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_integrated_grn_cascades
# ------------------------------------------------------------------------------

# Summarize directed two-step cascades in the integrated GRN without
# materializing every individual path.
#
# This function calculates the number of possible:
#
# 1. lncRNA -> TF -> PCG cascades
# 2. TF -> lncRNA -> PCG cascades
#
# Cascade counts are obtained from the number of compatible incoming and
# outgoing retained edges for each intermediate node.
#
# This avoids creating potentially very large global path tables.
#
# INPUT:
# - edges: integrated GRN edge table
#
# OUTPUT:
# - list containing:
#   - summary: total number of cascades by type
#   - lncRNA_TF_PCG_by_TF: cascade counts for each intermediate TF
#   - TF_lncRNA_PCG_by_lncRNA: cascade counts for each intermediate lncRNA
# ------------------------------------------------------------------------------

#' Summarize directed two-step cascades in the integrated GRN without
#'
#' materializing every individual path.
#' This function calculates the number of possible:
#' 1. lncRNA -> TF -> PCG cascades
#' 2. TF -> lncRNA -> PCG cascades
#' Cascade counts are obtained from the number of compatible incoming and
#' outgoing retained edges for each intermediate node.
#' This avoids creating potentially very large global path tables.
#'
#' @param edges integrated GRN edge table
#'
#' @return list containing: summary: total number of cascades by type lncRNA_TF_PCG_by_TF: cascade counts for each intermediate TF TF_lncRNA_PCG_by_lncRNA: cascade counts for each intermediate lncRNA
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
summarize_integrated_grn_cascades <- function(edges) {
  
  required_columns <- c(
    "regulator",
    "target",
    "edge_type"
  )
  
  if (!all(required_columns %in% colnames(edges))) {
    stop("Required integrated GRN columns are missing")
  }
  
  lnc_tf <- edges[
    edges$edge_type == "lncRNA_to_TF",
    c("regulator", "target"),
    drop = FALSE
  ]
  
  tf_pcg <- edges[
    edges$edge_type == "TF_to_PCG",
    c("regulator", "target"),
    drop = FALSE
  ]
  
  tf_lnc <- edges[
    edges$edge_type == "TF_to_lncRNA",
    c("regulator", "target"),
    drop = FALSE
  ]
  
  lnc_pcg <- edges[
    edges$edge_type == "lncRNA_to_PCG",
    c("regulator", "target"),
    drop = FALSE
  ]
  
  lnc_to_tf_counts <- table(lnc_tf$target)
  tf_to_pcg_counts <- table(tf_pcg$regulator)
  
  shared_tf <- intersect(
    names(lnc_to_tf_counts),
    names(tf_to_pcg_counts)
  )
  
  lncRNA_TF_PCG_by_TF <- data.frame(
    TF = shared_tf,
    incoming_lncRNA_edges = as.integer(lnc_to_tf_counts[shared_tf]),
    outgoing_PCG_edges = as.integer(tf_to_pcg_counts[shared_tf]),
    stringsAsFactors = FALSE
  )
  
  lncRNA_TF_PCG_by_TF$number_of_cascades <-
    lncRNA_TF_PCG_by_TF$incoming_lncRNA_edges *
    lncRNA_TF_PCG_by_TF$outgoing_PCG_edges
  
  tf_to_lnc_counts <- table(tf_lnc$target)
  lnc_to_pcg_counts <- table(lnc_pcg$regulator)
  
  shared_lnc <- intersect(
    names(tf_to_lnc_counts),
    names(lnc_to_pcg_counts)
  )
  
  TF_lncRNA_PCG_by_lncRNA <- data.frame(
    lncRNA = shared_lnc,
    incoming_TF_edges = as.integer(tf_to_lnc_counts[shared_lnc]),
    outgoing_PCG_edges = as.integer(lnc_to_pcg_counts[shared_lnc]),
    stringsAsFactors = FALSE
  )
  
  TF_lncRNA_PCG_by_lncRNA$number_of_cascades <-
    TF_lncRNA_PCG_by_lncRNA$incoming_TF_edges *
    TF_lncRNA_PCG_by_lncRNA$outgoing_PCG_edges
  
  summary <- data.frame(
    cascade_type = c(
      "lncRNA_to_TF_to_PCG",
      "TF_to_lncRNA_to_PCG"
    ),
    number_of_cascades = c(
      sum(lncRNA_TF_PCG_by_TF$number_of_cascades),
      sum(TF_lncRNA_PCG_by_lncRNA$number_of_cascades)
    ),
    stringsAsFactors = FALSE
  )
  
  return(list(
    summary = summary,
    lncRNA_TF_PCG_by_TF = lncRNA_TF_PCG_by_TF,
    TF_lncRNA_PCG_by_lncRNA = TF_lncRNA_PCG_by_lncRNA
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: build_integrated_grn_graph
# ------------------------------------------------------------------------------

# Build a directed igraph object from the integrated GRN.
#
# This function:
# - converts the complete integrated GRN into a directed graph
# - preserves edge attributes including GENIE3, Spearman, edge type, and
#   LncTarD annotations
# - preserves node attributes including gene class and degree information
# - does not filter or modify the integrated network
#
# The resulting graph can be reused for:
# - global network visualization
# - inbound and outbound neighborhood exploration
# - user-defined depth analysis
# - gene-centered path analysis
# - cycle and motif exploration
#
# INPUT:
# - edges: integrated GRN edge table
# - nodes: integrated GRN node table
#
# OUTPUT:
# - directed igraph object representing the complete integrated GRN
# ------------------------------------------------------------------------------

#' Build a directed igraph object from the integrated GRN
#'
#' This function:
#' - converts the complete integrated GRN into a directed graph
#' - preserves edge attributes including GENIE3, Spearman, edge type, and
#' LncTarD annotations
#' - preserves node attributes including gene class and degree information
#' - does not filter or modify the integrated network
#' The resulting graph can be reused for:
#' - global network visualization
#' - inbound and outbound neighborhood exploration
#' - user-defined depth analysis
#' - gene-centered path analysis
#' - cycle and motif exploration
#'
#' @param edges integrated GRN edge table
#' @param nodes integrated GRN node table
#'
#' @return directed igraph object representing the complete integrated GRN
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
build_integrated_grn_graph <- function(
    edges,
    nodes
) {
  
  required_edge_columns <- c(
    "regulator",
    "target"
  )
  
  required_node_columns <- c(
    "gene",
    "node_type"
  )
  
  if (!all(required_edge_columns %in% colnames(edges))) {
    stop("Required integrated GRN edge columns are missing")
  }
  
  if (!all(required_node_columns %in% colnames(nodes))) {
    stop("Required integrated GRN node columns are missing")
  }
  
  graph_edges <- edges
  
  graph_nodes <- nodes
  
  graph_nodes <- graph_nodes[
    ,
    c(
      "gene",
      setdiff(colnames(graph_nodes), "gene")
    ),
    drop = FALSE
  ]
  
  grn_graph <- igraph::graph_from_data_frame(
    d = graph_edges,
    directed = TRUE,
    vertices = graph_nodes
  )
  
  return(grn_graph)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_integrated_grn
# ------------------------------------------------------------------------------

# Visualize the complete integrated GRN without applying additional filtering.
#
# This function:
# - validates that the input object is an igraph network
# - uses all nodes and edges in the integrated GRN
# - arranges lncRNA, TF, and PCG nodes in separate vertical layers
# - distributes nodes across the available plotting area
# - uses transparent edges to reduce overplotting
# - distinguishes node and interaction classes
# - optionally displays legends
# - does not remove nodes or edges for visualization purposes
#
# INPUT:
# - graph: directed igraph object representing the integrated GRN
# - layout: optional precomputed layout matrix
# - show_legend: whether to display node and edge class legends
# - vertex_size: graphical size of network nodes
# - edge_width: graphical width of network edges
# - edge_alpha: transparency applied to network edges
# - arrow_size: graphical size of directed edge arrows
#
# OUTPUT:
# - list containing:
#   - layout: layout matrix used for visualization
#   - node_colors: colors assigned to network nodes
#   - edge_colors: colors assigned to network edges
# ------------------------------------------------------------------------------

#' Visualize the complete integrated GRN without applying additional filtering
#'
#' This function:
#' - validates that the input object is an igraph network
#' - uses all nodes and edges in the integrated GRN
#' - arranges lncRNA, TF, and PCG nodes in separate vertical layers
#' - distributes nodes across the available plotting area
#' - uses transparent edges to reduce overplotting
#' - distinguishes node and interaction classes
#' - optionally displays legends
#' - does not remove nodes or edges for visualization purposes
#'
#' @param graph directed igraph object representing the integrated GRN
#' @param layout optional precomputed layout matrix
#' @param show_legend whether to display node and edge class legends
#' @param vertex_size graphical size of network nodes
#' @param edge_width graphical width of network edges
#' @param edge_alpha transparency applied to network edges
#' @param arrow_size graphical size of directed edge arrows
#'
#' @return list containing: layout: layout matrix used for visualization node_colors: colors assigned to network nodes edge_colors: colors assigned to network edges
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_integrated_grn <- function(
    graph,
    layout = NULL,
    show_legend = TRUE,
    vertex_size = 0.8,
    edge_width = 0.08,
    edge_alpha = 0.04,
    arrow_size = 0
) {
  
  if (!inherits(graph, "igraph")) {
    stop("graph must be an igraph object")
  }
  
  node_type <- igraph::vertex_attr(graph, "node_type")
  edge_type <- igraph::edge_attr(graph, "edge_type")
  
  if (is.null(layout)) {
    
    layout <- matrix(
      NA_real_,
      nrow = igraph::vcount(graph),
      ncol = 2
    )
    
    lnc_index <- which(node_type == "lncRNA")
    tf_index <- which(node_type == "TF")
    pcg_index <- which(node_type == "PCG")
    
    layout[lnc_index, 1] <- -1
    layout[tf_index, 1] <- 0
    layout[pcg_index, 1] <- 1
    
    layout[lnc_index, 2] <- seq(
      -1,
      1,
      length.out = length(lnc_index)
    )
    
    layout[tf_index, 2] <- seq(
      -1,
      1,
      length.out = length(tf_index)
    )
    
    layout[pcg_index, 2] <- seq(
      -1,
      1,
      length.out = length(pcg_index)
    )
  }
  
  node_palette <- c(
    lncRNA = "#D55E00",
    TF = "#0072B2",
    PCG = "#999999"
  )
  
  edge_palette <- c(
    lncRNA_to_PCG = "#E69F00",
    lncRNA_to_TF = "#CC79A7",
    TF_to_lncRNA = "#56B4E9",
    TF_to_PCG = "#666666"
  )
  
  node_colors <- unname(node_palette[node_type])
  
  edge_colors <- grDevices::adjustcolor(
    unname(edge_palette[edge_type]),
    alpha.f = edge_alpha
  )
  
  graphics::par(
    mar = c(1, 1, 2, 1),
    xpd = NA
  )
  
  igraph::plot.igraph(
    graph,
    layout = layout,
    rescale = FALSE,
    xlim = c(-1.35, 1.35),
    ylim = c(-1.05, 1.05),
    vertex.label = NA,
    vertex.size = vertex_size,
    vertex.color = node_colors,
    vertex.frame.color = NA,
    edge.color = edge_colors,
    edge.width = edge_width,
    edge.arrow.size = arrow_size,
    edge.curved = 0
  )
  
  graphics::text(
    x = c(-1, 0, 1),
    y = c(1.03, 1.03, 1.03),
    labels = c("lncRNA", "TF", "PCG"),
    font = 2
  )
  
  if (show_legend) {
    
    graphics::legend(
      "bottom",
      inset = c(0, -0.12),
      legend = names(edge_palette),
      col = edge_palette,
      lty = 1,
      lwd = 2,
      horiz = TRUE,
      bty = "n",
      xpd = NA,
      title = "Edge type"
    )
  }
  
  invisible(
    list(
      layout = layout,
      node_colors = node_colors,
      edge_colors = edge_colors
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: explore_integrated_grn_gene
# ------------------------------------------------------------------------------

# Explore the neighborhood and directed paths of a selected gene in the
# integrated GRN.
#
# This function:
# - validates the selected gene
# - supports outbound, inbound, or bidirectional exploration
# - identifies all nodes reachable within a user-defined depth
# - assigns each reached node to its minimum network distance from the
#   selected gene
# - extracts the complete induced subnetwork containing the reached nodes
# - extracts a depth-layer network containing only edges connecting
#   consecutive exploration levels
# - preserves node and edge attributes from the integrated GRN
# - optionally enumerates simple paths up to the selected maximum depth
# - summarizes the number of nodes reached at each depth
# - summarizes node and edge composition in the selected neighborhood
#
# INPUT:
# - graph: directed igraph object representing the integrated GRN
# - gene: selected gene name
# - direction: exploration direction; one of "out", "in", or "all"
# - depth: maximum network depth
# - enumerate_paths: whether to enumerate simple paths
# - max_paths: maximum number of paths to retain in the returned path table
#
# OUTPUT:
# - list containing:
#   - gene: selected gene
#   - direction: selected exploration direction
#   - depth: selected maximum depth
#   - summary: gene-centered neighborhood summary
#   - depth_summary: number of reached nodes by minimum depth
#   - node_type_summary: reached nodes by gene class
#   - edge_type_summary: induced-subnetwork edges by interaction class
#   - layer_edge_type_summary: depth-layer edges by interaction class
#   - nodes: node table for the selected neighborhood
#   - edges: complete induced-subnetwork edge table
#   - layer_edges: edges connecting consecutive depth levels
#   - paths: optional simple-path table
#   - paths_truncated: whether the returned path table reached max_paths
#   - subgraph: complete induced igraph subnetwork
#   - layer_graph: igraph representation of the depth-layer network
# ------------------------------------------------------------------------------

#' Explore integrated grn gene
#'
#' Explore integrated grn gene.
#'
#' @param graph Integrated GRN graph generated by `build_integrated_grn_graph()`.
#' @param gene Single gene symbol to explore.
#' @param direction Neighborhood direction to explore.
#' @param depth Maximum graph distance used for the neighborhood search.
#' @param enumerate_paths Logical indicating whether simple paths should be enumerated.
#' @param max_paths Maximum number of paths returned when path enumeration is enabled.
#'
#' @return A list containing summaries, edge tables, optional path tables, and graph objects for the selected neighborhood.
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
explore_integrated_grn_gene <- function(
    graph,
    gene,
    direction = c("out", "in", "all"),
    depth = 2,
    enumerate_paths = TRUE,
    max_paths = 100000
) {
  
  if (!inherits(graph, "igraph")) {
    stop("graph must be an igraph object")
  }
  
  direction <- match.arg(direction)
  
  if (!is.character(gene) || length(gene) != 1L) {
    stop("gene must be a single character string")
  }
  
  if (!gene %in% igraph::V(graph)$name) {
    stop(paste0("Gene '", gene, "' is not present in the integrated GRN"))
  }
  
  if (
    length(depth) != 1L ||
    is.na(depth) ||
    depth < 1 ||
    depth != as.integer(depth)
  ) {
    stop("depth must be a positive integer")
  }
  
  if (
    length(max_paths) != 1L ||
    is.na(max_paths) ||
    max_paths < 1
  ) {
    stop("max_paths must be a positive integer")
  }
  
  neighborhood_vertices <- igraph::ego(
    graph,
    order = depth,
    nodes = gene,
    mode = direction,
    mindist = 0
  )[[1]]
  
  neighborhood_names <- igraph::V(graph)$name[
    neighborhood_vertices
  ]
  
  subgraph <- igraph::induced_subgraph(
    graph,
    vids = neighborhood_names
  )
  
  distances_from_gene <- igraph::distances(
    graph,
    v = gene,
    to = neighborhood_names,
    mode = direction,
    weights = NA
  )
  
  node_depth <- as.integer(
    distances_from_gene[1, ]
  )
  
  node_table <- data.frame(
    gene = neighborhood_names,
    depth = node_depth,
    stringsAsFactors = FALSE
  )
  
  vertex_attributes <- igraph::vertex_attr_names(graph)
  
  for (attribute_name in setdiff(vertex_attributes, "name")) {
    
    attribute_values <- igraph::vertex_attr(
      graph,
      attribute_name,
      index = neighborhood_names
    )
    
    node_table[[attribute_name]] <- attribute_values
  }
  
  edge_table <- igraph::as_data_frame(
    subgraph,
    what = "edges"
  )
  
  colnames(edge_table)[1:2] <- c(
    "regulator",
    "target"
  )
  
  node_depth_lookup <- setNames(
    node_table$depth,
    node_table$gene
  )
  
  edge_table$regulator_depth <- unname(
    node_depth_lookup[edge_table$regulator]
  )
  
  edge_table$target_depth <- unname(
    node_depth_lookup[edge_table$target]
  )
  
  if (direction == "out") {
    
    layer_edges <- edge_table[
      edge_table$target_depth ==
        edge_table$regulator_depth + 1,
      ,
      drop = FALSE
    ]
    
  } else if (direction == "in") {
    
    layer_edges <- edge_table[
      edge_table$regulator_depth ==
        edge_table$target_depth + 1,
      ,
      drop = FALSE
    ]
    
  } else {
    
    layer_edges <- edge_table[
      abs(
        edge_table$regulator_depth -
          edge_table$target_depth
      ) == 1,
      ,
      drop = FALSE
    ]
  }
  
  layer_graph <- igraph::graph_from_data_frame(
    d = layer_edges,
    directed = TRUE,
    vertices = node_table
  )
  
  depth_summary <- as.data.frame(
    table(node_table$depth),
    stringsAsFactors = FALSE
  )
  
  colnames(depth_summary) <- c(
    "depth",
    "number_of_nodes"
  )
  
  depth_summary$depth <- as.integer(
    as.character(depth_summary$depth)
  )
  
  if ("node_type" %in% colnames(node_table)) {
    
    node_type_summary <- as.data.frame(
      table(node_table$node_type),
      stringsAsFactors = FALSE
    )
    
    colnames(node_type_summary) <- c(
      "node_type",
      "number_of_nodes"
    )
    
  } else {
    
    node_type_summary <- NULL
  }
  
  if ("edge_type" %in% colnames(edge_table)) {
    
    edge_type_summary <- as.data.frame(
      table(edge_table$edge_type),
      stringsAsFactors = FALSE
    )
    
    colnames(edge_type_summary) <- c(
      "edge_type",
      "number_of_edges"
    )
    
  } else {
    
    edge_type_summary <- NULL
  }
  
  if ("edge_type" %in% colnames(layer_edges)) {
    
    layer_edge_type_summary <- as.data.frame(
      table(layer_edges$edge_type),
      stringsAsFactors = FALSE
    )
    
    colnames(layer_edge_type_summary) <- c(
      "edge_type",
      "number_of_edges"
    )
    
  } else {
    
    layer_edge_type_summary <- NULL
  }
  
  direct_out_edges <- igraph::degree(
    graph,
    v = gene,
    mode = "out"
  )
  
  direct_in_edges <- igraph::degree(
    graph,
    v = gene,
    mode = "in"
  )
  
  summary <- data.frame(
    metric = c(
      "selected_gene",
      "direction",
      "maximum_depth",
      "direct_out_degree",
      "direct_in_degree",
      "subnetwork_nodes",
      "induced_subnetwork_edges",
      "layer_network_edges"
    ),
    value = c(
      gene,
      direction,
      depth,
      direct_out_edges,
      direct_in_edges,
      igraph::vcount(subgraph),
      igraph::ecount(subgraph),
      nrow(layer_edges)
    ),
    stringsAsFactors = FALSE
  )
  
  paths <- NULL
  paths_truncated <- FALSE
  
  if (enumerate_paths) {
    
    path_list <- igraph::all_simple_paths(
      graph,
      from = gene,
      to = neighborhood_names[
        neighborhood_names != gene
      ],
      mode = direction,
      cutoff = depth
    )
    
    if (length(path_list) > max_paths) {
      path_list <- path_list[seq_len(max_paths)]
      paths_truncated <- TRUE
    }
    
    if (length(path_list) > 0L) {
      
      paths <- data.frame(
        path_id = seq_along(path_list),
        depth = lengths(path_list) - 1L,
        path = vapply(
          path_list,
          function(x) {
            paste(
              igraph::V(graph)$name[x],
              collapse = " -> "
            )
          },
          character(1)
        ),
        stringsAsFactors = FALSE
      )
      
    } else {
      
      paths <- data.frame(
        path_id = integer(0),
        depth = integer(0),
        path = character(0),
        stringsAsFactors = FALSE
      )
    }
  }
  
  return(list(
    gene = gene,
    direction = direction,
    depth = depth,
    summary = summary,
    depth_summary = depth_summary,
    node_type_summary = node_type_summary,
    edge_type_summary = edge_type_summary,
    layer_edge_type_summary = layer_edge_type_summary,
    nodes = node_table,
    edges = edge_table,
    layer_edges = layer_edges,
    paths = paths,
    paths_truncated = paths_truncated,
    subgraph = subgraph,
    layer_graph = layer_graph
  ))
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_integrated_grn_gene_neighborhood
# ------------------------------------------------------------------------------

# Visualize a gene-centered GRN neighborhood.
#
# This function:
# - visualizes the depth-layer graph returned by explore_integrated_grn_gene()
# - optionally restricts the visualization to selected interaction classes
# - highlights the selected gene
# - colors nodes according to gene class
# - colors edges according to interaction class
# - displays node labels for local gene-centered exploration
# - preserves the complete exploration object without modifying it
# - displays only edges connecting consecutive exploration levels
#
# IMPORTANT:
# Edge filtering is applied only to the visualization.
# The original exploration object, node tables, edge tables, and path results
# remain unchanged.
#
# INPUT:
# - exploration: result returned by explore_integrated_grn_gene()
# - show_labels: whether to display gene labels
# - edge_types_to_show: optional character vector of edge types to visualize;
#   NULL displays all available interaction classes
# - vertex_size: default node size
# - selected_gene_size: size of the selected gene node
# - edge_width: graphical edge width
# - arrow_size: directed edge arrow size
#
# OUTPUT:
# - invisibly returns a list containing:
#   - layout: layout used for visualization
#   - graph: graph actually displayed after optional edge-type selection
#   - displayed_edge_types: interaction classes included in the plot
# ------------------------------------------------------------------------------

#' Plot integrated grn gene neighborhood
#'
#' Plot integrated grn gene neighborhood.
#'
#' @param exploration Result returned by `explore_integrated_grn_gene()`.
#' @param show_labels Logical indicating whether node labels should be displayed.
#' @param edge_types_to_show Optional character vector of edge classes to display.
#' @param vertex_size Default size for plotted nodes.
#' @param selected_gene_size Size used to highlight the focal gene.
#' @param edge_width Width of plotted edges.
#' @param arrow_size Size of directed edge arrows.
#'
#' @return Invisibly returns a list containing the layout, plotted graph, and displayed edge types.
#' @export
#'
#' @examples
#' # See the package vignette for end-to-end workflow examples.
plot_integrated_grn_gene_neighborhood <- function(
    exploration,
    show_labels = TRUE,
    edge_types_to_show = NULL,
    vertex_size = 5,
    selected_gene_size = 9,
    edge_width = 0.7,
    arrow_size = 0.3
) {
  
  graph <- exploration$layer_graph
  selected_gene <- exploration$gene
  
  if (!inherits(graph, "igraph")) {
    stop("exploration$layer_graph must be an igraph object")
  }
  
  available_edge_types <- unique(
    igraph::edge_attr(graph, "edge_type")
  )
  
  if (!is.null(edge_types_to_show)) {
    
    invalid_edge_types <- setdiff(
      edge_types_to_show,
      available_edge_types
    )
    
    if (length(invalid_edge_types) > 0L) {
      stop(
        "Requested edge types are not present in the selected neighborhood: ",
        paste(invalid_edge_types, collapse = ", ")
      )
    }
    
    edges_to_keep <- which(
      igraph::edge_attr(graph, "edge_type") %in%
        edge_types_to_show
    )
    
    filtered_graph <- igraph::subgraph.edges(
      graph,
      eids = edges_to_keep,
      delete.vertices = TRUE
    )
    
    if (
      selected_gene %in% igraph::V(graph)$name &&
      !selected_gene %in% igraph::V(filtered_graph)$name
    ) {
      
      filtered_graph <- igraph::add_vertices(
        filtered_graph,
        nv = 1,
        name = selected_gene
      )
      
      selected_vertex <- igraph::V(graph)[
        igraph::V(graph)$name == selected_gene
      ]
      
      vertex_attributes <- setdiff(
        igraph::vertex_attr_names(graph),
        "name"
      )
      
      for (attribute_name in vertex_attributes) {
        
        igraph::vertex_attr(
          filtered_graph,
          attribute_name,
          index = selected_gene
        ) <- igraph::vertex_attr(
          graph,
          attribute_name,
          index = selected_vertex
        )
      }
    }
    
    graph <- filtered_graph
  }
  
  if (igraph::ecount(graph) == 0L) {
    warning("No edges are available for the requested visualization")
  }
  
  node_type <- igraph::vertex_attr(
    graph,
    "node_type"
  )
  
  edge_type <- igraph::edge_attr(
    graph,
    "edge_type"
  )
  
  node_palette <- c(
    lncRNA = "#D55E00",
    TF = "#0072B2",
    PCG = "#999999"
  )
  
  edge_palette <- c(
    lncRNA_to_PCG = "#E69F00",
    lncRNA_to_TF = "#CC79A7",
    TF_to_lncRNA = "#56B4E9",
    TF_to_PCG = "#666666"
  )
  
  node_colors <- unname(
    node_palette[node_type]
  )
  
  edge_colors <- unname(
    edge_palette[edge_type]
  )
  
  vertex_sizes <- rep(
    vertex_size,
    igraph::vcount(graph)
  )
  
  vertex_sizes[
    igraph::V(graph)$name == selected_gene
  ] <- selected_gene_size
  
  if (igraph::vcount(graph) > 1L) {
    
    layout <- igraph::layout_with_fr(graph)
    
  } else {
    
    layout <- matrix(
      c(0, 0),
      nrow = 1,
      ncol = 2
    )
  }
  
  igraph::plot.igraph(
    graph,
    layout = layout,
    vertex.label = if (show_labels) {
      igraph::V(graph)$name
    } else {
      NA
    },
    vertex.label.cex = 0.6,
    vertex.size = vertex_sizes,
    vertex.color = node_colors,
    vertex.frame.color = NA,
    edge.color = edge_colors,
    edge.width = edge_width,
    edge.arrow.size = arrow_size
  )
  
  displayed_node_types <- intersect(
    names(node_palette),
    unique(node_type)
  )
  
  if (length(displayed_node_types) > 0L) {
    
    graphics::legend(
      "topleft",
      legend = displayed_node_types,
      col = node_palette[displayed_node_types],
      pch = 19,
      bty = "n",
      title = "Node type"
    )
  }
  
  displayed_edge_types <- intersect(
    names(edge_palette),
    unique(edge_type)
  )
  
  if (length(displayed_edge_types) > 0L) {
    
    graphics::legend(
      "topright",
      legend = displayed_edge_types,
      col = edge_palette[displayed_edge_types],
      lty = 1,
      lwd = 2,
      bty = "n",
      title = "Edge type"
    )
  }
  
  invisible(
    list(
      layout = layout,
      graph = graph,
      displayed_edge_types = displayed_edge_types
    )
  )
}
