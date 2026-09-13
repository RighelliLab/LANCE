# ------------------------------------------------------------------------------
# FUNCTION: load_txt_counts_to_seurat
# ------------------------------------------------------------------------------
# Load a raw UMI count matrix and the corresponding cell annotation file and
# convert them into a Seurat object.
#
# The function:
# - checks the uniqueness and correspondence of gene and cell identifiers
# - converts the count matrix into a sparse matrix
# - creates a Seurat object
# - adds the provided cell annotations
# - calculates the percentage of mitochondrial counts for each cell
# - filters cells according to the mitochondrial percentage threshold
# - reports the number of cells before and after quality-control filtering
#
# INPUT:
# - counts_file: path to the raw UMI count matrix file
# - meta_file: path to the cell annotation file
# - save_path: optional path for saving the Seurat object
#
# OPTIONAL PARAMETERS:
# - mt_pattern: pattern for identifying mitochondrial genes (default: "^MT-")
# - mt_threshold: mitochondrial percentage threshold for cell filtering (default: 5)
# - min_cells: minimum number of cells required for a feature (default: 1)
# - min_features: minimum number of features required for a cell (default: 1)
#
# OUTPUT:
# - filtered Seurat object containing raw UMI counts and cell annotations
# ------------------------------------------------------------------------------

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
  
  # I file sono già organizzati come geni × cellule.
  if (anyDuplicated(rownames(counts))) {
    stop("Gene names in the counts matrix are not unique")
  }
  
  if (anyDuplicated(colnames(counts))) {
    stop("Cell barcodes in the counts matrix are not unique")
  }
  
  if (anyDuplicated(rownames(meta))) {
    stop("Cell barcodes in the metadata are not unique")
  }
  
  # Uniforma i barcode seguendo il comportamento del codice originale.
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
  
  # Riordina le annotazioni come le colonne della matrice.
  meta <- meta[colnames(counts), , drop = FALSE]
  
  dge <- Matrix::Matrix(
    as.matrix(counts),
    sparse = TRUE
  )
  
  n_cells_before_qc <- ncol(dge)
  
  obj <- Seurat::CreateSeuratObject(
    counts = dge,
    min.cells = min_cells,
    min.features = min_features
  )
  
  n_cells_after_seurat_filter <- ncol(obj)
  
  # CreateSeuratObject può eliminare geni/cellule in base alle soglie:
  # riallineiamo quindi i metadata alle cellule effettivamente conservate.
  meta <- meta[colnames(obj), , drop = FALSE]
  
  obj <- Seurat::AddMetaData(
    object = obj,
    metadata = meta
  )
  
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    obj,
    pattern = mt_pattern
  )
  
  n_cells_removed_mt <- sum(
    obj$percent.mt >= mt_threshold
  )
  
  obj <- subset(
    obj,
    subset = percent.mt < mt_threshold
  )
  
  n_cells_after_qc <- ncol(obj)
  
  cat("\nQC summary:\n")
  cat("Cells before QC:", n_cells_before_qc, "\n")
  cat(
    "Cells removed by Seurat filtering:",
    n_cells_before_qc - n_cells_after_seurat_filter,
    "\n"
  )
  cat(
    "Cells removed by mitochondrial filtering (percent.mt >=",
    mt_threshold,
    "%):",
    n_cells_removed_mt,
    "\n"
  )
  cat("Cells after QC:", n_cells_after_qc, "\n\n")
  
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

merge_seurat_objects <- function(
    obj_list,
    save_path = NULL
) {

  # Add dataset identifier
  for (name in names(obj_list)) {
    obj_list[[name]]$GSE <- name
  }

  if (length(obj_list) < 2L) {
    stop("obj_list must contain at least two Seurat objects")
  }
  if (is.null(names(obj_list)) || any(!nzchar(names(obj_list))) ||
      anyDuplicated(names(obj_list))) {
    stop("obj_list must have unique, non-empty names")
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
plot_cell_types <- function(obj, group.by) {

  cols <- scales::hue_pal()(length(unique(obj[[group.by]][,1])))

  DimPlot(
    obj,
    group.by = group.by,
    cols = cols,
    label = TRUE,
    label.box = TRUE,
    label.color = "white"
  ) +
    theme(legend.position = "none")
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
prepare_scrmra_input <- function(obj) {

  DefaultAssay(obj) <- "RNA"
  if (!"data" %in% SeuratObject::Layers(obj[["RNA"]])) {
    obj <- Seurat::NormalizeData(obj)
  }
  obj <- JoinLayers(obj)

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
#' Title
#'
#' @param obj
#' @param logcounts
#'
#' @returns
#' @export
#'
#' @examples
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
# - prints diagnostic counts for quality control
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
  
  # In caso di ambiguità, la categoria lncRNA ha la priorità.
  # I geni ambigui vengono quindi mantenuti tra gli lncRNA.
  lnc_features <- unique(lnc_raw)
  
  # Gli stessi geni vengono esclusi dai protein-coding,
  # così i due insiemi finali rimangono disgiunti.
  prot_features <- setdiff(
    unique(prot_raw),
    ambiguous_genes
  )
  
  if (length(intersect(lnc_features, prot_features)) > 0L) {
    stop(
      "lncRNA and protein-coding feature sets are not disjoint"
    )
  }
  
  cat("\n================ DEBUG ================\n")
  cat("lnc raw:", length(lnc_raw), "\n")
  cat("prot raw:", length(prot_raw), "\n")
  cat(
    "ambiguous assigned to lncRNA:",
    length(ambiguous_genes),
    "\n"
  )
  cat("lnc final:", length(lnc_features), "\n")
  cat(
    "protein-coding final:",
    length(prot_features),
    "\n"
  )
  cat("======================================\n")
  
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
  out <- aggregate(pct.exp ~ features.plot, df, function(x) Gini(x))

  colnames(out) <- c("gene", "gini")
  return(out)
}

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

plot_gini_results <- function(GiniCoef_clean, lnc_n, prot_n) {
  
  ggplot(GiniCoef_clean, aes(x = type, y = gini, fill = type)) +
    geom_boxplot() +
    ggpubr::stat_compare_means(hjust = 0) +
    theme_bw() +
    scale_fill_manual(
      values = c(
        "lncRNA" = "#0072B2",
        "protein_coding" = "#D55E00"
      )
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      legend.position = "bottom"
    ) +
    annotate(
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
compute_gene_gini_single <- function(dotplot_data, gene) {

  df <- dotplot_data[, c("features.plot", "pct.exp")]

  out <- aggregate(pct.exp ~ features.plot, df, Gini)

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

normalize_grn_pseudobulk <- function(pseudobulk_counts) {
  
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for pseudobulk normalization")
  }
  
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

plot_grn_expression_distribution <- function(
    expression_statistics
) {
  
  mean_logcpm_plot <- ggplot(
    expression_statistics,
    aes(x = mean_logcpm)
  ) +
    geom_histogram(
      bins = 50,
      color = "black",
      fill = "#4E79A7"
    ) +
    labs(
      x = "Mean logCPM",
      y = "Number of genes",
      title = "Distribution of Mean Gene Expression"
    ) +
    theme_classic()
  
  detected_pseudobulk_plot <- ggplot(
    expression_statistics,
    aes(x = n_detected_pseudobulk)
  ) +
    geom_histogram(
      binwidth = 1,
      color = "black",
      fill = "#E15759"
    ) +
    labs(
      x = "Detected pseudobulk samples",
      y = "Number of genes",
      title = "Distribution of Gene Detection"
    ) +
    theme_classic()
  
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

plot_focused_grn_hub_lnc_network <- function(
    focused_network,
    seed = 1234
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for GRN visualization")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for GRN visualization")
  }
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for non-overlapping node labels")
  }
  
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
      title = "CRC Regulatory Hub-lncRNA Network"
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

plot_grn_lnc_target_hub_network <- function(
    target_hub_network,
    seed = 1234
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for GRN visualization")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for GRN visualization")
  }
  
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required for non-overlapping node labels")
  }
  
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
      title = "CRC lncRNA Target Hub Network"
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
      title = "Gene Biotype Composition of the Final CRC GRN"
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
      title = "Regulatory Hubs in the Final CRC GRN"
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
      title = "lncRNA Connectivity of Regulatory Hubs in the Final CRC GRN"
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
  
  if (!requireNamespace("GENIE3", quietly = TRUE)) {
    stop("Package 'GENIE3' is required for GRN inference")
  }
  
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
# quantile threshold and an absolute Spearman correlation threshold, while
# calculating a combined edge score for downstream network visualization.
#
# This function:
# - validates the required GENIE3, Spearman, and interaction-class columns
# - converts GENIE3 importance scores into class-specific percentiles
# - calculates a combined edge score from normalized GENIE3 importance and
#   absolute Spearman correlation
# - calculates the GENIE3 filtering threshold from the specified quantile
# - retains interactions with GENIE3 importance at or above the threshold
# - retains interactions with absolute Spearman correlation at or above the
#   specified threshold
# - preserves the original GENIE3 and signed Spearman values
# - preserves all original edge annotations
# - summarizes the applied filtering criteria
# - summarizes the distribution of the combined score among retained edges
#
# GENIE3 NORMALIZATION:
# - GENIE3 importance is normalized independently within each interaction class
# - rank-based percentile normalization is used
# - the minimum GENIE3 rank is mapped to 0
# - the maximum GENIE3 rank is mapped to 1
# - tied GENIE3 values receive their average rank
#
# COMBINED SCORE:
# - combined_score = genie3_percentile * abs(spearman_rho)
# - the score ranges from 0 to 1
# - high values require both high relative GENIE3 importance and strong
#   absolute Spearman correlation
# - the sign of Spearman is not included in the combined score and remains
#   available separately in spearman_rho
#
# IMPORTANT:
# The combined score is not used for edge filtering.
# Edge retention continues to depend on two independent criteria:
# - GENIE3 class-specific quantile threshold
# - absolute Spearman correlation threshold
#
# INPUT:
# - edges:
#     interaction table generated by infer_integrated_grn_relationships()
#
# - genie3_quantile:
#     class-specific GENIE3 quantile used for edge filtering
#
# - spearman_threshold:
#     minimum absolute Spearman correlation used for edge filtering
#
# OUTPUT:
# - list containing:
#
#   - edges:
#       filtered interaction table containing the original edge attributes plus:
#       - genie3_percentile
#       - combined_score
#
#   - summary:
#       filtering summary containing the applied GENIE3 and Spearman thresholds
#       and the number and fraction of retained interactions
#
#   - score_summary:
#       summary of the GENIE3 percentile and combined-score distributions among
#       the retained interactions
#
# ------------------------------------------------------------------------------

filter_integrated_grn_relationships <- function(
    edges,
    genie3_quantile = 0.95,
    spearman_threshold = 0.5
) {
  
  required_columns <- c(
    "genie3_weight",
    "spearman_rho",
    "edge_type"
  )
  
  if (!all(required_columns %in% colnames(edges))) {
    stop(
      "Required GRN columns are missing from the edge table"
    )
  }
  
  if (nrow(edges) == 0L) {
    stop(
      "The edge table is empty"
    )
  }
  
  if (
    length(genie3_quantile) != 1L ||
    is.na(genie3_quantile) ||
    genie3_quantile < 0 ||
    genie3_quantile > 1
  ) {
    stop(
      "genie3_quantile must be a single value between 0 and 1"
    )
  }
  
  if (
    length(spearman_threshold) != 1L ||
    is.na(spearman_threshold) ||
    spearman_threshold < 0 ||
    spearman_threshold > 1
  ) {
    stop(
      "spearman_threshold must be a single value between 0 and 1"
    )
  }
  
  if (anyNA(edges$genie3_weight)) {
    stop(
      "GENIE3 weights contain missing values"
    )
  }
  
  if (anyNA(edges$spearman_rho)) {
    stop(
      "Spearman correlations contain missing values"
    )
  }
  
  edge_types <- unique(
    edges$edge_type
  )
  
  if (length(edge_types) != 1L) {
    stop(
      "filter_integrated_grn_relationships() must receive exactly one edge class"
    )
  }
  
  # Normalize GENIE3 importance within the complete interaction class.
  
  genie3_rank <- rank(
    edges$genie3_weight,
    ties.method = "average",
    na.last = "keep"
  )
  
  number_of_edges <- sum(
    !is.na(edges$genie3_weight)
  )
  
  if (number_of_edges > 1L) {
    
    edges$genie3_percentile <- (
      genie3_rank - 1
    ) / (
      number_of_edges - 1
    )
    
  } else {
    
    edges$genie3_percentile <- 1
  }
  
  # Combine relative GENIE3 importance with absolute Spearman correlation.
  
  edges$combined_score <-
    edges$genie3_percentile *
    abs(edges$spearman_rho)
  
  # Calculate the class-specific GENIE3 filtering threshold.
  
  genie3_threshold <- stats::quantile(
    edges$genie3_weight,
    probs = genie3_quantile,
    na.rm = TRUE,
    names = FALSE
  )
  
  # Apply the two independent filtering criteria.
  
  keep <-
    edges$genie3_weight >= genie3_threshold &
    abs(edges$spearman_rho) >= spearman_threshold
  
  filtered_edges <- edges[
    keep & !is.na(keep),
    ,
    drop = FALSE
  ]
  
  # Filtering summary.
  
  summary <- data.frame(
    edge_type = edge_types,
    original_edges = nrow(edges),
    genie3_quantile = genie3_quantile,
    genie3_threshold = genie3_threshold,
    spearman_threshold = spearman_threshold,
    retained_edges = nrow(filtered_edges),
    retained_fraction = nrow(filtered_edges) / nrow(edges),
    stringsAsFactors = FALSE
  )
  
  # Combined-score summary for retained interactions.
  
  score_summary <- data.frame(
    edge_type = edge_types,
    retained_edges = nrow(filtered_edges),
    genie3_percentile_min = min(
      filtered_edges$genie3_percentile
    ),
    genie3_percentile_median = stats::median(
      filtered_edges$genie3_percentile
    ),
    genie3_percentile_max = max(
      filtered_edges$genie3_percentile
    ),
    combined_score_min = min(
      filtered_edges$combined_score
    ),
    combined_score_median = stats::median(
      filtered_edges$combined_score
    ),
    combined_score_max = max(
      filtered_edges$combined_score
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      edges = filtered_edges,
      summary = summary,
      score_summary = score_summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: integrate_filtered_grn_interaction_classes
# ------------------------------------------------------------------------------

# Integrate independently filtered GRN interaction classes into a single
# regulatory network.
#
# This function:
# - validates the required columns in each filtered interaction table
# - preserves original GENIE3 and Spearman edge scores
# - preserves class-specific GENIE3 percentile scores
# - preserves the combined GENIE3-Spearman edge score calculated in STEP 7B
# - combines all filtered interaction classes
# - checks for missing regulator, target, score, and annotation values
# - checks that normalized and combined edge scores are within the expected range
# - checks for duplicated regulator-target-edge_type combinations
# - constructs a node table from all regulators and targets
# - assigns each node to the lncRNA, TF, or PCG gene class
# - calculates in-degree, out-degree, and total degree in the integrated network
# - preserves all filtered edges without applying additional selection
#
# IMPORTANT:
# No filtering, normalization, or score calculation is performed here.
# `genie3_percentile` and `combined_score` must already have been calculated
# independently for each interaction class in STEP 7B.
#
# INPUT:
# - edge_tables:
#     named list containing the four filtered interaction tables
#
# - gene_sets:
#     integrated GRN gene sets containing:
#       - lncRNA
#       - TF
#       - PCG
#
# OUTPUT:
# - list containing:
#
#   - edges:
#       integrated filtered interaction table containing:
#       - regulator
#       - target
#       - genie3_weight
#       - genie3_percentile
#       - spearman_rho
#       - combined_score
#       - regulator_type
#       - target_type
#       - edge_type
#
#   - nodes:
#       integrated node annotation table containing gene class and degree
#       statistics
#
#   - summary:
#       integration and validation summary
#
# ------------------------------------------------------------------------------

integrate_filtered_grn_interaction_classes <- function(
    edge_tables,
    gene_sets
) {
  
  required_columns <- c(
    "regulator",
    "target",
    "genie3_weight",
    "genie3_percentile",
    "spearman_rho",
    "combined_score",
    "regulator_type",
    "target_type",
    "edge_type"
  )
  
  if (!is.list(edge_tables) || length(edge_tables) == 0L) {
    stop("edge_tables must be a non-empty list")
  }
  
  if (is.null(names(edge_tables)) || any(names(edge_tables) == "")) {
    stop("edge_tables must be a named list")
  }
  
  for (i in seq_along(edge_tables)) {
    
    if (!all(required_columns %in% colnames(edge_tables[[i]]))) {
      
      missing_columns <- setdiff(
        required_columns,
        colnames(edge_tables[[i]])
      )
      
      stop(
        paste0(
          "Required columns are missing from interaction table ",
          names(edge_tables)[i],
          ": ",
          paste(missing_columns, collapse = ", ")
        )
      )
    }
  }
  
  integrated_edges <- do.call(
    rbind,
    lapply(
      edge_tables,
      function(x) {
        x[, required_columns, drop = FALSE]
      }
    )
  )
  
  rownames(integrated_edges) <- NULL
  
  duplicated_edges <- duplicated(
    integrated_edges[
      ,
      c(
        "regulator",
        "target",
        "edge_type"
      )
    ]
  )
  
  invalid_genie3_percentiles <- sum(
    !is.na(integrated_edges$genie3_percentile) &
      (
        integrated_edges$genie3_percentile < 0 |
          integrated_edges$genie3_percentile > 1
      )
  )
  
  invalid_combined_scores <- sum(
    !is.na(integrated_edges$combined_score) &
      (
        integrated_edges$combined_score < 0 |
          integrated_edges$combined_score > 1
      )
  )
  
  node_names <- unique(
    c(
      integrated_edges$regulator,
      integrated_edges$target
    )
  )
  
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
  
  in_degree <- table(
    integrated_edges$target
  )
  
  out_degree <- table(
    integrated_edges$regulator
  )
  
  nodes <- data.frame(
    gene = node_names,
    node_type = node_type,
    in_degree = as.integer(
      in_degree[node_names]
    ),
    out_degree = as.integer(
      out_degree[node_names]
    ),
    stringsAsFactors = FALSE
  )
  
  nodes$in_degree[
    is.na(nodes$in_degree)
  ] <- 0L
  
  nodes$out_degree[
    is.na(nodes$out_degree)
  ] <- 0L
  
  nodes$total_degree <-
    nodes$in_degree +
    nodes$out_degree
  
  summary <- data.frame(
    metric = c(
      "number_of_edges",
      "number_of_nodes",
      "duplicated_edges",
      "missing_regulators",
      "missing_targets",
      "missing_genie3_weights",
      "missing_genie3_percentiles",
      "missing_spearman_values",
      "missing_combined_scores",
      "invalid_genie3_percentiles",
      "invalid_combined_scores",
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
      sum(is.na(integrated_edges$genie3_percentile)),
      sum(is.na(integrated_edges$spearman_rho)),
      sum(is.na(integrated_edges$combined_score)),
      invalid_genie3_percentiles,
      invalid_combined_scores,
      sum(is.na(integrated_edges$regulator_type)),
      sum(is.na(integrated_edges$target_type)),
      sum(is.na(integrated_edges$edge_type)),
      sum(is.na(nodes$node_type))
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      edges = integrated_edges,
      nodes = nodes,
      summary = summary
    )
  )
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
# - preserves the original GENIE3 importance score
# - preserves the signed Spearman correlation coefficient
# - preserves the class-specific GENIE3 percentile
# - preserves the combined edge score calculated in STEP 7B
# - does not recalculate or modify any previously computed edge score
#
# IMPORTANT:
# LncTarD is used only as an external annotation resource.
# Database support does not determine whether an edge is retained in the
# integrated GRN.
#
# Matching is performed on the directed regulator-target pair. Therefore,
# LncTarD support is assigned only when both regulator and target match in
# the same direction.
#
# INPUT:
# - integrated_edges: filtered and integrated GRN edge table containing:
#   - regulator
#   - target
#   - genie3_weight
#   - genie3_percentile
#   - spearman_rho
#   - combined_score
#   - regulator_type
#   - target_type
#   - edge_type
#
# - lnctard: LncTarD interaction table containing the regulator, target,
#   regulator type, target type, and record identifier required for
#   directed interaction matching
#
# OUTPUT:
# - list containing:
#
#   - edges:
#       complete integrated GRN edge table with the additional LncTarD
#       annotation attributes:
#       - lnctard_supported
#       - lnctard_records
#       - lnctard_rids
#
#   - lnc_edges:
#       complete subset of lncRNA-regulated GRN interactions with LncTarD
#       annotations and all previously calculated edge-score attributes
#
#   - supported_edges:
#       subset of lncRNA-regulated interactions supported by at least one
#       LncTarD record, preserving all GRN edge-score attributes
#
#   - summary:
#       summary of lncRNA-regulated interactions and LncTarD database support
#
# ------------------------------------------------------------------------------

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

# Build a directed igraph object from the complete integrated GRN.
#
# This function:
# - validates the required edge and node columns
# - constructs a directed graph from the integrated regulator-target edge table
# - preserves all edge attributes contained in the integrated GRN
# - preserves the original GENIE3 importance score
# - preserves the class-specific GENIE3 percentile
# - preserves the signed Spearman correlation coefficient
# - preserves the combined GENIE3-Spearman edge score
# - preserves interaction-class and regulator/target-type annotations
# - preserves LncTarD annotations when present
# - preserves node attributes including gene class and degree information
# - does not filter, normalize, or modify the integrated network
#
# IMPORTANT:
# igraph::graph_from_data_frame() transfers every column of the edge table
# other than regulator and target to the resulting graph as an edge attribute.
#
# The resulting graph therefore retains:
# - genie3_weight
# - genie3_percentile
# - spearman_rho
# - combined_score
# - regulator_type
# - target_type
# - edge_type
# - lnctard_supported
# - lnctard_records
# - lnctard_rids
#
# The resulting graph can be reused for:
# - global network visualization
# - edge-score-based visualization
# - inbound and outbound neighborhood exploration
# - user-defined depth analysis
# - gene-centered path analysis
# - cycle and motif exploration
#
# INPUT:
# - edges:
#     complete integrated GRN edge table containing all score and annotation
#     attributes
#
# - nodes:
#     integrated GRN node table containing gene class and degree information
#
# OUTPUT:
# - directed igraph object representing the complete integrated GRN with all
#   node and edge attributes preserved
#
# ------------------------------------------------------------------------------
build_integrated_grn_graph <- function(
    edges,
    nodes
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network analysis")
  }
  
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
# - colors nodes according to gene class
# - colors edges according to regulatory interaction class
# - represents edge weight using the combined GENIE3-Spearman score
# - scales graphical edge width according to combined_score
# - uses edge transparency to reduce overplotting
# - optionally displays legends
# - does not remove nodes or edges for visualization purposes
#
# EDGE REPRESENTATION:
# - edge color:
#     interaction class (`edge_type`)
#
# - edge width:
#     proportional to `combined_score`
#
# - combined_score:
#     genie3_percentile * abs(spearman_rho)
#
# IMPORTANT:
# Graphical edge-width scaling affects only visualization.
# The original `combined_score`, `genie3_weight`, `genie3_percentile`, and
# signed `spearman_rho` attributes stored in the graph are not modified.
#
# INPUT:
# - graph:
#     directed igraph object representing the integrated GRN
#
# - layout:
#     optional precomputed layout matrix
#
# - show_legend:
#     whether to display the edge-class legend
#
# - vertex_size:
#     graphical size of network nodes
#
# - edge_width_range:
#     numeric vector of length two defining the minimum and maximum graphical
#     edge widths used to represent combined_score
#
# - edge_alpha:
#     transparency applied to network edges
#
# - arrow_size:
#     graphical size of directed edge arrows
#
# OUTPUT:
# - invisibly returns a list containing:
#
#   - layout:
#       layout matrix used for visualization
#
#   - node_colors:
#       colors assigned to network nodes
#
#   - edge_colors:
#       colors assigned to network edges
#
#   - edge_widths:
#       graphical widths derived from combined_score
#
# ------------------------------------------------------------------------------

plot_integrated_grn <- function(
    graph,
    layout = NULL,
    show_legend = TRUE,
    vertex_size = 0.8,
    edge_width_range = c(0.03, 0.15),
    edge_alpha = 0.04,
    arrow_size = 0
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network visualization")
  }
  
  if (!inherits(graph, "igraph")) {
    stop("graph must be an igraph object")
  }
  
  required_edge_attributes <- c(
    "edge_type",
    "combined_score"
  )
  
  if (!all(required_edge_attributes %in% igraph::edge_attr_names(graph))) {
    stop(
      "Required edge attributes are missing from the integrated GRN graph"
    )
  }
  
  if (
    length(edge_width_range) != 2L ||
    any(is.na(edge_width_range)) ||
    edge_width_range[1] < 0 ||
    edge_width_range[2] <= edge_width_range[1]
  ) {
    stop(
      "edge_width_range must contain two increasing non-negative values"
    )
  }
  
  node_type <- igraph::vertex_attr(
    graph,
    "node_type"
  )
  
  edge_type <- igraph::edge_attr(
    graph,
    "edge_type"
  )
  
  combined_score <- igraph::edge_attr(
    graph,
    "combined_score"
  )
  
  if (anyNA(combined_score)) {
    stop(
      "combined_score contains missing values"
    )
  }
  
  if (is.null(layout)) {
    
    layout <- matrix(
      NA_real_,
      nrow = igraph::vcount(graph),
      ncol = 2
    )
    
    lnc_index <- which(
      node_type == "lncRNA"
    )
    
    tf_index <- which(
      node_type == "TF"
    )
    
    pcg_index <- which(
      node_type == "PCG"
    )
    
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
  
  node_colors <- unname(
    node_palette[node_type]
  )
  
  edge_colors <- grDevices::adjustcolor(
    unname(
      edge_palette[edge_type]
    ),
    alpha.f = edge_alpha
  )
  
  score_min <- min(
    combined_score,
    na.rm = TRUE
  )
  
  score_max <- max(
    combined_score,
    na.rm = TRUE
  )
  
  if (score_max > score_min) {
    
    edge_widths <-
      edge_width_range[1] +
      (
        (combined_score - score_min) /
          (score_max - score_min)
      ) *
      (
        edge_width_range[2] -
          edge_width_range[1]
      )
    
  } else {
    
    edge_widths <- rep(
      mean(edge_width_range),
      igraph::ecount(graph)
    )
  }
  
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
    edge.width = edge_widths,
    edge.arrow.size = arrow_size,
    edge.curved = 0
  )
  
  graphics::text(
    x = c(-1, 0, 1),
    y = c(1.03, 1.03, 1.03),
    labels = c(
      "lncRNA",
      "TF",
      "PCG"
    ),
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
      edge_colors = edge_colors,
      edge_widths = edge_widths
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

explore_integrated_grn_gene <- function(
    graph,
    gene,
    direction = c("out", "in", "all"),
    depth = 2,
    enumerate_paths = TRUE,
    max_paths = 100000
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network exploration")
  }
  
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

# FUNCTION: plot_integrated_grn_gene_neighborhood

# Visualize a gene-centered integrated GRN neighborhood.
#
# This function:
# - visualizes the depth-layer graph returned by explore_integrated_grn_gene()
# - optionally restricts the visualization to selected interaction classes
# - optionally filters edges according to LncTarD support
# - optionally applies a minimum absolute Spearman correlation threshold
# - highlights the selected gene
# - colors nodes according to gene class
# - colors edges according to interaction class
# - scales edge width according to the combined GENIE3-Spearman score
# - represents positive Spearman correlations with solid lines
# - represents negative Spearman correlations with dashed lines
# - displays node labels for local gene-centered exploration
# - preserves the complete exploration object without modifying it
# - displays only edges connecting consecutive exploration levels
# - calculates a force-directed Fruchterman-Reingold layout when no layout is supplied
# - allows a previously calculated layout to be reused
# - returns the exact scores and annotations of the displayed edges
#
# IMPORTANT:
# All edge filtering performed by this function is applied only to the
# visualization. The original exploration object, node tables, edge tables,
# and path results remain unchanged.
#
# Edge aesthetics:
# - edge color = interaction class
# - edge width = combined_score
# - edge line type = sign of spearman_rho
#
# INPUT:
# - exploration: result returned by explore_integrated_grn_gene()
# - show_labels: whether to display gene labels
# - edge_types_to_show: optional character vector of interaction classes;
#   NULL displays all available interaction classes
# - lnctard_filter: LncTarD visualization filter; one of "all", "supported",
#   or "unsupported"
# - min_abs_rho: optional minimum absolute Spearman correlation threshold;
#   NULL applies no additional Spearman filtering
# - vertex_size: default node size
# - selected_gene_size: size of the selected gene node
# - edge_width_range: minimum and maximum graphical edge width used to
#   represent combined_score
# - arrow_size: directed edge arrow size
# - label_cex: gene-label size
# - layout: optional precomputed two-column layout matrix;
#   NULL calculates a Fruchterman-Reingold layout
# - layout_niter: number of Fruchterman-Reingold layout iterations
# - layout_scale: multiplicative factor used to expand a newly calculated layout
# - plot_title: optional title displayed above the network
#
# OUTPUT:
# - invisibly returns a list containing:
#   - layout: layout used for visualization
#   - graph: graph actually displayed after optional visualization filtering
#   - displayed_edge_types: interaction classes included in the plot
#   - displayed_edges: number of displayed edges
#   - displayed_nodes: number of displayed nodes
#   - displayed_edge_scores: exact attributes of displayed edges
#   - edge_widths: graphical edge widths derived from combined_score

plot_integrated_grn_gene_neighborhood <- function(
    exploration,
    show_labels = TRUE,
    edge_types_to_show = NULL,
    lnctard_filter = c("all", "supported", "unsupported"),
    min_abs_rho = NULL,
    vertex_size = 5,
    selected_gene_size = 11,
    edge_width_range = c(0.25, 1.6),
    arrow_size = 0.15,
    label_cex = 0.50,
    layout = NULL,
    layout_niter = 3000,
    layout_scale = 2,
    plot_title = NULL
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for network visualization")
  }
  
  lnctard_filter <- match.arg(lnctard_filter)
  
  graph <- exploration$layer_graph
  selected_gene <- exploration$gene
  
  if (!inherits(graph, "igraph")) {
    stop("exploration$layer_graph must be an igraph object")
  }
  
  required_edge_attributes <- c(
    "edge_type",
    "genie3_weight",
    "genie3_percentile",
    "spearman_rho",
    "combined_score",
    "lnctard_supported"
  )
  
  missing_edge_attributes <- setdiff(
    required_edge_attributes,
    igraph::edge_attr_names(graph)
  )
  
  if (length(missing_edge_attributes) > 0L) {
    stop(
      "Required edge attributes are missing from the graph: ",
      paste(missing_edge_attributes, collapse = ", ")
    )
  }
  
  available_edge_types <- unique(
    igraph::edge_attr(graph, "edge_type")
  )
  
  edges_to_keep <- rep(
    TRUE,
    igraph::ecount(graph)
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
    
    edges_to_keep <- edges_to_keep &
      igraph::edge_attr(graph, "edge_type") %in% edge_types_to_show
  }
  
  if (lnctard_filter == "supported") {
    edges_to_keep <- edges_to_keep &
      igraph::edge_attr(graph, "lnctard_supported")
  }
  
  if (lnctard_filter == "unsupported") {
    edges_to_keep <- edges_to_keep &
      !igraph::edge_attr(graph, "lnctard_supported")
  }
  
  if (!is.null(min_abs_rho)) {
    
    if (
      length(min_abs_rho) != 1L ||
      !is.numeric(min_abs_rho) ||
      is.na(min_abs_rho) ||
      min_abs_rho < 0 ||
      min_abs_rho > 1
    ) {
      stop("min_abs_rho must be NULL or a numeric value between 0 and 1")
    }
    
    edges_to_keep <- edges_to_keep &
      abs(igraph::edge_attr(graph, "spearman_rho")) >= min_abs_rho
  }
  
  filtered_graph <- igraph::subgraph_from_edges(
    graph,
    eids = which(edges_to_keep),
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
  
  spearman_rho <- igraph::edge_attr(
    graph,
    "spearman_rho"
  )
  
  combined_score <- igraph::edge_attr(
    graph,
    "combined_score"
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
  
  node_colors <- grDevices::adjustcolor(
    unname(node_palette[node_type]),
    alpha.f = 0.70
  )
  
  edge_colors <- grDevices::adjustcolor(
    unname(edge_palette[edge_type]),
    alpha.f = 0.25
  )
  
  vertex_sizes <- rep(
    vertex_size,
    igraph::vcount(graph)
  )
  
  vertex_sizes[
    igraph::V(graph)$name == selected_gene
  ] <- selected_gene_size
  
  if (igraph::ecount(graph) > 0L) {
    
    score_min <- min(
      combined_score,
      na.rm = TRUE
    )
    
    score_max <- max(
      combined_score,
      na.rm = TRUE
    )
    
    if (
      is.finite(score_min) &&
      is.finite(score_max) &&
      score_max > score_min
    ) {
      
      edge_widths <- edge_width_range[1] +
        (
          (combined_score - score_min) /
            (score_max - score_min)
        ) *
        (
          edge_width_range[2] -
            edge_width_range[1]
        )
      
    } else {
      
      edge_widths <- rep(
        mean(edge_width_range),
        igraph::ecount(graph)
      )
    }
    
  } else {
    
    edge_widths <- numeric(0)
  }
  
  edge_lty <- ifelse(
    spearman_rho >= 0,
    1,
    2
  )
  
  if (is.null(layout)) {
    
    if (igraph::vcount(graph) > 1L) {
      
      layout <- igraph::layout_with_fr(
        graph,
        niter = layout_niter
      )
      
      layout <- layout * layout_scale
      
    } else {
      
      layout <- matrix(
        c(0, 0),
        nrow = 1,
        ncol = 2
      )
    }
    
  } else {
    
    if (
      !is.matrix(layout) ||
      nrow(layout) != igraph::vcount(graph) ||
      ncol(layout) != 2L
    ) {
      stop(
        "Provided layout must be a two-column matrix with one row per displayed node"
      )
    }
  }
  
  old_par <- graphics::par(
    mar = c(1, 1, 4, 1),
    xpd = NA
  )
  
  on.exit(
    graphics::par(old_par),
    add = TRUE
  )
  
  vertex_labels <- rep(
    NA_character_,
    igraph::vcount(graph)
  )
  
  if (show_labels) {
    vertex_labels[
      igraph::V(graph)$name == selected_gene
    ] <- selected_gene
  }
  
  igraph::plot.igraph(
    graph,
    layout = layout,
    vertex.label = vertex_labels,
    vertex.label.cex = label_cex,
    vertex.label.color = "#00008B",
    vertex.size = vertex_sizes,
    vertex.color = node_colors,
    vertex.frame.color = NA,
    edge.color = edge_colors,
    edge.width = edge_widths,
    edge.lty = edge_lty,
    edge.arrow.size = arrow_size,
    rescale = TRUE
  )
  
  if (!is.null(plot_title)) {
    
    graphics::title(
      main = plot_title,
      cex.main = 0.9,
      font.main = 2,
      line = 1
    )
  }
  
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
      pt.cex = 1.1,
      bty = "n",
      title = "Node type",
      cex = 0.85
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
      title = "Edge type",
      cex = 0.85
    )
  }
  
  if (igraph::ecount(graph) > 0L) {
    
    graphics::legend(
      "right",
      legend = c(
        "Positive rho",
        "Negative rho"
      ),
      col = "#666666",
      lty = c(1, 2),
      lwd = 2,
      bty = "n",
      title = "Spearman sign",
      cex = 0.85
    )
  }
  
  if (igraph::ecount(graph) > 0L) {
    
    displayed_edge_scores <- igraph::as_data_frame(
      graph,
      what = "edges"
    )
    
    colnames(displayed_edge_scores)[
      colnames(displayed_edge_scores) == "from"
    ] <- "regulator"
    
    colnames(displayed_edge_scores)[
      colnames(displayed_edge_scores) == "to"
    ] <- "target"
    
    displayed_edge_scores$abs_spearman_rho <- abs(
      displayed_edge_scores$spearman_rho
    )
    
  } else {
    
    displayed_edge_scores <- data.frame()
  }
  
  invisible(
    list(
      layout = layout,
      graph = graph,
      displayed_edge_types = displayed_edge_types,
      displayed_edges = igraph::ecount(graph),
      displayed_nodes = igraph::vcount(graph),
      displayed_edge_scores = displayed_edge_scores,
      edge_widths = edge_widths
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_integrated_grn_topology
# ------------------------------------------------------------------------------

# Perform an extended topological characterization of the complete integrated GRN.
#
# This function:
# - validates the directed integrated GRN
# - calculates in-degree, out-degree, and total degree for every node
# - classifies nodes according to their directed network role
# - calculates network density
# - calculates weakly and strongly connected components
# - summarizes degree distributions
# - preserves the original graph without modifying nodes or edges
#
# NODE ROLES:
# - source_only:
#     node with outgoing interactions but no incoming interactions
# - target_only:
#     node with incoming interactions but no outgoing interactions
# - regulator_and_target:
#     node with both incoming and outgoing interactions
# - isolated:
#     node with no incoming or outgoing interactions
#
# INPUT:
# - graph:
#     directed igraph object representing the complete integrated GRN
#
# OUTPUT:
# - list containing:
#
#   - network_summary:
#       global network statistics
#
#   - node_connectivity_summary:
#       summary statistics for in-degree, out-degree, and total degree
#
#   - node_table:
#       node-level table containing:
#       - gene
#       - node_type
#       - in_degree
#       - out_degree
#       - total_degree
#       - network_role
#
# ------------------------------------------------------------------------------

summarize_integrated_grn_topology <- function(
    graph
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for integrated GRN topology analysis"
    )
  }
  
  if (!inherits(graph, "igraph")) {
    stop(
      "graph must be an igraph object"
    )
  }
  
  if (!igraph::is_directed(graph)) {
    stop(
      "The integrated GRN must be a directed graph"
    )
  }
  
  if (!"node_type" %in% igraph::vertex_attr_names(graph)) {
    stop(
      "The integrated GRN graph does not contain the 'node_type' vertex attribute"
    )
  }
  
  node_names <- igraph::V(graph)$name
  
  if (is.null(node_names)) {
    stop(
      "The integrated GRN graph does not contain node names"
    )
  }
  
  in_degree <- igraph::degree(
    graph,
    mode = "in",
    loops = FALSE
  )
  
  out_degree <- igraph::degree(
    graph,
    mode = "out",
    loops = FALSE
  )
  
  total_degree <- igraph::degree(
    graph,
    mode = "all",
    loops = FALSE
  )
  
  network_role <- ifelse(
    out_degree > 0 & in_degree == 0,
    "source_only",
    ifelse(
      in_degree > 0 & out_degree == 0,
      "target_only",
      ifelse(
        in_degree > 0 & out_degree > 0,
        "regulator_and_target",
        "isolated"
      )
    )
  )
  
  node_table <- data.frame(
    gene = node_names,
    node_type = igraph::vertex_attr(
      graph,
      "node_type"
    ),
    in_degree = as.integer(in_degree),
    out_degree = as.integer(out_degree),
    total_degree = as.integer(total_degree),
    network_role = network_role,
    stringsAsFactors = FALSE
  )
  
  weak_components <- igraph::components(
    graph,
    mode = "weak"
  )
  
  strong_components <- igraph::components(
    graph,
    mode = "strong"
  )
  
  network_summary <- data.frame(
    metric = c(
      "number_of_nodes",
      "number_of_edges",
      "network_density",
      "weakly_connected_components",
      "strongly_connected_components",
      "source_only_nodes",
      "target_only_nodes",
      "regulator_and_target_nodes",
      "isolated_nodes"
    ),
    value = c(
      igraph::vcount(graph),
      igraph::ecount(graph),
      igraph::edge_density(
        graph,
        loops = FALSE
      ),
      weak_components$no,
      strong_components$no,
      sum(node_table$network_role == "source_only"),
      sum(node_table$network_role == "target_only"),
      sum(node_table$network_role == "regulator_and_target"),
      sum(node_table$network_role == "isolated")
    ),
    stringsAsFactors = FALSE
  )
  
  node_connectivity_summary <- data.frame(
    metric = c(
      "mean_in_degree",
      "median_in_degree",
      "maximum_in_degree",
      "mean_out_degree",
      "median_out_degree",
      "maximum_out_degree",
      "mean_total_degree",
      "median_total_degree",
      "maximum_total_degree"
    ),
    value = c(
      mean(node_table$in_degree),
      stats::median(node_table$in_degree),
      max(node_table$in_degree),
      mean(node_table$out_degree),
      stats::median(node_table$out_degree),
      max(node_table$out_degree),
      mean(node_table$total_degree),
      stats::median(node_table$total_degree),
      max(node_table$total_degree)
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      network_summary = network_summary,
      node_connectivity_summary = node_connectivity_summary,
      node_table = node_table
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: evaluate_integrated_grn_hub_thresholds
# ------------------------------------------------------------------------------

# Evaluate candidate degree thresholds for hub identification in the integrated GRN.
#
# This function:
# - evaluates multiple high degree-distribution quantiles
# - evaluates regulatory hubs using out-degree
# - evaluates target hubs using in-degree
# - evaluates global hubs using total degree
# - calculates the degree cutoff associated with each candidate quantile
# - reports the number and proportion of nodes retained at each threshold
# - does not classify or remove nodes
#
# Regulatory thresholds are calculated only among nodes with out-degree > 0.
# Target thresholds are calculated only among nodes with in-degree > 0.
# Global thresholds are calculated using all nodes in the integrated GRN.
#
# INPUT:
# - node_table:
#     node-level table generated by summarize_integrated_grn_topology()
#
# - quantiles:
#     numeric vector of candidate degree quantiles to evaluate
#
# OUTPUT:
# - data.frame containing:
#     - hub_type
#     - quantile
#     - degree_cutoff
#     - eligible_nodes
#     - hub_nodes
#     - hub_fraction
#
# ------------------------------------------------------------------------------

evaluate_integrated_grn_hub_thresholds <- function(
    node_table,
    quantiles = c(0.90, 0.95, 0.99, 0.995)
) {
  
  required_columns <- c(
    "gene",
    "node_type",
    "in_degree",
    "out_degree",
    "total_degree"
  )
  
  if (!all(required_columns %in% colnames(node_table))) {
    stop(
      "node_table must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (
    anyNA(quantiles) ||
    any(quantiles <= 0) ||
    any(quantiles >= 1)
  ) {
    stop(
      "quantiles must contain values strictly between 0 and 1"
    )
  }
  
  hub_metrics <- list(
    regulatory = node_table$out_degree[
      node_table$out_degree > 0
    ],
    target = node_table$in_degree[
      node_table$in_degree > 0
    ],
    global = node_table$total_degree
  )
  
  results <- list()
  result_index <- 1L
  
  for (hub_type in names(hub_metrics)) {
    
    values <- hub_metrics[[hub_type]]
    
    if (length(values) == 0L) {
      stop(
        "No eligible nodes found for hub type: ",
        hub_type
      )
    }
    
    for (hub_quantile in quantiles) {
      
      degree_cutoff <- as.numeric(
        stats::quantile(
          values,
          probs = hub_quantile,
          na.rm = TRUE,
          names = FALSE
        )
      )
      
      hub_nodes <- sum(
        values >= degree_cutoff
      )
      
      results[[result_index]] <- data.frame(
        hub_type = hub_type,
        quantile = hub_quantile,
        degree_cutoff = degree_cutoff,
        eligible_nodes = length(values),
        hub_nodes = hub_nodes,
        hub_fraction = hub_nodes / length(values),
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
# FUNCTION: identify_integrated_grn_hubs
# ------------------------------------------------------------------------------

# Identify highly connected hubs in the complete integrated GRN.
#
# This function:
# - identifies regulatory hubs using high out-degree
# - identifies target hubs using high in-degree
# - identifies global hubs using high total degree
# - calculates degree cutoffs independently for the three hub definitions
# - preserves node gene-class annotations and connectivity statistics
# - allows the same node to belong to more than one hub category
#
# HUB DEFINITIONS:
# - regulatory hub:
#     out-degree at or above the selected quantile among nodes with out-degree > 0
#
# - target hub:
#     in-degree at or above the selected quantile among nodes with in-degree > 0
#
# - global hub:
#     total degree at or above the selected quantile among all GRN nodes
#
# INPUT:
# - node_table:
#     node-level table generated by summarize_integrated_grn_topology()
#
# - hub_quantile:
#     degree-distribution quantile used to identify highly connected nodes
#
# OUTPUT:
# - list containing:
#
#   - regulatory_hubs:
#       regulatory hubs ordered by decreasing out-degree
#
#   - target_hubs:
#       target hubs ordered by decreasing in-degree
#
#   - global_hubs:
#       global hubs ordered by decreasing total degree
#
#   - thresholds:
#       degree cutoffs and number of hubs identified for each hub category
#
# ------------------------------------------------------------------------------

identify_integrated_grn_hubs <- function(
    node_table,
    hub_quantile = 0.99
) {
  
  required_columns <- c(
    "gene",
    "node_type",
    "in_degree",
    "out_degree",
    "total_degree",
    "network_role"
  )
  
  if (!all(required_columns %in% colnames(node_table))) {
    stop(
      "node_table must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (
    length(hub_quantile) != 1L ||
    is.na(hub_quantile) ||
    hub_quantile <= 0 ||
    hub_quantile >= 1
  ) {
    stop(
      "hub_quantile must be a single value strictly between 0 and 1"
    )
  }
  
  regulatory_nodes <- node_table[
    node_table$out_degree > 0,
    ,
    drop = FALSE
  ]
  
  target_nodes <- node_table[
    node_table$in_degree > 0,
    ,
    drop = FALSE
  ]
  
  if (nrow(regulatory_nodes) == 0L) {
    stop(
      "No regulatory nodes with out-degree > 0 were found"
    )
  }
  
  if (nrow(target_nodes) == 0L) {
    stop(
      "No target nodes with in-degree > 0 were found"
    )
  }
  
  regulatory_cutoff <- as.numeric(
    stats::quantile(
      regulatory_nodes$out_degree,
      probs = hub_quantile,
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  target_cutoff <- as.numeric(
    stats::quantile(
      target_nodes$in_degree,
      probs = hub_quantile,
      na.rm = TRUE,
      names = FALSE
    )
  )
  
  global_cutoff <- as.numeric(
    stats::quantile(
      node_table$total_degree,
      probs = hub_quantile,
      na.rm = TRUE,
      names = FALSE
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
  
  global_hubs <- node_table[
    node_table$total_degree >= global_cutoff,
    ,
    drop = FALSE
  ]
  
  regulatory_hubs <- regulatory_hubs[
    order(
      -regulatory_hubs$out_degree,
      -regulatory_hubs$total_degree
    ),
    ,
    drop = FALSE
  ]
  
  target_hubs <- target_hubs[
    order(
      -target_hubs$in_degree,
      -target_hubs$total_degree
    ),
    ,
    drop = FALSE
  ]
  
  global_hubs <- global_hubs[
    order(
      -global_hubs$total_degree
    ),
    ,
    drop = FALSE
  ]
  
  thresholds <- data.frame(
    hub_type = c(
      "regulatory",
      "target",
      "global"
    ),
    quantile = rep(
      hub_quantile,
      3
    ),
    degree_metric = c(
      "out_degree",
      "in_degree",
      "total_degree"
    ),
    degree_cutoff = c(
      regulatory_cutoff,
      target_cutoff,
      global_cutoff
    ),
    hub_nodes = c(
      nrow(regulatory_hubs),
      nrow(target_hubs),
      nrow(global_hubs)
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
# FUNCTION: summarize_integrated_grn_hub_composition
# ------------------------------------------------------------------------------

# Summarize the gene-class composition of integrated GRN hubs.
#
# This function:
# - evaluates regulatory, target, and global hub tables
# - counts hubs belonging to each integrated GRN gene class
# - calculates gene-class proportions within each hub category
# - preserves zero-count gene classes for direct comparison
#
# EXPECTED GENE CLASSES:
# - TF
# - lncRNA
# - PCG
#
# INPUT:
# - hubs:
#     list generated by identify_integrated_grn_hubs()
#
# OUTPUT:
# - data.frame containing:
#     - hub_type
#     - node_type
#     - n_hubs
#     - total_hubs
#     - proportion
#
# ------------------------------------------------------------------------------

summarize_integrated_grn_hub_composition <- function(
    hubs
) {
  
  required_hub_tables <- c(
    "regulatory_hubs",
    "target_hubs",
    "global_hubs"
  )
  
  if (!all(required_hub_tables %in% names(hubs))) {
    stop(
      "hubs must contain: ",
      paste(required_hub_tables, collapse = ", ")
    )
  }
  
  expected_node_types <- c(
    "TF",
    "lncRNA",
    "PCG"
  )
  
  hub_tables <- list(
    regulatory = hubs$regulatory_hubs,
    target = hubs$target_hubs,
    global = hubs$global_hubs
  )
  
  results <- lapply(
    names(hub_tables),
    function(hub_type) {
      
      current_hubs <- hub_tables[[hub_type]]
      
      if (!"node_type" %in% colnames(current_hubs)) {
        stop(
          "Hub table does not contain 'node_type': ",
          hub_type
        )
      }
      
      hub_counts <- table(
        factor(
          current_hubs$node_type,
          levels = expected_node_types
        )
      )
      
      total_hubs <- nrow(current_hubs)
      
      data.frame(
        hub_type = hub_type,
        node_type = expected_node_types,
        n_hubs = as.integer(hub_counts),
        total_hubs = total_hubs,
        proportion = if (total_hubs > 0L) {
          as.integer(hub_counts) / total_hubs
        } else {
          0
        },
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
# FUNCTION: plot_integrated_grn_degree_distributions
# ------------------------------------------------------------------------------

# Visualize in-degree, out-degree, and total-degree distributions of the integrated GRN.
#
# This function:
# - converts node-level degree measurements to long format
# - plots the three degree distributions separately using facets
# - displays node counts on a logarithmic scale
# - adds the selected hub cutoff for each degree metric
#
# INPUT:
# - node_table:
#     node-level table generated by summarize_integrated_grn_topology()
#
# - hub_thresholds:
#     threshold table generated by identify_integrated_grn_hubs()
#
# OUTPUT:
# - ggplot object containing the degree-distribution visualization
#
# ------------------------------------------------------------------------------

plot_integrated_grn_degree_distributions <- function(
    node_table,
    hub_thresholds
) {
  
  required_columns <- c(
    "in_degree",
    "out_degree",
    "total_degree"
  )
  
  if (!all(required_columns %in% colnames(node_table))) {
    stop(
      "node_table must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  required_threshold_columns <- c(
    "degree_metric",
    "degree_cutoff"
  )
  
  if (!all(required_threshold_columns %in% colnames(hub_thresholds))) {
    stop(
      "hub_thresholds must contain: ",
      paste(required_threshold_columns, collapse = ", ")
    )
  }
  
  degree_data <- data.frame(
    degree_metric = rep(
      c("in_degree", "out_degree", "total_degree"),
      each = nrow(node_table)
    ),
    degree = c(
      node_table$in_degree,
      node_table$out_degree,
      node_table$total_degree
    ),
    stringsAsFactors = FALSE
  )
  
  degree_data$degree_metric <- factor(
    degree_data$degree_metric,
    levels = c(
      "in_degree",
      "out_degree",
      "total_degree"
    ),
    labels = c(
      "In-degree",
      "Out-degree",
      "Total degree"
    )
  )
  
  threshold_data <- hub_thresholds[
    hub_thresholds$degree_metric %in%
      c("in_degree", "out_degree", "total_degree"),
    c("degree_metric", "degree_cutoff"),
    drop = FALSE
  ]
  
  threshold_data$degree_metric <- factor(
    threshold_data$degree_metric,
    levels = c(
      "in_degree",
      "out_degree",
      "total_degree"
    ),
    labels = c(
      "In-degree",
      "Out-degree",
      "Total degree"
    )
  )
  
  p <- ggplot2::ggplot(
    degree_data,
    ggplot2::aes(x = degree)
  ) +
    ggplot2::geom_histogram(
      bins = 60
    ) +
    ggplot2::geom_vline(
      data = threshold_data,
      ggplot2::aes(xintercept = degree_cutoff),
      linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::facet_wrap(
      ~ degree_metric,
      scales = "free_x",
      ncol = 1
    ) +
    ggplot2::scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10)
    ) +
    ggplot2::labs(
      title = "Integrated GRN degree distributions",
      x = "Degree",
      y = "Number of nodes (pseudo-log scale)"
      ) +
    ggplot2::theme_classic()
  
  return(p)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_integrated_grn_hub_composition
# ------------------------------------------------------------------------------

# Visualize gene-class composition across integrated GRN hub categories.
#
# This function:
# - uses the hub-composition summary generated in STEP 13D
# - compares TF, lncRNA, and PCG representation
# - displays absolute hub counts for each hub category
#
# INPUT:
# - hub_composition:
#     data.frame generated by summarize_integrated_grn_hub_composition()
#
# OUTPUT:
# - ggplot object showing hub composition by gene class
#
# ------------------------------------------------------------------------------

plot_integrated_grn_hub_composition <- function(
    hub_composition
) {
  
  required_columns <- c(
    "hub_type",
    "node_type",
    "n_hubs"
  )
  
  if (!all(required_columns %in% colnames(hub_composition))) {
    stop(
      "hub_composition must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  hub_composition$hub_type <- factor(
    hub_composition$hub_type,
    levels = c(
      "regulatory",
      "target",
      "global"
    ),
    labels = c(
      "Regulatory",
      "Target",
      "Global"
    )
  )
  
  p <- ggplot2::ggplot(
    hub_composition,
    ggplot2::aes(
      x = hub_type,
      y = n_hubs,
      fill = node_type
    )
  ) +
    ggplot2::geom_col(
      position = "stack"
    ) +
    ggplot2::labs(
      title = "Integrated GRN hub composition",
      x = "Hub category",
      y = "Number of hubs",
      fill = "Gene type"
    ) +
    ggplot2::theme_classic()
  
  return(p)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_integrated_grn_regulatory_hubs
# ------------------------------------------------------------------------------

# Visualize all regulatory hubs identified in the integrated GRN.
#
# This function:
# - uses the complete regulatory-hub table generated in STEP 13C
# - orders all regulatory hubs according to decreasing out-degree
# - displays the out-degree of every identified regulatory hub
# - distinguishes regulatory hubs according to gene class
# - reduces bar width and gene-label size to improve readability
# - does not select, filter, or exclude any identified hub
#
# INPUT:
# - regulatory_hubs:
#     complete regulatory-hub table generated by
#     identify_integrated_grn_hubs()
#
# OUTPUT:
# - ggplot object showing all identified regulatory hubs
#
# ------------------------------------------------------------------------------

plot_integrated_grn_regulatory_hubs <- function(
    regulatory_hubs
) {
  
  required_columns <- c(
    "gene",
    "node_type",
    "out_degree"
  )
  
  if (!all(required_columns %in% colnames(regulatory_hubs))) {
    stop(
      "regulatory_hubs must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (nrow(regulatory_hubs) == 0L) {
    stop(
      "No regulatory hubs are available for visualization"
    )
  }
  
  plot_data <- regulatory_hubs[
    order(
      -regulatory_hubs$out_degree
    ),
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
      y = out_degree,
      fill = node_type
    )
  ) +
    ggplot2::geom_col(
      width = 0.65
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Integrated GRN regulatory hubs",
      x = "Gene",
      y = "Out-degree",
      fill = "Gene type"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(
        size = 7
      ),
      axis.title.y = ggplot2::element_text(
        margin = ggplot2::margin(r = 8)
      )
    )
  
  return(plot)
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_integrated_grn_for_leiden
# ------------------------------------------------------------------------------

# Prepare the complete integrated GRN for Leiden community detection.
#
# This function:
# - validates the complete directed integrated GRN
# - verifies the presence of the selected edge-weight attribute
# - verifies that all selected edge weights are finite and strictly positive
# - converts the directed GRN to an undirected graph using mode = "each"
# - preserves one undirected edge for every directed edge in the original GRN
# - preserves all vertices and vertex attributes
# - assigns the selected GRN edge score to the standard Leiden "weight" attribute
# - checks whether node and edge counts are preserved after conversion
# - reports self-loops and multiple edges without removing them
# - does not modify the original directed GRN
#
# IMPORTANT:
# The undirected graph is generated exclusively for community detection.
# Direction-dependent GRN analyses must continue to use the original directed graph.
#
# INPUT:
# - graph:
#     directed igraph object representing the complete integrated GRN
#
# - weight_attribute:
#     name of the positive edge attribute to use as Leiden interaction strength
#
# OUTPUT:
# - list containing:
#
#   - graph:
#       complete undirected graph prepared for Leiden community detection
#
#   - summary:
#       graph-conversion and integrity summary
#
#   - weight_summary:
#       summary statistics of the Leiden edge weights
#
# ------------------------------------------------------------------------------

prepare_integrated_grn_for_leiden <- function(
    graph,
    weight_attribute = "combined_score"
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for Leiden community detection"
    )
  }
  
  if (!inherits(graph, "igraph")) {
    stop(
      "graph must be an igraph object"
    )
  }
  
  if (!igraph::is_directed(graph)) {
    stop(
      "The input integrated GRN must be directed"
    )
  }
  
  if (!weight_attribute %in% igraph::edge_attr_names(graph)) {
    stop(
      "The requested Leiden weight attribute is missing from the graph: ",
      weight_attribute
    )
  }
  
  edge_weights <- igraph::edge_attr(
    graph,
    weight_attribute
  )
  
  if (!is.numeric(edge_weights)) {
    stop(
      "The Leiden edge-weight attribute must be numeric"
    )
  }
  
  if (anyNA(edge_weights)) {
    stop(
      "The Leiden edge-weight attribute contains missing values"
    )
  }
  
  if (any(!is.finite(edge_weights))) {
    stop(
      "The Leiden edge-weight attribute contains non-finite values"
    )
  }
  
  if (any(edge_weights <= 0)) {
    stop(
      "Leiden edge weights must be strictly positive"
    )
  }
  
  original_nodes <- igraph::vcount(graph)
  original_edges <- igraph::ecount(graph)
  
  leiden_graph <- igraph::as_undirected(
    graph,
    mode = "each"
  )
  
  leiden_graph <- igraph::set_edge_attr(
    leiden_graph,
    name = "weight",
    value = igraph::edge_attr(
      leiden_graph,
      weight_attribute
    )
  )
  
  leiden_nodes <- igraph::vcount(
    leiden_graph
  )
  
  leiden_edges <- igraph::ecount(
    leiden_graph
  )
  
  if (leiden_nodes != original_nodes) {
    stop(
      "Node count changed during Leiden graph preparation"
    )
  }
  
  if (leiden_edges != original_edges) {
    stop(
      "Edge count changed during Leiden graph preparation"
    )
  }
  
  leiden_weights <- igraph::edge_attr(
    leiden_graph,
    "weight"
  )
  
  summary <- data.frame(
    metric = c(
      "original_directed_nodes",
      "original_directed_edges",
      "leiden_undirected_nodes",
      "leiden_undirected_edges",
      "nodes_preserved",
      "edges_preserved",
      "leiden_graph_directed",
      "self_loops",
      "multiple_edges",
      "missing_weights",
      "nonpositive_weights"
    ),
    value = c(
      original_nodes,
      original_edges,
      leiden_nodes,
      leiden_edges,
      leiden_nodes == original_nodes,
      leiden_edges == original_edges,
      igraph::is_directed(leiden_graph),
      sum(igraph::which_loop(leiden_graph)),
      sum(igraph::which_multiple(leiden_graph)),
      sum(is.na(leiden_weights)),
      sum(leiden_weights <= 0)
    ),
    stringsAsFactors = FALSE
  )
  
  weight_summary <- data.frame(
    metric = c(
      "minimum_weight",
      "first_quartile_weight",
      "median_weight",
      "mean_weight",
      "third_quartile_weight",
      "maximum_weight"
    ),
    value = c(
      min(leiden_weights),
      as.numeric(stats::quantile(leiden_weights, 0.25)),
      stats::median(leiden_weights),
      mean(leiden_weights),
      as.numeric(stats::quantile(leiden_weights, 0.75)),
      max(leiden_weights)
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      graph = leiden_graph,
      summary = summary,
      weight_summary = weight_summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: run_integrated_grn_leiden
# ------------------------------------------------------------------------------

# Run Leiden community detection on the prepared integrated GRN.
#
# This function:
# - validates the undirected Leiden graph
# - verifies the presence and validity of the edge-weight attribute
# - runs weighted Leiden community detection
# - fixes the random seed for reproducibility
# - extracts the complete community membership vector
# - calculates community sizes
# - calculates modularity of the resulting partition
# - preserves all nodes and does not filter or remove communities
#
# INPUT:
# - graph:
#     undirected igraph object generated by prepare_integrated_grn_for_leiden()
#
# - objective_function:
#     Leiden objective function; supported values are "modularity" and "CPM"
#
# - resolution:
#     Leiden resolution parameter
#
# - beta:
#     parameter controlling randomness during the Leiden refinement step
#
# - n_iterations:
#     number of Leiden iterations
#
# - seed:
#     random seed used for reproducibility
#
# OUTPUT:
# - list containing:
#
#   - communities:
#       igraph communities object returned by cluster_leiden()
#
#   - membership:
#       complete node-level community membership vector
#
#   - community_sizes:
#       data.frame containing the number of nodes in every detected community
#
#   - summary:
#       summary of Leiden parameters and resulting partition
#
# ------------------------------------------------------------------------------

run_integrated_grn_leiden <- function(
    graph,
    objective_function = "modularity",
    resolution = 1,
    beta = 0.01,
    n_iterations = 2,
    seed = 1
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for Leiden community detection"
    )
  }
  
  if (!inherits(graph, "igraph")) {
    stop(
      "graph must be an igraph object"
    )
  }
  
  if (igraph::is_directed(graph)) {
    stop(
      "The Leiden graph must be undirected"
    )
  }
  
  if (!objective_function %in% c("modularity", "CPM")) {
    stop(
      "objective_function must be either 'modularity' or 'CPM'"
    )
  }
  
  if (
    length(resolution) != 1L ||
    is.na(resolution) ||
    !is.finite(resolution) ||
    resolution <= 0
  ) {
    stop(
      "resolution must be a single positive finite value"
    )
  }
  
  if (
    length(beta) != 1L ||
    is.na(beta) ||
    !is.finite(beta) ||
    beta < 0
  ) {
    stop(
      "beta must be a single non-negative finite value"
    )
  }
  
  if (
    length(n_iterations) != 1L ||
    is.na(n_iterations) ||
    n_iterations < 1
  ) {
    stop(
      "n_iterations must be a positive integer"
    )
  }
  
  if (
    length(seed) != 1L ||
    is.na(seed)
  ) {
    stop(
      "seed must be a single non-missing value"
    )
  }
  
  if (!"weight" %in% igraph::edge_attr_names(graph)) {
    stop(
      "The Leiden graph does not contain a 'weight' edge attribute"
    )
  }
  
  edge_weights <- igraph::edge_attr(
    graph,
    "weight"
  )
  
  if (
    anyNA(edge_weights) ||
    any(!is.finite(edge_weights)) ||
    any(edge_weights <= 0)
  ) {
    stop(
      "Leiden edge weights must be finite, non-missing, and strictly positive"
    )
  }
  
  set.seed(seed)
  
  leiden_communities <- igraph::cluster_leiden(
    graph = graph,
    objective_function = objective_function,
    weights = edge_weights,
    resolution = resolution,
    beta = beta,
    n_iterations = n_iterations
  )
  
  membership <- igraph::membership(
    leiden_communities
  )
  
  if (length(membership) != igraph::vcount(graph)) {
    stop(
      "Leiden membership does not contain all graph nodes"
    )
  }
  
  community_size_vector <- igraph::sizes(
    leiden_communities
  )
  
  community_sizes <- data.frame(
    community = as.integer(
      names(community_size_vector)
    ),
    n_nodes = as.integer(
      community_size_vector
    ),
    stringsAsFactors = FALSE
  )
  
  community_sizes <- community_sizes[
    order(
      -community_sizes$n_nodes,
      community_sizes$community
    ),
    ,
    drop = FALSE
  ]
  
  rownames(community_sizes) <- NULL
  
  partition_modularity <- igraph::modularity(
    graph,
    membership = membership,
    weights = edge_weights,
    resolution = resolution
  )
  
  summary <- data.frame(
    metric = c(
      "objective_function",
      "resolution",
      "beta",
      "n_iterations",
      "seed",
      "number_of_nodes",
      "number_of_edges",
      "number_of_communities",
      "modularity",
      "minimum_community_size",
      "median_community_size",
      "mean_community_size",
      "maximum_community_size"
    ),
    value = c(
      objective_function,
      resolution,
      beta,
      n_iterations,
      seed,
      igraph::vcount(graph),
      igraph::ecount(graph),
      length(community_size_vector),
      partition_modularity,
      min(community_size_vector),
      stats::median(community_size_vector),
      mean(community_size_vector),
      max(community_size_vector)
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      communities = leiden_communities,
      membership = membership,
      community_sizes = community_sizes,
      summary = summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: assign_integrated_grn_leiden_membership
# ------------------------------------------------------------------------------

# Assign Leiden community membership to the complete integrated GRN node table.
#
# This function:
# - validates the node-level topology table
# - validates the Leiden membership vector
# - matches community assignments to genes using vertex names
# - preserves all node-level topology and gene-class annotations
# - checks that every node receives exactly one community assignment
# - checks for missing and duplicated gene identifiers
# - does not filter, merge, renumber, or remove communities
#
# INPUT:
# - node_table:
#     complete node-level table generated by
#     summarize_integrated_grn_topology()
#
# - membership:
#     Leiden membership vector generated by run_integrated_grn_leiden()
#
# OUTPUT:
# - list containing:
#
#   - nodes:
#       complete GRN node table with an additional leiden_community column
#
#   - summary:
#       integrity summary of the community-assignment procedure
#
# ------------------------------------------------------------------------------

assign_integrated_grn_leiden_membership <- function(
    node_table,
    membership
) {
  
  required_columns <- c(
    "gene",
    "node_type",
    "in_degree",
    "out_degree",
    "total_degree",
    "network_role"
  )
  
  if (!all(required_columns %in% colnames(node_table))) {
    stop(
      "node_table must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (anyDuplicated(node_table$gene) > 0L) {
    stop(
      "node_table contains duplicated gene identifiers"
    )
  }
  
  if (is.null(names(membership))) {
    stop(
      "Leiden membership vector must contain node names"
    )
  }
  
  if (anyDuplicated(names(membership)) > 0L) {
    stop(
      "Leiden membership contains duplicated node names"
    )
  }
  
  if (anyNA(membership)) {
    stop(
      "Leiden membership contains missing community assignments"
    )
  }
  
  membership_index <- match(
    node_table$gene,
    names(membership)
  )
  
  if (anyNA(membership_index)) {
    stop(
      "Some GRN nodes are missing from the Leiden membership vector"
    )
  }
  
  leiden_community <- membership[
    membership_index
  ]
  
  if (length(leiden_community) != nrow(node_table)) {
    stop(
      "Leiden community assignment does not contain all GRN nodes"
    )
  }
  
  nodes <- node_table
  
  nodes$leiden_community <- as.integer(
    leiden_community
  )
  
  if (anyNA(nodes$leiden_community)) {
    stop(
      "Missing Leiden community assignments were generated"
    )
  }
  
  summary <- data.frame(
    metric = c(
      "number_of_nodes",
      "membership_nodes",
      "assigned_nodes",
      "nodes_preserved",
      "duplicated_genes",
      "missing_community_assignments",
      "number_of_communities"
    ),
    value = c(
      nrow(node_table),
      length(membership),
      nrow(nodes),
      nrow(nodes) == nrow(node_table),
      sum(duplicated(nodes$gene)),
      sum(is.na(nodes$leiden_community)),
      length(unique(nodes$leiden_community))
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      nodes = nodes,
      summary = summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_integrated_grn_leiden_communities
# ------------------------------------------------------------------------------

# Summarize the size and basic structure of Leiden communities in the integrated GRN.
#
# This function:
# - uses the complete node-level Leiden assignment
# - uses the original directed integrated GRN edge table
# - reports node counts and node-type composition for every community
# - counts intra-community and inter-community edges
# - calculates the proportion of edges remaining within each community
# - summarizes node-level in-degree, out-degree, and total-degree statistics
# - does not filter, merge, renumber, or remove communities
#
# INPUT:
# - node_assignment:
#     complete node table containing leiden_community generated in STEP 14C
#
# - edges:
#     complete directed integrated GRN edge table
#
# OUTPUT:
# - list containing:
#
#   - community_summary:
#       one row per Leiden community with structural statistics
#
#   - edge_summary:
#       detailed counts of intra-community and inter-community edges
#
# ------------------------------------------------------------------------------

summarize_integrated_grn_leiden_communities <- function(
    node_assignment,
    edges
) {
  
  required_node_columns <- c(
    "gene",
    "node_type",
    "in_degree",
    "out_degree",
    "total_degree",
    "leiden_community"
  )
  
  if (!all(required_node_columns %in% colnames(node_assignment))) {
    stop(
      "node_assignment must contain: ",
      paste(required_node_columns, collapse = ", ")
    )
  }
  
  required_edge_columns <- c(
    "regulator",
    "target"
  )
  
  if (!all(required_edge_columns %in% colnames(edges))) {
    stop(
      "edges must contain: ",
      paste(required_edge_columns, collapse = ", ")
    )
  }
  
  if (anyDuplicated(node_assignment$gene) > 0L) {
    stop(
      "node_assignment contains duplicated gene identifiers"
    )
  }
  
  if (anyNA(node_assignment$leiden_community)) {
    stop(
      "node_assignment contains missing Leiden community assignments"
    )
  }
  
  regulator_index <- match(
    edges$regulator,
    node_assignment$gene
  )
  
  target_index <- match(
    edges$target,
    node_assignment$gene
  )
  
  if (anyNA(regulator_index) || anyNA(target_index)) {
    stop(
      "Some GRN edge genes are missing from the Leiden node assignment"
    )
  }
  
  regulator_community <- node_assignment$leiden_community[
    regulator_index
  ]
  
  target_community <- node_assignment$leiden_community[
    target_index
  ]
  
  intra_community <- regulator_community == target_community
  
  community_ids <- sort(
    unique(
      node_assignment$leiden_community
    )
  )
  
  total_nodes <- nrow(
    node_assignment
  )
  
  community_summary_list <- lapply(
    community_ids,
    function(current_community) {
      
      current_nodes <- node_assignment[
        node_assignment$leiden_community == current_community,
        ,
        drop = FALSE
      ]
      
      n_nodes <- nrow(
        current_nodes
      )
      
      n_tf <- sum(
        current_nodes$node_type == "TF",
        na.rm = TRUE
      )
      
      n_lncrna <- sum(
        current_nodes$node_type == "lncRNA",
        na.rm = TRUE
      )
      
      n_pcg <- sum(
        current_nodes$node_type == "PCG",
        na.rm = TRUE
      )
      
      intra_edges <- sum(
        intra_community &
          regulator_community == current_community
      )
      
      outgoing_inter_edges <- sum(
        !intra_community &
          regulator_community == current_community
      )
      
      incoming_inter_edges <- sum(
        !intra_community &
          target_community == current_community
      )
      
      total_incident_edges <- (
        intra_edges +
          outgoing_inter_edges +
          incoming_inter_edges
      )
      
      intra_edge_fraction <- if (
        total_incident_edges > 0L
      ) {
        intra_edges / total_incident_edges
      } else {
        NA_real_
      }
      
      data.frame(
        community = current_community,
        n_nodes = n_nodes,
        node_fraction = n_nodes / total_nodes,
        n_TF = n_tf,
        n_lncRNA = n_lncrna,
        n_PCG = n_pcg,
        intra_edges = intra_edges,
        outgoing_inter_edges = outgoing_inter_edges,
        incoming_inter_edges = incoming_inter_edges,
        intra_edge_fraction = intra_edge_fraction,
        mean_in_degree = mean(current_nodes$in_degree),
        mean_out_degree = mean(current_nodes$out_degree),
        mean_total_degree = mean(current_nodes$total_degree),
        stringsAsFactors = FALSE
      )
    }
  )
  
  community_summary <- do.call(
    rbind,
    community_summary_list
  )
  
  rownames(
    community_summary
  ) <- NULL
  
  edge_summary <- data.frame(
    regulator_community = regulator_community,
    target_community = target_community,
    intra_community = intra_community,
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      community_summary = community_summary,
      edge_summary = edge_summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_leiden_gene_type_composition
# ------------------------------------------------------------------------------

# Summarize TF, lncRNA, and PCG composition across Leiden communities.
#
# This function:
# - uses the complete Leiden node-assignment table
# - counts TF, lncRNA, and PCG genes within every community
# - calculates within-community gene-class proportions
# - retains zero-count gene classes explicitly
# - verifies that all nodes are represented exactly once in the summary
# - does not filter, merge, renumber, or remove communities
#
# GENE CLASSES:
# - TF
# - lncRNA
# - PCG
#
# INPUT:
# - node_assignment:
#     complete node-level table generated in STEP 14C containing:
#     - gene
#     - node_type
#     - leiden_community
#
# OUTPUT:
# - list containing:
#
#   - composition:
#       long-format table containing:
#       - community
#       - node_type
#       - n_nodes
#       - community_total
#       - proportion
#
#   - summary:
#       integrity summary of the community-composition analysis
#
# ------------------------------------------------------------------------------

summarize_leiden_gene_type_composition <- function(
    node_assignment
) {
  
  required_columns <- c(
    "gene",
    "node_type",
    "leiden_community"
  )
  
  if (!all(required_columns %in% colnames(node_assignment))) {
    stop(
      "node_assignment must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (anyDuplicated(node_assignment$gene) > 0L) {
    stop(
      "node_assignment contains duplicated gene identifiers"
    )
  }
  
  if (anyNA(node_assignment$leiden_community)) {
    stop(
      "node_assignment contains missing Leiden community assignments"
    )
  }
  
  expected_node_types <- c(
    "TF",
    "lncRNA",
    "PCG"
  )
  
  unexpected_node_types <- setdiff(
    unique(node_assignment$node_type),
    expected_node_types
  )
  
  unexpected_node_types <- unexpected_node_types[
    !is.na(unexpected_node_types)
  ]
  
  if (length(unexpected_node_types) > 0L) {
    stop(
      "Unexpected node types found: ",
      paste(unexpected_node_types, collapse = ", ")
    )
  }
  
  if (anyNA(node_assignment$node_type)) {
    stop(
      "node_assignment contains missing node_type values"
    )
  }
  
  community_ids <- sort(
    unique(
      node_assignment$leiden_community
    )
  )
  
  composition_list <- lapply(
    community_ids,
    function(current_community) {
      
      current_nodes <- node_assignment[
        node_assignment$leiden_community == current_community,
        ,
        drop = FALSE
      ]
      
      community_total <- nrow(
        current_nodes
      )
      
      node_type_counts <- table(
        factor(
          current_nodes$node_type,
          levels = expected_node_types
        )
      )
      
      data.frame(
        community = current_community,
        node_type = expected_node_types,
        n_nodes = as.integer(node_type_counts),
        community_total = community_total,
        proportion = as.integer(node_type_counts) / community_total,
        stringsAsFactors = FALSE
      )
    }
  )
  
  composition <- do.call(
    rbind,
    composition_list
  )
  
  rownames(composition) <- NULL
  
  summarized_nodes <- sum(
    composition$n_nodes
  )
  
  summary <- data.frame(
    metric = c(
      "number_of_nodes",
      "summarized_nodes",
      "nodes_preserved",
      "number_of_communities",
      "number_of_gene_classes",
      "missing_node_types",
      "missing_community_assignments"
    ),
    value = c(
      nrow(node_assignment),
      summarized_nodes,
      summarized_nodes == nrow(node_assignment),
      length(community_ids),
      length(expected_node_types),
      sum(is.na(node_assignment$node_type)),
      sum(is.na(node_assignment$leiden_community))
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      composition = composition,
      summary = summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: characterize_leiden_intra_community_structure
# ------------------------------------------------------------------------------

# Characterize the internal directed-network structure of each Leiden community.
#
# This function:
# - uses the complete directed integrated GRN
# - uses the complete Leiden node assignment generated in STEP 14C
# - constructs the node-induced directed subgraph of every Leiden community
# - retains all nodes assigned to each community
# - retains all directed edges connecting nodes within the same community
# - calculates community-specific directed-network topology
# - calculates in-degree, out-degree, and total degree using only intra-community edges
# - does not filter nodes, edges, or communities
#
# INPUT:
# - graph:
#     complete directed igraph object representing the integrated GRN
#
# - node_assignment:
#     complete node-level table generated in STEP 14C containing:
#     - gene
#     - leiden_community
#
# OUTPUT:
# - list containing:
#
#   - community_summary:
#       one row per Leiden community containing internal topological statistics
#
#   - node_table:
#       complete node-level table containing community-specific in-degree,
#       out-degree, total degree, and intra-community network role
#
# ------------------------------------------------------------------------------

characterize_leiden_intra_community_structure <- function(
    graph,
    node_assignment
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for intra-community network analysis"
    )
  }
  
  if (!inherits(graph, "igraph")) {
    stop(
      "graph must be an igraph object"
    )
  }
  
  if (!igraph::is_directed(graph)) {
    stop(
      "The integrated GRN must be directed"
    )
  }
  
  required_columns <- c(
    "gene",
    "leiden_community"
  )
  
  if (!all(required_columns %in% colnames(node_assignment))) {
    stop(
      "node_assignment must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (anyDuplicated(node_assignment$gene) > 0L) {
    stop(
      "node_assignment contains duplicated gene identifiers"
    )
  }
  
  if (anyNA(node_assignment$leiden_community)) {
    stop(
      "node_assignment contains missing Leiden community assignments"
    )
  }
  
  graph_genes <- igraph::V(graph)$name
  
  if (is.null(graph_genes)) {
    stop(
      "Graph vertices must contain gene names"
    )
  }
  
  if (!setequal(graph_genes, node_assignment$gene)) {
    stop(
      "Graph vertices and Leiden node assignments do not contain the same genes"
    )
  }
  
  community_ids <- sort(
    unique(node_assignment$leiden_community)
  )
  
  community_summary_list <- vector(
    "list",
    length(community_ids)
  )
  
  node_table_list <- vector(
    "list",
    length(community_ids)
  )
  
  for (i in seq_along(community_ids)) {
    
    current_community <- community_ids[i]
    
    current_genes <- node_assignment$gene[
      node_assignment$leiden_community == current_community
    ]
    
    community_graph <- igraph::induced_subgraph(
      graph,
      vids = current_genes
    )
    
    in_degree <- igraph::degree(
      community_graph,
      mode = "in",
      loops = FALSE
    )
    
    out_degree <- igraph::degree(
      community_graph,
      mode = "out",
      loops = FALSE
    )
    
    total_degree <- igraph::degree(
      community_graph,
      mode = "all",
      loops = FALSE
    )
    
    network_role <- ifelse(
      in_degree == 0 & out_degree == 0,
      "isolated",
      ifelse(
        in_degree == 0 & out_degree > 0,
        "source_only",
        ifelse(
          in_degree > 0 & out_degree == 0,
          "target_only",
          "regulator_and_target"
        )
      )
    )
    
    weak_components <- igraph::components(
      community_graph,
      mode = "weak"
    )$no
    
    strong_components <- igraph::components(
      community_graph,
      mode = "strong"
    )$no
    
    community_summary_list[[i]] <- data.frame(
      community = current_community,
      n_nodes = igraph::vcount(community_graph),
      intra_edges = igraph::ecount(community_graph),
      density = igraph::edge_density(
        community_graph,
        loops = FALSE
      ),
      weakly_connected_components = weak_components,
      strongly_connected_components = strong_components,
      source_only_nodes = sum(network_role == "source_only"),
      target_only_nodes = sum(network_role == "target_only"),
      regulator_and_target_nodes = sum(network_role == "regulator_and_target"),
      isolated_nodes = sum(network_role == "isolated"),
      mean_in_degree = mean(in_degree),
      median_in_degree = stats::median(in_degree),
      maximum_in_degree = max(in_degree),
      mean_out_degree = mean(out_degree),
      median_out_degree = stats::median(out_degree),
      maximum_out_degree = max(out_degree),
      mean_total_degree = mean(total_degree),
      median_total_degree = stats::median(total_degree),
      maximum_total_degree = max(total_degree),
      stringsAsFactors = FALSE
    )
    
    node_table_list[[i]] <- data.frame(
      gene = names(in_degree),
      leiden_community = current_community,
      intra_in_degree = as.numeric(in_degree),
      intra_out_degree = as.numeric(out_degree),
      intra_total_degree = as.numeric(total_degree),
      intra_network_role = network_role,
      stringsAsFactors = FALSE
    )
  }
  
  community_summary <- do.call(
    rbind,
    community_summary_list
  )
  
  node_table <- do.call(
    rbind,
    node_table_list
  )
  
  rownames(community_summary) <- NULL
  rownames(node_table) <- NULL
  
  if (sum(community_summary$n_nodes) != igraph::vcount(graph)) {
    stop(
      "Community subgraphs do not preserve all integrated GRN nodes"
    )
  }
  
  if (nrow(node_table) != igraph::vcount(graph)) {
    stop(
      "Intra-community node table does not contain all integrated GRN nodes"
    )
  }
  
  return(
    list(
      community_summary = community_summary,
      node_table = node_table
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: identify_leiden_community_hubs
# ------------------------------------------------------------------------------

# Identify highly connected nodes independently within each Leiden community.
#
# This function:
# - uses intra-community node connectivity calculated in STEP 15B
# - evaluates each Leiden community independently
# - identifies regulatory hubs using intra-community out-degree
# - identifies target hubs using intra-community in-degree
# - identifies global community hubs using intra-community total degree
# - calculates community-specific degree cutoffs
# - preserves all node annotations
# - allows the same node to belong to more than one hub category
# - does not remove or merge communities
#
# HUB DEFINITIONS:
# - regulatory hub:
#     intra-community out-degree at or above the selected quantile among nodes
#     with intra_out_degree > 0 in the corresponding community
#
# - target hub:
#     intra-community in-degree at or above the selected quantile among nodes
#     with intra_in_degree > 0 in the corresponding community
#
# - global community hub:
#     intra-community total degree at or above the selected quantile among all
#     nodes in the corresponding community
#
# INPUT:
# - intra_node_table:
#     node-level table generated by
#     characterize_leiden_intra_community_structure()
#
# - node_assignment:
#     complete Leiden node-assignment table generated in STEP 14C,
#     used to restore gene-class annotations
#
# - hub_quantile:
#     degree-distribution quantile used to identify community-specific hubs
#
# OUTPUT:
# - list containing:
#
#   - regulatory_hubs:
#       all community-specific regulatory hubs
#
#   - target_hubs:
#       all community-specific target hubs
#
#   - global_hubs:
#       all community-specific global hubs
#
#   - thresholds:
#       community-specific cutoffs and hub counts for all three hub definitions
#
# ------------------------------------------------------------------------------

identify_leiden_community_hubs <- function(
    intra_node_table,
    node_assignment,
    hub_quantile = 0.99
) {
  
  required_intra_columns <- c(
    "gene",
    "leiden_community",
    "intra_in_degree",
    "intra_out_degree",
    "intra_total_degree",
    "intra_network_role"
  )
  
  if (!all(required_intra_columns %in% colnames(intra_node_table))) {
    stop(
      "intra_node_table must contain: ",
      paste(required_intra_columns, collapse = ", ")
    )
  }
  
  required_assignment_columns <- c(
    "gene",
    "node_type"
  )
  
  if (!all(required_assignment_columns %in% colnames(node_assignment))) {
    stop(
      "node_assignment must contain: ",
      paste(required_assignment_columns, collapse = ", ")
    )
  }
  
  if (
    length(hub_quantile) != 1L ||
    is.na(hub_quantile) ||
    !is.finite(hub_quantile) ||
    hub_quantile <= 0 ||
    hub_quantile >= 1
  ) {
    stop(
      "hub_quantile must be a single value strictly between 0 and 1"
    )
  }
  
  if (anyDuplicated(intra_node_table$gene) > 0L) {
    stop(
      "intra_node_table contains duplicated gene identifiers"
    )
  }
  
  if (anyDuplicated(node_assignment$gene) > 0L) {
    stop(
      "node_assignment contains duplicated gene identifiers"
    )
  }
  
  node_type_index <- match(
    intra_node_table$gene,
    node_assignment$gene
  )
  
  if (anyNA(node_type_index)) {
    stop(
      "Some intra-community genes are missing from node_assignment"
    )
  }
  
  nodes <- intra_node_table
  
  nodes$node_type <- node_assignment$node_type[
    node_type_index
  ]
  
  community_ids <- sort(
    unique(nodes$leiden_community)
  )
  
  regulatory_hub_list <- list()
  target_hub_list <- list()
  global_hub_list <- list()
  threshold_list <- list()
  
  result_index <- 1L
  
  for (current_community in community_ids) {
    
    current_nodes <- nodes[
      nodes$leiden_community == current_community,
      ,
      drop = FALSE
    ]
    
    regulatory_nodes <- current_nodes[
      current_nodes$intra_out_degree > 0,
      ,
      drop = FALSE
    ]
    
    target_nodes <- current_nodes[
      current_nodes$intra_in_degree > 0,
      ,
      drop = FALSE
    ]
    
    if (nrow(regulatory_nodes) == 0L) {
      stop(
        "No regulatory nodes found in Leiden community ",
        current_community
      )
    }
    
    if (nrow(target_nodes) == 0L) {
      stop(
        "No target nodes found in Leiden community ",
        current_community
      )
    }
    
    regulatory_cutoff <- as.numeric(
      stats::quantile(
        regulatory_nodes$intra_out_degree,
        probs = hub_quantile,
        na.rm = TRUE,
        names = FALSE
      )
    )
    
    target_cutoff <- as.numeric(
      stats::quantile(
        target_nodes$intra_in_degree,
        probs = hub_quantile,
        na.rm = TRUE,
        names = FALSE
      )
    )
    
    global_cutoff <- as.numeric(
      stats::quantile(
        current_nodes$intra_total_degree,
        probs = hub_quantile,
        na.rm = TRUE,
        names = FALSE
      )
    )
    
    regulatory_hubs <- regulatory_nodes[
      regulatory_nodes$intra_out_degree >= regulatory_cutoff,
      ,
      drop = FALSE
    ]
    
    target_hubs <- target_nodes[
      target_nodes$intra_in_degree >= target_cutoff,
      ,
      drop = FALSE
    ]
    
    global_hubs <- current_nodes[
      current_nodes$intra_total_degree >= global_cutoff,
      ,
      drop = FALSE
    ]
    
    regulatory_hubs <- regulatory_hubs[
      order(
        -regulatory_hubs$intra_out_degree,
        -regulatory_hubs$intra_total_degree
      ),
      ,
      drop = FALSE
    ]
    
    target_hubs <- target_hubs[
      order(
        -target_hubs$intra_in_degree,
        -target_hubs$intra_total_degree
      ),
      ,
      drop = FALSE
    ]
    
    global_hubs <- global_hubs[
      order(
        -global_hubs$intra_total_degree
      ),
      ,
      drop = FALSE
    ]
    
    regulatory_hub_list[[result_index]] <- regulatory_hubs
    target_hub_list[[result_index]] <- target_hubs
    global_hub_list[[result_index]] <- global_hubs
    
    threshold_list[[result_index]] <- data.frame(
      community = current_community,
      hub_type = c(
        "regulatory",
        "target",
        "global"
      ),
      quantile = rep(
        hub_quantile,
        3
      ),
      degree_metric = c(
        "intra_out_degree",
        "intra_in_degree",
        "intra_total_degree"
      ),
      degree_cutoff = c(
        regulatory_cutoff,
        target_cutoff,
        global_cutoff
      ),
      eligible_nodes = c(
        nrow(regulatory_nodes),
        nrow(target_nodes),
        nrow(current_nodes)
      ),
      hub_nodes = c(
        nrow(regulatory_hubs),
        nrow(target_hubs),
        nrow(global_hubs)
      ),
      stringsAsFactors = FALSE
    )
    
    result_index <- result_index + 1L
  }
  
  regulatory_hubs <- do.call(
    rbind,
    regulatory_hub_list
  )
  
  target_hubs <- do.call(
    rbind,
    target_hub_list
  )
  
  global_hubs <- do.call(
    rbind,
    global_hub_list
  )
  
  thresholds <- do.call(
    rbind,
    threshold_list
  )
  
  rownames(regulatory_hubs) <- NULL
  rownames(target_hubs) <- NULL
  rownames(global_hubs) <- NULL
  rownames(thresholds) <- NULL
  
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
# FUNCTION: characterize_leiden_tf_lncrna_profiles
# ------------------------------------------------------------------------------

# Characterize TF and lncRNA genes across Leiden communities and annotate
# community-specific hub status.
#
# This function:
# - uses the complete intra-community node table generated in STEP 15B
# - restores node-type annotations from the complete Leiden node assignment
# - retains all TF and lncRNA genes
# - annotates regulatory, target, and global community-hub membership
# - summarizes hub representation by community and gene class
# - does not apply top-N selection or additional filtering
#
# INPUT:
# - intra_node_table:
#     complete intra-community node table generated in STEP 15B
#
# - node_assignment:
#     complete Leiden node-assignment table generated in STEP 14C
#
# - community_hubs:
#     hub-analysis object generated in STEP 15C
#
# OUTPUT:
# - list containing:
#
#   - genes:
#       complete TF and lncRNA node-level table with hub-status annotations
#
#   - summary:
#       community-level TF/lncRNA hub summary
#
# ------------------------------------------------------------------------------

characterize_leiden_tf_lncrna_profiles <- function(
    intra_node_table,
    node_assignment,
    community_hubs
) {
  
  required_intra_columns <- c(
    "gene",
    "leiden_community",
    "intra_in_degree",
    "intra_out_degree",
    "intra_total_degree",
    "intra_network_role"
  )
  
  if (!all(required_intra_columns %in% colnames(intra_node_table))) {
    stop(
      "intra_node_table must contain: ",
      paste(required_intra_columns, collapse = ", ")
    )
  }
  
  required_assignment_columns <- c(
    "gene",
    "node_type"
  )
  
  if (!all(required_assignment_columns %in% colnames(node_assignment))) {
    stop(
      "node_assignment must contain: ",
      paste(required_assignment_columns, collapse = ", ")
    )
  }
  
  required_hub_objects <- c(
    "regulatory_hubs",
    "target_hubs",
    "global_hubs"
  )
  
  if (!all(required_hub_objects %in% names(community_hubs))) {
    stop(
      "community_hubs must contain: ",
      paste(required_hub_objects, collapse = ", ")
    )
  }
  
  if (anyDuplicated(intra_node_table$gene) > 0L) {
    stop(
      "intra_node_table contains duplicated gene identifiers"
    )
  }
  
  if (anyDuplicated(node_assignment$gene) > 0L) {
    stop(
      "node_assignment contains duplicated gene identifiers"
    )
  }
  
  node_type_index <- match(
    intra_node_table$gene,
    node_assignment$gene
  )
  
  if (anyNA(node_type_index)) {
    stop(
      "Some intra-community genes are missing from node_assignment"
    )
  }
  
  genes <- intra_node_table
  
  genes$node_type <- node_assignment$node_type[
    node_type_index
  ]
  
  genes <- genes[
    genes$node_type %in% c(
      "TF",
      "lncRNA"
    ),
    ,
    drop = FALSE
  ]
  
  genes$is_regulatory_hub <- genes$gene %in%
    community_hubs$regulatory_hubs$gene
  
  genes$is_target_hub <- genes$gene %in%
    community_hubs$target_hubs$gene
  
  genes$is_global_hub <- genes$gene %in%
    community_hubs$global_hubs$gene
  
  genes <- genes[
    order(
      genes$leiden_community,
      genes$node_type,
      -genes$intra_total_degree,
      genes$gene
    ),
    ,
    drop = FALSE
  ]
  
  rownames(genes) <- NULL
  
  community_ids <- sort(
    unique(
      genes$leiden_community
    )
  )
  
  gene_types <- c(
    "TF",
    "lncRNA"
  )
  
  summary_list <- list()
  result_index <- 1L
  
  for (current_community in community_ids) {
    
    for (current_type in gene_types) {
      
      current_genes <- genes[
        genes$leiden_community == current_community &
          genes$node_type == current_type,
        ,
        drop = FALSE
      ]
      
      summary_list[[result_index]] <- data.frame(
        community = current_community,
        node_type = current_type,
        n_genes = nrow(current_genes),
        regulatory_hubs = sum(current_genes$is_regulatory_hub),
        target_hubs = sum(current_genes$is_target_hub),
        global_hubs = sum(current_genes$is_global_hub),
        stringsAsFactors = FALSE
      )
      
      result_index <- result_index + 1L
    }
  }
  
  summary <- do.call(
    rbind,
    summary_list
  )
  
  rownames(summary) <- NULL
  
  return(
    list(
      genes = genes,
      summary = summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_leiden_functional_gene_sets
# ------------------------------------------------------------------------------

# Prepare community-specific gene sets for downstream functional enrichment.
#
# This function:
# - uses the complete Leiden node assignment
# - separates genes by Leiden community
# - preserves TF, lncRNA, and PCG gene classes
# - creates protein-coding enrichment gene sets using TF + PCG genes
# - creates separate lncRNA gene sets for downstream functional association
# - does not apply degree, hub, score, top-N, or community-size filtering
#
# INPUT:
# - node_assignment:
#     complete Leiden node-assignment table generated in STEP 14C
#
# OUTPUT:
# - list containing:
#
#   - enrichment_gene_sets:
#       named list containing all TF + PCG genes for each Leiden community
#
#   - tf_gene_sets:
#       named list containing all TF genes for each Leiden community
#
#   - lncrna_gene_sets:
#       named list containing all lncRNA genes for each Leiden community
#
#   - pcg_gene_sets:
#       named list containing all PCG genes for each Leiden community
#
#   - summary:
#       number of genes retained in each gene set for every community
#
# ------------------------------------------------------------------------------

prepare_leiden_functional_gene_sets <- function(
    node_assignment
) {
  
  required_columns <- c(
    "gene",
    "node_type",
    "leiden_community"
  )
  
  if (!all(required_columns %in% colnames(node_assignment))) {
    stop(
      "node_assignment must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (anyDuplicated(node_assignment$gene) > 0L) {
    stop(
      "node_assignment contains duplicated gene identifiers"
    )
  }
  
  if (anyNA(node_assignment$leiden_community)) {
    stop(
      "node_assignment contains missing Leiden community assignments"
    )
  }
  
  if (anyNA(node_assignment$node_type)) {
    stop(
      "node_assignment contains missing node-type annotations"
    )
  }
  
  expected_node_types <- c(
    "TF",
    "lncRNA",
    "PCG"
  )
  
  unexpected_node_types <- setdiff(
    unique(node_assignment$node_type),
    expected_node_types
  )
  
  if (length(unexpected_node_types) > 0L) {
    stop(
      "Unexpected node types found: ",
      paste(unexpected_node_types, collapse = ", ")
    )
  }
  
  community_ids <- sort(
    unique(node_assignment$leiden_community)
  )
  
  enrichment_gene_sets <- list()
  tf_gene_sets <- list()
  lncrna_gene_sets <- list()
  pcg_gene_sets <- list()
  summary_list <- list()
  
  for (i in seq_along(community_ids)) {
    
    current_community <- community_ids[i]
    
    current_nodes <- node_assignment[
      node_assignment$leiden_community == current_community,
      ,
      drop = FALSE
    ]
    
    tf_genes <- current_nodes$gene[
      current_nodes$node_type == "TF"
    ]
    
    lncrna_genes <- current_nodes$gene[
      current_nodes$node_type == "lncRNA"
    ]
    
    pcg_genes <- current_nodes$gene[
      current_nodes$node_type == "PCG"
    ]
    
    enrichment_genes <- current_nodes$gene[
      current_nodes$node_type %in% c(
        "TF",
        "PCG"
      )
    ]
    
    community_name <- paste0(
      "community_",
      current_community
    )
    
    enrichment_gene_sets[[community_name]] <- enrichment_genes
    tf_gene_sets[[community_name]] <- tf_genes
    lncrna_gene_sets[[community_name]] <- lncrna_genes
    pcg_gene_sets[[community_name]] <- pcg_genes
    
    summary_list[[i]] <- data.frame(
      community = current_community,
      total_nodes = nrow(current_nodes),
      n_TF = length(tf_genes),
      n_lncRNA = length(lncrna_genes),
      n_PCG = length(pcg_genes),
      n_enrichment_genes = length(enrichment_genes),
      stringsAsFactors = FALSE
    )
  }
  
  summary <- do.call(
    rbind,
    summary_list
  )
  
  rownames(summary) <- NULL
  
  if (sum(summary$total_nodes) != nrow(node_assignment)) {
    stop(
      "Community gene sets do not preserve all Leiden-assigned nodes"
    )
  }
  
  if (
    sum(summary$n_TF) +
    sum(summary$n_lncRNA) +
    sum(summary$n_PCG) != nrow(node_assignment)
  ) {
    stop(
      "Gene-class-specific gene sets do not preserve all nodes"
    )
  }
  
  return(
    list(
      enrichment_gene_sets = enrichment_gene_sets,
      tf_gene_sets = tf_gene_sets,
      lncrna_gene_sets = lncrna_gene_sets,
      pcg_gene_sets = pcg_gene_sets,
      summary = summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: perform_leiden_go_enrichment
# ------------------------------------------------------------------------------

# Perform Gene Ontology enrichment independently for each Leiden community.
#
# This function:
# - uses the complete TF + PCG gene sets prepared for functional analysis
# - defines the enrichment universe from all TF + PCG genes represented in the GRN
# - maps gene symbols to Entrez Gene identifiers
# - performs GO enrichment independently for every community
# - supports Biological Process (BP), Molecular Function (MF),
#   and Cellular Component (CC)
# - allows one or multiple GO ontologies to be analysed in the same call
# - applies multiple-testing correction independently within each enrichment
# - retains the complete enrichment output
# - annotates statistical significance without filtering the stored results
# - records gene-mapping statistics for every community
#
# INPUT:
# - enrichment_gene_sets:
#     named list containing TF + PCG gene symbols for each Leiden community
#
# - universe_genes:
#     complete set of TF + PCG gene symbols represented in the GRN
#
# - ontologies:
#     one or more GO ontologies among:
#       "BP" = Biological Process
#       "MF" = Molecular Function
#       "CC" = Cellular Component
#
# - p_adjust_method:
#     multiple-testing correction method passed to enrichGO()
#
# - significance_cutoff:
#     adjusted p-value threshold used to annotate significant GO terms
#
# OUTPUT:
# - list containing:
#   - results_by_ontology
#   - combined_results
#   - mapping_summary
#   - universe_mapping
#
# ------------------------------------------------------------------------------

perform_leiden_go_enrichment <- function(
    enrichment_gene_sets,
    universe_genes,
    ontologies = c("BP"),
    p_adjust_method = "BH",
    significance_cutoff = 0.05
) {
  
  # --------------------------------------------------------------------------
  # Validate input
  # --------------------------------------------------------------------------
  
  if (!is.list(enrichment_gene_sets)) {
    stop("enrichment_gene_sets must be a list")
  }
  
  if (length(enrichment_gene_sets) == 0L) {
    stop("enrichment_gene_sets is empty")
  }
  
  if (is.null(names(enrichment_gene_sets))) {
    stop("enrichment_gene_sets must be a named list")
  }
  
  if (length(universe_genes) == 0L) {
    stop("universe_genes is empty")
  }
  
  valid_ontologies <- c("BP", "MF", "CC")
  
  ontologies <- unique(
    toupper(ontologies)
  )
  
  invalid_ontologies <- setdiff(
    ontologies,
    valid_ontologies
  )
  
  if (length(invalid_ontologies) > 0L) {
    stop(
      "Invalid GO ontology: ",
      paste(invalid_ontologies, collapse = ", "),
      ". Allowed values are: BP, MF, CC."
    )
  }
  
  if (
    !is.numeric(significance_cutoff) ||
    length(significance_cutoff) != 1L ||
    significance_cutoff <= 0 ||
    significance_cutoff >= 1
  ) {
    stop(
      "significance_cutoff must be a numeric value between 0 and 1"
    )
  }
  
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package 'clusterProfiler' is required")
  }
  
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Package 'org.Hs.eg.db' is required")
  }
  
  # --------------------------------------------------------------------------
  # Prepare enrichment universe
  # --------------------------------------------------------------------------
  
  universe_genes <- unique(
    universe_genes[
      !is.na(universe_genes) &
        universe_genes != ""
    ]
  )
  
  universe_mapping <- clusterProfiler::bitr(
    universe_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  
  universe_mapping <- unique(
    universe_mapping
  )
  
  universe_entrez <- unique(
    universe_mapping$ENTREZID
  )
  
  if (length(universe_entrez) == 0L) {
    stop(
      "No universe genes could be mapped to Entrez Gene identifiers"
    )
  }
  
  # --------------------------------------------------------------------------
  # Map community gene sets once
  # --------------------------------------------------------------------------
  
  community_mappings <- list()
  community_entrez_sets <- list()
  mapping_summary_list <- list()
  
  for (community_name in names(enrichment_gene_sets)) {
    
    community_genes <- unique(
      enrichment_gene_sets[[community_name]][
        !is.na(enrichment_gene_sets[[community_name]]) &
          enrichment_gene_sets[[community_name]] != ""
      ]
    )
    
    community_mapping <- clusterProfiler::bitr(
      community_genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
    
    community_mapping <- unique(
      community_mapping
    )
    
    community_entrez <- unique(
      community_mapping$ENTREZID
    )
    
    community_entrez <- intersect(
      community_entrez,
      universe_entrez
    )
    
    community_mappings[[community_name]] <- community_mapping
    community_entrez_sets[[community_name]] <- community_entrez
    
    mapping_summary_list[[community_name]] <- data.frame(
      community = community_name,
      input_genes = length(community_genes),
      mapped_symbols = length(
        unique(community_mapping$SYMBOL)
      ),
      mapped_entrez = length(community_entrez),
      unmapped_symbols = length(
        setdiff(
          community_genes,
          unique(community_mapping$SYMBOL)
        )
      ),
      mapping_fraction =
        length(unique(community_mapping$SYMBOL)) /
        length(community_genes),
      stringsAsFactors = FALSE
    )
  }
  
  mapping_summary <- do.call(
    rbind,
    mapping_summary_list
  )
  
  rownames(mapping_summary) <- NULL
  
  # --------------------------------------------------------------------------
  # Perform enrichment independently for each requested ontology
  # --------------------------------------------------------------------------
  
  results_by_ontology <- list()
  all_combined_results <- list()
  
  for (ontology in ontologies) {
    
    enrichment_results <- list()
    
    for (community_name in names(community_entrez_sets)) {
      
      community_entrez <-
        community_entrez_sets[[community_name]]
      
      if (length(community_entrez) == 0L) {
        
        enrichment_results[[community_name]] <- data.frame()
        
        next
      }
      
      enrichment <- clusterProfiler::enrichGO(
        gene = community_entrez,
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ontology,
        pAdjustMethod = p_adjust_method,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
      )
      
      enrichment_table <- as.data.frame(
        enrichment
      )
      
      if (nrow(enrichment_table) > 0L) {
        
        enrichment_table$community <- community_name
        
        enrichment_table$ontology <- ontology
        
        enrichment_table$significant <- (
          enrichment_table$p.adjust <= significance_cutoff
        )
        
        enrichment_table <- enrichment_table[
          order(
            enrichment_table$p.adjust,
            enrichment_table$pvalue
          ),
          ,
          drop = FALSE
        ]
      }
      
      enrichment_results[[community_name]] <-
        enrichment_table
    }
    
    nonempty_results <- enrichment_results[
      vapply(
        enrichment_results,
        nrow,
        integer(1)
      ) > 0L
    ]
    
    if (length(nonempty_results) > 0L) {
      
      combined_results <- do.call(
        rbind,
        nonempty_results
      )
      
      rownames(combined_results) <- NULL
      
    } else {
      
      combined_results <- data.frame()
    }
    
    results_by_ontology[[ontology]] <- list(
      enrichment_results = enrichment_results,
      combined_results = combined_results
    )
    
    if (nrow(combined_results) > 0L) {
      all_combined_results[[ontology]] <- combined_results
    }
  }
  
  # --------------------------------------------------------------------------
  # Combine results from all requested ontologies
  # --------------------------------------------------------------------------
  
  if (length(all_combined_results) > 0L) {
    
    combined_results <- do.call(
      rbind,
      all_combined_results
    )
    
    rownames(combined_results) <- NULL
    
  } else {
    
    combined_results <- data.frame()
  }
  
  # --------------------------------------------------------------------------
  # Return
  # --------------------------------------------------------------------------
  
  return(
    list(
      ontologies = ontologies,
      results_by_ontology = results_by_ontology,
      combined_results = combined_results,
      mapping_summary = mapping_summary,
      universe_mapping = universe_mapping
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_leiden_go_enrichment
# ------------------------------------------------------------------------------

# Summarize Gene Ontology enrichment results across Leiden communities.
#
# This function:
# - works with BP, MF, or CC enrichment results
# - preserves all Leiden communities
# - counts total and significant GO terms per community
# - calculates the proportion of significant GO terms
# - identifies the GO term with the lowest adjusted p-value
# - does not apply additional filtering or top-N selection
#
# ------------------------------------------------------------------------------

summarize_leiden_go_enrichment <- function(
    combined_results,
    mapping_summary,
    ontology
) {
  
  ontology <- toupper(ontology)
  
  valid_ontologies <- c("BP", "MF", "CC")
  
  if (!ontology %in% valid_ontologies) {
    stop(
      "ontology must be one of: ",
      paste(valid_ontologies, collapse = ", ")
    )
  }
  
  required_result_columns <- c(
    "community",
    "ID",
    "Description",
    "p.adjust",
    "significant"
  )
  
  if (!all(required_result_columns %in% colnames(combined_results))) {
    stop(
      "combined_results must contain: ",
      paste(required_result_columns, collapse = ", ")
    )
  }
  
  if (!"community" %in% colnames(mapping_summary)) {
    stop("mapping_summary must contain a 'community' column")
  }
  
  communities <- sort(
    unique(mapping_summary$community)
  )
  
  summary_list <- lapply(
    communities,
    function(comm) {
      
      community_results <- combined_results[
        combined_results$community == comm,
        ,
        drop = FALSE
      ]
      
      n_terms <- nrow(community_results)
      
      n_significant <- sum(
        community_results$significant,
        na.rm = TRUE
      )
      
      proportion_significant <- if (n_terms > 0L) {
        n_significant / n_terms
      } else {
        0
      }
      
      valid_p <- which(
        !is.na(community_results$p.adjust)
      )
      
      if (length(valid_p) > 0L) {
        
        best_index <- valid_p[
          which.min(
            community_results$p.adjust[valid_p]
          )
        ]
        
        best_go_id <- community_results$ID[best_index]
        best_go_description <- community_results$Description[best_index]
        best_p_adjust <- community_results$p.adjust[best_index]
        
      } else {
        
        best_go_id <- NA_character_
        best_go_description <- NA_character_
        best_p_adjust <- NA_real_
      }
      
      data.frame(
        ontology = ontology,
        community = comm,
        total_go_terms = n_terms,
        significant_go_terms = n_significant,
        proportion_significant = proportion_significant,
        best_go_id = best_go_id,
        best_go_description = best_go_description,
        best_p_adjust = best_p_adjust,
        stringsAsFactors = FALSE
      )
    }
  )
  
  summary_table <- do.call(
    rbind,
    summary_list
  )
  
  rownames(summary_table) <- NULL
  
  return(summary_table)
}

# ------------------------------------------------------------------------------
# FUNCTION: get_go_background_size
# ------------------------------------------------------------------------------

# Extract the effective GO background size from enrichment results.
#
# clusterProfiler reports the ontology-specific enrichment background through
# the denominator of the BgRatio column. Because the effective annotation
# background may differ among GO ontologies, this function retrieves that
# value directly from the enrichment output rather than assuming that the
# same background is shared by BP, MF, and CC.
#
# Parameters:
# - enrichment_results: data frame containing GO enrichment results for one
#   ontology
#
# Returns:
# - an integer representing the effective background size
# - NA_integer_ if no valid BgRatio value is available
#
# ------------------------------------------------------------------------------

get_go_background_size <- function(
    enrichment_results
) {
  
  if (!is.data.frame(enrichment_results)) {
    stop("enrichment_results must be a data.frame")
  }
  
  if (
    nrow(enrichment_results) == 0L ||
    !"BgRatio" %in% colnames(enrichment_results)
  ) {
    return(NA_integer_)
  }
  
  valid_background_ratios <- enrichment_results$BgRatio[
    !is.na(enrichment_results$BgRatio)
  ]
  
  if (length(valid_background_ratios) == 0L) {
    return(NA_integer_)
  }
  
  background_ratio <- as.character(
    valid_background_ratios[1]
  )
  
  if (!grepl("/", background_ratio, fixed = TRUE)) {
    return(NA_integer_)
  }
  
  as.integer(
    sub(
      ".*/",
      "",
      background_ratio
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: compute_go_coverage
# ------------------------------------------------------------------------------

# Compute gene-mapping and ontology coverage statistics for GO enrichment.
#
# This function:
# - starts from the gene-mapping summary generated during GO enrichment
# - evaluates, for each Leiden community, how many mapped genes are effectively
#   represented in the selected GO ontology
# - extracts the tested-gene denominator from GeneRatio
# - extracts the ontology-specific background size from BgRatio
# - quantifies losses occurring during SYMBOL-to-ENTREZ conversion and
#   ontology annotation
# - works with any GO ontology (BP, MF, or CC)
#
# Parameters:
# - enrichment_results: data frame containing GO enrichment results for one
#   ontology
# - mapping_summary: data frame summarizing SYMBOL-to-ENTREZ mapping for each
#   community
# - ontology: ontology label, typically "BP", "MF", or "CC"
#
# Returns:
# - a data frame containing one row per community with mapping and coverage
#   statistics
#
# ------------------------------------------------------------------------------

compute_go_coverage <- function(
    enrichment_results,
    mapping_summary,
    ontology
) {
  
  if (!is.data.frame(enrichment_results)) {
    stop("enrichment_results must be a data.frame")
  }
  
  if (!is.data.frame(mapping_summary)) {
    stop("mapping_summary must be a data.frame")
  }
  
  ontology <- toupper(ontology)
  
  if (!ontology %in% c("BP", "MF", "CC")) {
    stop("ontology must be one of: BP, MF, CC")
  }
  
  coverage <- mapping_summary
  coverage$ontology <- ontology
  
  coverage$go_tested_genes <- vapply(
    coverage$community,
    function(current_community) {
      
      current_results <- enrichment_results[
        enrichment_results$community == current_community,
        ,
        drop = FALSE
      ]
      
      if (
        nrow(current_results) == 0L ||
        !"GeneRatio" %in% colnames(current_results)
      ) {
        return(NA_integer_)
      }
      
      as.integer(
        sub(
          ".*/",
          "",
          as.character(current_results$GeneRatio[1])
        )
      )
    },
    integer(1)
  )
  
  coverage$go_background_genes <- vapply(
    coverage$community,
    function(current_community) {
      
      current_results <- enrichment_results[
        enrichment_results$community == current_community,
        ,
        drop = FALSE
      ]
      
      if (
        nrow(current_results) == 0L ||
        !"BgRatio" %in% colnames(current_results)
      ) {
        return(NA_integer_)
      }
      
      as.integer(
        sub(
          ".*/",
          "",
          as.character(current_results$BgRatio[1])
        )
      )
    },
    integer(1)
  )
  
  coverage$symbol_to_entrez_loss <- (
    coverage$input_genes -
      coverage$mapped_symbols
  )
  
  coverage$entrez_to_go_loss <- (
    coverage$mapped_entrez -
      coverage$go_tested_genes
  )
  
  coverage$total_loss <- (
    coverage$input_genes -
      coverage$go_tested_genes
  )
  
  coverage$go_coverage_fraction <- (
    coverage$go_tested_genes /
      coverage$input_genes
  )
  
  coverage
}


# ------------------------------------------------------------------------------
# FUNCTION: compute_go_coverage_by_ontology
# ------------------------------------------------------------------------------

# Compute GO coverage statistics for multiple ontologies.
#
# This function:
# - iterates over the requested GO ontologies
# - applies compute_go_coverage() independently to each ontology
# - preserves ontology-specific enrichment results and coverage statistics
# - returns a named list indexed by ontology
#
# Parameters:
# - go_results: named list containing one enrichment result table per ontology
# - mapping_summary: data frame summarizing SYMBOL-to-ENTREZ mapping for each
#   Leiden community
# - ontologies: character vector specifying the GO ontologies to analyse,
#   typically any combination of "BP", "MF", and "CC"
#
# Returns:
# - a named list containing one GO coverage data frame per ontology
#
# ------------------------------------------------------------------------------

compute_go_coverage_by_ontology <- function(
    go_results,
    mapping_summary,
    ontologies = names(go_results)
) {
  
  if (!is.list(go_results)) {
    stop("go_results must be a list")
  }
  
  if (!is.data.frame(mapping_summary)) {
    stop("mapping_summary must be a data.frame")
  }
  
  ontologies <- unique(
    toupper(ontologies)
  )
  
  valid_ontologies <- c(
    "BP",
    "MF",
    "CC"
  )
  
  if (!all(ontologies %in% valid_ontologies)) {
    stop(
      "ontologies must contain only: ",
      paste(valid_ontologies, collapse = ", ")
    )
  }
  
  missing_ontologies <- setdiff(
    ontologies,
    names(go_results)
  )
  
  if (length(missing_ontologies) > 0L) {
    stop(
      "Missing enrichment results for: ",
      paste(missing_ontologies, collapse = ", ")
    )
  }
  
  coverage_results <- lapply(
    ontologies,
    function(current_ontology) {
      
      compute_go_coverage(
        enrichment_results =
          go_results[[current_ontology]],
        mapping_summary =
          mapping_summary,
        ontology =
          current_ontology
      )
    }
  )
  
  names(coverage_results) <- ontologies
  
  coverage_results
}

# ------------------------------------------------------------------------------
# FUNCTION: extract_go_results_by_ontology
# ------------------------------------------------------------------------------

# Extract ontology-specific result tables from a multi-ontology GO analysis.
#
# Parameters:
# - go_enrichment: output generated by perform_leiden_go_enrichment()
#
# Returns:
# - named list containing one combined enrichment table per ontology
#
# ------------------------------------------------------------------------------

extract_go_results_by_ontology <- function(
    go_enrichment
) {
  
  ontologies <- go_enrichment$ontologies
  
  results <- lapply(
    ontologies,
    function(current_ontology) {
      
      go_enrichment$
        results_by_ontology[[current_ontology]]$
        combined_results
    }
  )
  
  names(results) <- ontologies
  
  results
}

# ------------------------------------------------------------------------------
# FUNCTION: summarize_go_enrichment_by_ontology
# ------------------------------------------------------------------------------

# Summarize Leiden GO enrichment results independently for each ontology.
#
# Parameters:
# - go_results: named list containing ontology-specific enrichment tables
# - mapping_summary: community-level SYMBOL-to-ENTREZ mapping summary
#
# Returns:
# - named list containing one enrichment summary per ontology
#
# ------------------------------------------------------------------------------

summarize_go_enrichment_by_ontology <- function(
    go_results,
    mapping_summary
) {
  
  ontologies <- names(go_results)
  
  summaries <- lapply(
    ontologies,
    function(current_ontology) {
      
      summarize_leiden_go_enrichment(
        combined_results =
          go_results[[current_ontology]],
        mapping_summary =
          mapping_summary,
        ontology =
          current_ontology
      )
    }
  )
  
  names(summaries) <- ontologies
  
  summaries
}

# ------------------------------------------------------------------------------
# FUNCTION: extract_gene_community_go_results
# ------------------------------------------------------------------------------

# Extract GO enrichment results for the Leiden community containing a
# selected gene.
#
# This function:
# - identifies the Leiden community containing the selected gene
# - extracts ontology-specific enrichment results for that community
# - retains both complete and significant enrichment results
# - orders significant terms by adjusted p-value
#
# Parameters:
# - gene: gene symbol of interest
# - node_assignment: Leiden node-assignment data frame
# - go_results: named list containing ontology-specific enrichment tables
#
# Returns:
# - list containing the selected community, complete ontology-specific
#   results, and significant ontology-specific results
#
# ------------------------------------------------------------------------------

extract_gene_community_go_results <- function(
    gene,
    node_assignment,
    go_results
) {
  
  gene_community <- unique(
    node_assignment$leiden_community[
      node_assignment$gene == gene
    ]
  )
  
  if (length(gene_community) == 0L) {
    stop("Selected gene was not found in node_assignment")
  }
  
  if (length(gene_community) != 1L) {
    stop("Selected gene is associated with more than one Leiden community")
  }
  
  community_label <- paste0(
    "community_",
    gene_community
  )
  
  community_results <- lapply(
    go_results,
    function(current_results) {
      
      current_results[
        current_results$community == community_label,
        ,
        drop = FALSE
      ]
    }
  )
  
  significant_results <- lapply(
    community_results,
    function(current_results) {
      
      if (nrow(current_results) == 0L) {
        return(current_results)
      }
      
      significant <- current_results[
        current_results$significant,
        ,
        drop = FALSE
      ]
      
      significant[
        order(significant$p.adjust),
        ,
        drop = FALSE
      ]
    }
  )
  
  list(
    gene = gene,
    community = gene_community,
    community_label = community_label,
    results = community_results,
    significant_results = significant_results
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: compute_go_universe_coverage
# ------------------------------------------------------------------------------

# Compute ontology-specific coverage statistics for the enrichment universe.
#
# Parameters:
# - universe_genes: input enrichment-universe gene symbols
# - universe_mapping: SYMBOL-to-ENTREZ mapping table for the universe
# - go_results: named list containing ontology-specific enrichment tables
#
# Returns:
# - data frame containing universe coverage statistics for each ontology
#
# ------------------------------------------------------------------------------

compute_go_universe_coverage <- function(
    universe_genes,
    universe_mapping,
    go_results
) {
  
  ontologies <- names(go_results)
  
  mapped_entrez <- length(
    unique(
      universe_mapping$ENTREZID[
        !is.na(universe_mapping$ENTREZID)
      ]
    )
  )
  
  background_sizes <- vapply(
    ontologies,
    function(current_ontology) {
      
      get_go_background_size(
        go_results[[current_ontology]]
      )
    },
    integer(1)
  )
  
  output <- data.frame(
    ontology = ontologies,
    input_universe_genes = length(universe_genes),
    mapped_universe_entrez = mapped_entrez,
    go_background_genes = as.integer(background_sizes),
    stringsAsFactors = FALSE
  )
  
  output$symbol_to_entrez_loss <- (
    output$input_universe_genes -
      output$mapped_universe_entrez
  )
  
  output$entrez_to_go_loss <- (
    output$mapped_universe_entrez -
      output$go_background_genes
  )
  
  output$total_loss <- (
    output$input_universe_genes -
      output$go_background_genes
  )
  
  output$go_coverage_fraction <- (
    output$go_background_genes /
      output$input_universe_genes
  )
  
  output
}

# ------------------------------------------------------------------------------
# FUNCTION: save_go_enrichment_outputs
# ------------------------------------------------------------------------------

# Save general and ontology-specific GO enrichment outputs.
#
# This function:
# - saves the complete multi-ontology GO enrichment object
# - saves functional gene sets, enrichment summaries, and coverage results
# - saves complete and significant GO results for the community containing
#   the selected gene
# - exports general analysis tables
# - exports ontology-specific BP, MF, or CC result tables
# - preserves separate filename prefixes for the Leiden analysis and the
#   selected-gene community analysis
#
# INPUT:
# - go_enrichment:
#     output generated by perform_leiden_go_enrichment()
#
# - go_results:
#     named list containing complete ontology-specific enrichment tables
#
# - go_summary:
#     named list containing ontology-specific enrichment summaries
#
# - go_coverage:
#     named list containing ontology-specific GO coverage tables
#
# - universe_coverage:
#     ontology-specific enrichment-universe coverage table
#
# - gene_go_results:
#     output generated by extract_gene_community_go_results()
#
# - functional_gene_sets:
#     output generated by prepare_leiden_functional_gene_sets()
#
# - output_dir:
#     directory where RDS objects are saved
#
# - tables_dir:
#     directory where TSV tables are saved
#
# - prefix:
#     filename prefix used for the Leiden-community GO analysis
#
# - gene_prefix:
#     filename prefix used for the selected-gene community GO analysis
#
# OUTPUT:
# - invisibly returns TRUE after all outputs have been successfully saved
#
# ------------------------------------------------------------------------------

save_go_enrichment_outputs <- function(
    go_enrichment,
    go_results,
    go_summary,
    go_coverage,
    universe_coverage,
    gene_go_results,
    functional_gene_sets,
    output_dir,
    tables_dir,
    prefix,
    gene_prefix
) {
  
  # Validate input
  
  if (!is.list(go_enrichment)) {
    stop("go_enrichment must be a list")
  }
  
  if (!is.list(go_results) || is.null(names(go_results))) {
    stop("go_results must be a named list")
  }
  
  if (!is.list(go_summary)) {
    stop("go_summary must be a list")
  }
  
  if (!is.list(go_coverage)) {
    stop("go_coverage must be a list")
  }
  
  if (!is.data.frame(universe_coverage)) {
    stop("universe_coverage must be a data.frame")
  }
  
  if (!is.list(gene_go_results)) {
    stop("gene_go_results must be a list")
  }
  
  if (!is.list(functional_gene_sets)) {
    stop("functional_gene_sets must be a list")
  }
  
  if (
    length(prefix) != 1L ||
    is.na(prefix) ||
    !nzchar(prefix)
  ) {
    stop("prefix must be a single non-empty character string")
  }
  
  if (
    length(gene_prefix) != 1L ||
    is.na(gene_prefix) ||
    !nzchar(gene_prefix)
  ) {
    stop("gene_prefix must be a single non-empty character string")
  }
  
  if (!dir.exists(output_dir)) {
    stop(
      "Output directory does not exist: ",
      output_dir
    )
  }
  
  if (!dir.exists(tables_dir)) {
    stop(
      "Tables directory does not exist: ",
      tables_dir
    )
  }
  
  # Combine ontology-specific coverage and summary tables
  
  coverage_combined <- do.call(
    rbind,
    go_coverage
  )
  
  rownames(coverage_combined) <- NULL
  
  
  summary_combined <- do.call(
    rbind,
    go_summary
  )
  
  rownames(summary_combined) <- NULL
  
  # Save general R objects
  
  saveRDS(
    functional_gene_sets,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_functional_gene_sets.rds"
      )
    )
  )
  
  
  saveRDS(
    go_enrichment,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_go_enrichment.rds"
      )
    )
  )
  
  
  saveRDS(
    go_summary,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_go_summary.rds"
      )
    )
  )
  
  
  saveRDS(
    go_coverage,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_go_coverage.rds"
      )
    )
  )
  
  
  saveRDS(
    gene_go_results$results,
    file.path(
      output_dir,
      paste0(
        gene_prefix,
        "_go_results.rds"
      )
    )
  )
  
  
  saveRDS(
    gene_go_results$significant_results,
    file.path(
      output_dir,
      paste0(
        gene_prefix,
        "_go_significant.rds"
      )
    )
  )
  
  # Save ontology-specific R objects
  
  for (current_ontology in names(go_results)) {
    
    ontology_label <- tolower(
      current_ontology
    )
    
    
    saveRDS(
      go_enrichment$
        results_by_ontology[[current_ontology]],
      file.path(
        output_dir,
        paste0(
          prefix,
          "_go_",
          ontology_label,
          "_enrichment.rds"
        )
      )
    )
    
    
    saveRDS(
      go_summary[[current_ontology]],
      file.path(
        output_dir,
        paste0(
          prefix,
          "_go_",
          ontology_label,
          "_summary.rds"
        )
      )
    )
    
    
    saveRDS(
      go_coverage[[current_ontology]],
      file.path(
        output_dir,
        paste0(
          prefix,
          "_go_",
          ontology_label,
          "_coverage.rds"
        )
      )
    )
  }
  
  # Export general tables

  write.table(
    functional_gene_sets$summary,
    file = file.path(
      tables_dir,
      paste0(
        prefix,
        "_functional_gene_sets_summary.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  write.table(
    go_enrichment$mapping_summary,
    file = file.path(
      tables_dir,
      paste0(
        prefix,
        "_go_mapping_summary.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  write.table(
    coverage_combined,
    file = file.path(
      tables_dir,
      paste0(
        prefix,
        "_go_coverage.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  write.table(
    universe_coverage,
    file = file.path(
      tables_dir,
      paste0(
        prefix,
        "_go_universe_coverage.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  write.table(
    go_enrichment$combined_results,
    file = file.path(
      tables_dir,
      paste0(
        prefix,
        "_go_complete_results.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  write.table(
    summary_combined,
    file = file.path(
      tables_dir,
      paste0(
        prefix,
        "_go_summary.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Export ontology-specific tables
  
  for (current_ontology in names(go_results)) {
    
    ontology_label <- tolower(
      current_ontology
    )
    
    
    write.table(
      go_results[[current_ontology]],
      file = file.path(
        tables_dir,
        paste0(
          prefix,
          "_go_",
          ontology_label,
          "_complete_results.tsv"
        )
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    
    write.table(
      go_coverage[[current_ontology]],
      file = file.path(
        tables_dir,
        paste0(
          prefix,
          "_go_",
          ontology_label,
          "_coverage.tsv"
        )
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    
    write.table(
      go_summary[[current_ontology]],
      file = file.path(
        tables_dir,
        paste0(
          prefix,
          "_go_",
          ontology_label,
          "_summary.tsv"
        )
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    
    write.table(
      gene_go_results$
        significant_results[[current_ontology]],
      file = file.path(
        tables_dir,
        paste0(
          gene_prefix,
          "_go_",
          ontology_label,
          "_significant.tsv"
        )
      ),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# FUNCTION: extract_integrated_grn_downstream_subgraph
# ------------------------------------------------------------------------------

# Extract the complete downstream regulatory subgraph originating from a
# selected gene in the integrated GRN.
#
# This function:
# - uses the complete directed integrated GRN
# - verifies that the selected gene is represented in the network
# - identifies every node reachable from the selected gene by following
#   outgoing regulatory edges
# - does not impose a maximum downstream depth
# - preserves the selected gene as the root node
# - constructs the complete node-induced downstream subgraph
# - preserves all original node and edge attributes
# - calculates the shortest directed distance from the selected gene to every
#   downstream node
# - does not apply score, correlation, gene-type, or distance filtering
#
# IMPORTANT:
# The returned graph is an induced subgraph. Therefore, after identifying all
# nodes reachable downstream from the selected gene, every original GRN edge
# connecting two retained nodes is preserved. The resulting graph is not only
# a shortest-path tree but the complete downstream regulatory subnetwork.
#
# INPUT:
# - graph:
#     complete directed igraph object representing the integrated GRN
#
# - gene:
#     gene symbol used as the root of the downstream analysis
#
# OUTPUT:
# - list containing:
#
#   - graph:
#       complete directed induced subgraph containing the selected gene and all
#       nodes reachable from it through outgoing regulatory interactions
#
#   - nodes:
#       complete node table containing all original vertex attributes together
#       with:
#       - downstream_distance: shortest directed distance from the selected gene
#       - is_selected_gene: whether the node is the root gene
#
#   - edges:
#       complete edge table containing all original edge attributes for
#       interactions connecting nodes retained in the downstream subgraph
#
#   - summary:
#       general summary of the extracted downstream network, including node
#       count, downstream-gene count, edge count, maximum downstream distance,
#       and gene-type composition
#
# ------------------------------------------------------------------------------

extract_integrated_grn_downstream_subgraph <- function(
    graph,
    gene
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for downstream subgraph analysis"
    )
  }
  
  if (!inherits(graph, "igraph")) {
    stop(
      "graph must be an igraph object"
    )
  }
  
  if (!igraph::is_directed(graph)) {
    stop(
      "The integrated GRN must be directed"
    )
  }
  
  if (
    length(gene) != 1L ||
    is.na(gene) ||
    !nzchar(gene)
  ) {
    stop(
      "gene must be a single non-empty gene symbol"
    )
  }
  
  graph_genes <- igraph::V(graph)$name
  
  if (is.null(graph_genes)) {
    stop(
      "Graph vertices must contain gene names"
    )
  }
  
  if (!gene %in% graph_genes) {
    stop(
      "Selected gene is not present in the integrated GRN: ",
      gene
    )
  }
  
  downstream_vertices <- igraph::subcomponent(
    graph,
    v = gene,
    mode = "out"
  )
  
  downstream_genes <- igraph::V(graph)$name[
    downstream_vertices
  ]
  
  downstream_graph <- igraph::induced_subgraph(
    graph,
    vids = downstream_vertices
  )
  
  downstream_distances <- igraph::distances(
    graph,
    v = gene,
    to = downstream_genes,
    mode = "out"
  )
  
  downstream_distances <- as.numeric(
    downstream_distances[1, ]
  )
  
  if (any(!is.finite(downstream_distances))) {
    stop(
      "Non-finite distances were found among downstream-reachable nodes"
    )
  }
  
  nodes <- igraph::as_data_frame(
    downstream_graph,
    what = "vertices"
  )
  
  distance_index <- match(
    nodes$name,
    downstream_genes
  )
  
  if (anyNA(distance_index)) {
    stop(
      "Some downstream graph nodes could not be matched to their distances"
    )
  }
  
  nodes$downstream_distance <- downstream_distances[
    distance_index
  ]
  
  nodes$is_selected_gene <- nodes$name == gene
  
  if ("gene" %in% colnames(nodes)) {
    nodes$gene <- nodes$name
    nodes$name <- NULL
  } else {
    colnames(nodes)[colnames(nodes) == "name"] <- "gene"
  }
  
  nodes <- nodes[
    order(
      nodes$downstream_distance,
      nodes$gene
    ),
    ,
    drop = FALSE
  ]
  
  rownames(nodes) <- NULL
  
  edges <- igraph::as_data_frame(
    downstream_graph,
    what = "edges"
  )
  
  colnames(edges)[
    colnames(edges) == "from"
  ] <- "regulator"
  
  colnames(edges)[
    colnames(edges) == "to"
  ] <- "target"
  
  rownames(edges) <- NULL
  
  maximum_downstream_distance <- max(
    nodes$downstream_distance
  )
  
  node_type_counts <- if ("node_type" %in% colnames(nodes)) {
    table(
      factor(
        nodes$node_type,
        levels = c(
          "TF",
          "lncRNA",
          "PCG"
        )
      )
    )
  } else {
    c(
      TF = NA_integer_,
      lncRNA = NA_integer_,
      PCG = NA_integer_
    )
  }
  
  summary <- data.frame(
    metric = c(
      "selected_gene",
      "total_nodes",
      "downstream_genes",
      "total_edges",
      "maximum_downstream_distance",
      "TF_nodes",
      "lncRNA_nodes",
      "PCG_nodes"
    ),
    value = c(
      gene,
      igraph::vcount(downstream_graph),
      igraph::vcount(downstream_graph) - 1L,
      igraph::ecount(downstream_graph),
      maximum_downstream_distance,
      node_type_counts["TF"],
      node_type_counts["lncRNA"],
      node_type_counts["PCG"]
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      graph = downstream_graph,
      nodes = nodes,
      edges = edges,
      summary = summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: characterize_downstream_subgraph
# ------------------------------------------------------------------------------

# Characterize the structure and composition of a complete gene-centered
# downstream subgraph.
#
# This function:
# - compares downstream-subgraph size with the complete integrated GRN
# - calculates the proportion of global nodes and edges represented
# - summarizes nodes according to downstream distance from the selected gene
# - reports TF, lncRNA, and PCG composition at every downstream distance
# - calculates downstream-subgraph in-degree and out-degree for retained nodes
# - identifies nodes that can further propagate regulatory interactions
# - summarizes positive, negative, zero, and missing Spearman correlations
# - does not filter or remove nodes or edges
#
# INPUT:
# - downstream_subgraph:
#     list generated by extract_integrated_grn_downstream_subgraph(),
#     containing:
#     - graph
#     - nodes
#     - edges
#     - summary
#
# - full_graph:
#     complete directed igraph object representing the integrated GRN
#
# OUTPUT:
# - list containing:
#
#   - summary:
#       global comparison between the downstream subgraph and the complete GRN,
#       including node and edge retention proportions
#
#   - distance_summary:
#       one row per downstream distance containing:
#       - total number of nodes
#       - TF, lncRNA, and PCG counts
#       - nodes with outgoing interactions
#       - nodes without outgoing interactions
#       - mean and maximum downstream-subgraph out-degree
#
#   - edge_sign_summary:
#       number and proportion of downstream edges associated with positive,
#       negative, zero, or missing Spearman correlations
#
# ------------------------------------------------------------------------------

characterize_downstream_subgraph <- function(
    downstream_subgraph,
    full_graph
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for downstream subgraph characterization"
    )
  }
  
  required_objects <- c(
    "graph",
    "nodes",
    "edges",
    "summary"
  )
  
  if (!all(required_objects %in% names(downstream_subgraph))) {
    stop(
      "downstream_subgraph must contain: ",
      paste(required_objects, collapse = ", ")
    )
  }
  
  downstream_graph <- downstream_subgraph$graph
  nodes <- downstream_subgraph$nodes
  edges <- downstream_subgraph$edges
  
  if (!inherits(downstream_graph, "igraph")) {
    stop(
      "downstream_subgraph$graph must be an igraph object"
    )
  }
  
  if (!inherits(full_graph, "igraph")) {
    stop(
      "full_graph must be an igraph object"
    )
  }
  
  if (!igraph::is_directed(downstream_graph)) {
    stop(
      "The downstream subgraph must be directed"
    )
  }
  
  if (!igraph::is_directed(full_graph)) {
    stop(
      "The complete integrated GRN must be directed"
    )
  }
  
  required_node_columns <- c(
    "gene",
    "node_type",
    "downstream_distance"
  )
  
  if (!all(required_node_columns %in% colnames(nodes))) {
    stop(
      "downstream node table must contain: ",
      paste(required_node_columns, collapse = ", ")
    )
  }
  
  if (!"spearman_rho" %in% colnames(edges)) {
    stop(
      "downstream edge table must contain 'spearman_rho'"
    )
  }
  
  if (anyDuplicated(nodes$gene) > 0L) {
    stop(
      "Downstream node table contains duplicated genes"
    )
  }
  
  downstream_nodes <- igraph::vcount(
    downstream_graph
  )
  
  downstream_edges <- igraph::ecount(
    downstream_graph
  )
  
  full_nodes <- igraph::vcount(
    full_graph
  )
  
  full_edges <- igraph::ecount(
    full_graph
  )
  
  selected_gene <- nodes$gene[
    nodes$downstream_distance == 0
  ]
  
  if (length(selected_gene) != 1L) {
    stop(
      "Exactly one downstream root gene must have distance 0"
    )
  }
  
  summary <- data.frame(
    metric = c(
      "selected_gene",
      "global_nodes",
      "downstream_nodes",
      "downstream_node_fraction",
      "global_edges",
      "downstream_edges",
      "downstream_edge_fraction",
      "maximum_downstream_distance"
    ),
    value = c(
      selected_gene,
      full_nodes,
      downstream_nodes,
      downstream_nodes / full_nodes,
      full_edges,
      downstream_edges,
      downstream_edges / full_edges,
      max(nodes$downstream_distance)
    ),
    stringsAsFactors = FALSE
  )
  
  downstream_out_degree <- igraph::degree(
    downstream_graph,
    mode = "out",
    loops = FALSE
  )
  
  downstream_in_degree <- igraph::degree(
    downstream_graph,
    mode = "in",
    loops = FALSE
  )
  
  out_degree_index <- match(
    nodes$gene,
    names(downstream_out_degree)
  )
  
  in_degree_index <- match(
    nodes$gene,
    names(downstream_in_degree)
  )
  
  if (anyNA(out_degree_index) || anyNA(in_degree_index)) {
    stop(
      "Some downstream nodes could not be matched to graph degree statistics"
    )
  }
  
  nodes$downstream_out_degree <- as.numeric(
    downstream_out_degree[out_degree_index]
  )
  
  nodes$downstream_in_degree <- as.numeric(
    downstream_in_degree[in_degree_index]
  )
  
  distances <- sort(
    unique(nodes$downstream_distance)
  )
  
  distance_summary_list <- lapply(
    distances,
    function(current_distance) {
      
      current_nodes <- nodes[
        nodes$downstream_distance == current_distance,
        ,
        drop = FALSE
      ]
      
      data.frame(
        downstream_distance = current_distance,
        n_nodes = nrow(current_nodes),
        n_TF = sum(current_nodes$node_type == "TF"),
        n_lncRNA = sum(current_nodes$node_type == "lncRNA"),
        n_PCG = sum(current_nodes$node_type == "PCG"),
        propagating_nodes = sum(
          current_nodes$downstream_out_degree > 0
        ),
        terminal_nodes = sum(
          current_nodes$downstream_out_degree == 0
        ),
        mean_out_degree = mean(
          current_nodes$downstream_out_degree
        ),
        maximum_out_degree = max(
          current_nodes$downstream_out_degree
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  
  distance_summary <- do.call(
    rbind,
    distance_summary_list
  )
  
  rownames(distance_summary) <- NULL
  
  edge_sign <- ifelse(
    is.na(edges$spearman_rho),
    "missing",
    ifelse(
      edges$spearman_rho > 0,
      "positive",
      ifelse(
        edges$spearman_rho < 0,
        "negative",
        "zero"
      )
    )
  )
  
  expected_signs <- c(
    "positive",
    "negative",
    "zero",
    "missing"
  )
  
  edge_sign_counts <- table(
    factor(
      edge_sign,
      levels = expected_signs
    )
  )
  
  edge_sign_summary <- data.frame(
    correlation_sign = expected_signs,
    n_edges = as.integer(edge_sign_counts),
    proportion = as.integer(edge_sign_counts) / nrow(edges),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      summary = summary,
      distance_summary = distance_summary,
      edge_sign_summary = edge_sign_summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_gene_centered_local_subgraph
# ------------------------------------------------------------------------------

# Prepare a local gene-centered downstream subgraph for community detection.
#
# This function:
# - uses the complete downstream subgraph generated in STEP 17A
# - retains the selected root gene and all nodes whose shortest directed
#   distance from the root is less than or equal to a selected maximum distance
# - constructs the complete induced directed subgraph of the retained nodes
# - preserves all original node and edge attributes
# - preserves the downstream-distance annotation generated in STEP 17A
# - does not apply gene-type, score, correlation, or database-support filtering
#
# IMPORTANT:
# max_distance defines the structural extent of the local gene-centered
# network. In the CASC19 analysis, max_distance = 2 is selected according to
# the propagation pattern characterized in STEP 17B, before the strong network
# expansion observed at downstream distance 3.
#
# INPUT:
# - downstream_subgraph:
#     complete downstream-subgraph object generated by
#     extract_integrated_grn_downstream_subgraph(), containing:
#     - graph
#     - nodes
#     - edges
#     - summary
#
# - max_distance:
#     maximum shortest directed distance from the selected root gene to retain
#     in the local subgraph
#
# OUTPUT:
# - list containing:
#
#   - graph:
#       directed induced igraph object containing all nodes up to the selected
#       downstream distance and every edge connecting retained nodes
#
#   - nodes:
#       complete node table for the retained local subgraph, preserving original
#       annotations and downstream distance
#
#   - edges:
#       complete edge table for interactions connecting retained nodes
#
#   - summary:
#       local-subgraph summary containing root gene, selected maximum distance,
#       retained node and edge counts, and gene-type composition
#
# ------------------------------------------------------------------------------

prepare_gene_centered_local_subgraph <- function(
    downstream_subgraph,
    max_distance
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for gene-centered subgraph preparation"
    )
  }
  
  required_objects <- c(
    "graph",
    "nodes",
    "edges",
    "summary"
  )
  
  if (!all(required_objects %in% names(downstream_subgraph))) {
    stop(
      "downstream_subgraph must contain: ",
      paste(required_objects, collapse = ", ")
    )
  }
  
  graph <- downstream_subgraph$graph
  nodes <- downstream_subgraph$nodes
  
  if (!inherits(graph, "igraph")) {
    stop(
      "downstream_subgraph$graph must be an igraph object"
    )
  }
  
  if (!igraph::is_directed(graph)) {
    stop(
      "The downstream subgraph must be directed"
    )
  }
  
  required_node_columns <- c(
    "gene",
    "node_type",
    "downstream_distance",
    "is_selected_gene"
  )
  
  if (!all(required_node_columns %in% colnames(nodes))) {
    stop(
      "downstream node table must contain: ",
      paste(required_node_columns, collapse = ", ")
    )
  }
  
  if (
    length(max_distance) != 1L ||
    is.na(max_distance) ||
    !is.finite(max_distance) ||
    max_distance < 0 ||
    max_distance != as.integer(max_distance)
  ) {
    stop(
      "max_distance must be a single non-negative integer"
    )
  }
  
  selected_gene <- nodes$gene[
    nodes$is_selected_gene
  ]
  
  if (length(selected_gene) != 1L) {
    stop(
      "Exactly one selected root gene must be present"
    )
  }
  
  local_nodes <- nodes[
    nodes$downstream_distance <= max_distance,
    ,
    drop = FALSE
  ]
  
  if (nrow(local_nodes) == 0L) {
    stop(
      "No nodes were retained in the local downstream subgraph"
    )
  }
  
  if (!selected_gene %in% local_nodes$gene) {
    stop(
      "The selected root gene was not retained"
    )
  }
  
  local_graph <- igraph::induced_subgraph(
    graph,
    vids = local_nodes$gene
  )
  
  graph_node_names <- igraph::V(local_graph)$name
  
  node_index <- match(
    graph_node_names,
    local_nodes$gene
  )
  
  if (anyNA(node_index)) {
    stop(
      "Some local graph nodes could not be matched to the downstream node table"
    )
  }
  
  local_nodes <- local_nodes[
    node_index,
    ,
    drop = FALSE
  ]
  
  rownames(local_nodes) <- NULL
  
  local_edges <- igraph::as_data_frame(
    local_graph,
    what = "edges"
  )
  
  colnames(local_edges)[
    colnames(local_edges) == "from"
  ] <- "regulator"
  
  colnames(local_edges)[
    colnames(local_edges) == "to"
  ] <- "target"
  
  rownames(local_edges) <- NULL
  
  node_type_counts <- table(
    factor(
      local_nodes$node_type,
      levels = c(
        "TF",
        "lncRNA",
        "PCG"
      )
    )
  )
  
  summary <- data.frame(
    metric = c(
      "selected_gene",
      "maximum_retained_distance",
      "retained_nodes",
      "retained_edges",
      "TF_nodes",
      "lncRNA_nodes",
      "PCG_nodes"
    ),
    value = c(
      selected_gene,
      max_distance,
      igraph::vcount(local_graph),
      igraph::ecount(local_graph),
      node_type_counts["TF"],
      node_type_counts["lncRNA"],
      node_type_counts["PCG"]
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      graph = local_graph,
      nodes = local_nodes,
      edges = local_edges,
      summary = summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: characterize_gene_centered_local_nodes
# ------------------------------------------------------------------------------

# Calculate node-level topology within a gene-centered local subgraph.
#
# This function:
# - uses the directed local gene-centered subgraph
# - calculates local in-degree, out-degree, and total degree
# - assigns each node a local network role
# - preserves gene type, downstream distance, and root-gene annotation
# - prepares a node table compatible with
#   assign_integrated_grn_leiden_membership()
# - does not filter or remove nodes
#
# INPUT:
# - local_subgraph:
#     local gene-centered subgraph object generated by
#     prepare_gene_centered_local_subgraph(), containing:
#     - graph
#     - nodes
#     - edges
#     - summary
#
# OUTPUT:
# - data.frame containing one row per local-network node with:
#     - gene
#     - node_type
#     - in_degree
#     - out_degree
#     - total_degree
#     - network_role
#     - downstream_distance
#     - is_selected_gene
#
# ------------------------------------------------------------------------------

characterize_gene_centered_local_nodes <- function(
    local_subgraph
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for local node characterization"
    )
  }
  
  if (!all(c("graph", "nodes") %in% names(local_subgraph))) {
    stop(
      "local_subgraph must contain 'graph' and 'nodes'"
    )
  }
  
  graph <- local_subgraph$graph
  nodes <- local_subgraph$nodes
  
  if (!inherits(graph, "igraph")) {
    stop(
      "local_subgraph$graph must be an igraph object"
    )
  }
  
  if (!igraph::is_directed(graph)) {
    stop(
      "The local gene-centered graph must be directed"
    )
  }
  
  required_columns <- c(
    "gene",
    "node_type",
    "downstream_distance",
    "is_selected_gene"
  )
  
  if (!all(required_columns %in% colnames(nodes))) {
    stop(
      "local node table must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  in_degree <- igraph::degree(
    graph,
    mode = "in",
    loops = FALSE
  )
  
  out_degree <- igraph::degree(
    graph,
    mode = "out",
    loops = FALSE
  )
  
  total_degree <- igraph::degree(
    graph,
    mode = "all",
    loops = FALSE
  )
  
  node_index <- match(
    nodes$gene,
    names(in_degree)
  )
  
  if (anyNA(node_index)) {
    stop(
      "Some local nodes could not be matched to graph vertices"
    )
  }
  
  nodes$in_degree <- as.numeric(
    in_degree[node_index]
  )
  
  nodes$out_degree <- as.numeric(
    out_degree[node_index]
  )
  
  nodes$total_degree <- as.numeric(
    total_degree[node_index]
  )
  
  nodes$network_role <- ifelse(
    nodes$in_degree == 0 & nodes$out_degree == 0,
    "isolated",
    ifelse(
      nodes$in_degree == 0,
      "source_only",
      ifelse(
        nodes$out_degree == 0,
        "target_only",
        "regulator_and_target"
      )
    )
  )
  
  nodes <- nodes[
    ,
    c(
      "gene",
      "node_type",
      "in_degree",
      "out_degree",
      "total_degree",
      "network_role",
      "downstream_distance",
      "is_selected_gene"
    ),
    drop = FALSE
  ]
  
  rownames(nodes) <- NULL
  
  return(nodes)
}

# ------------------------------------------------------------------------------
# FUNCTION: annotate_cancer_related_go_bp
# ------------------------------------------------------------------------------

# Annotate significant GO Biological Process enrichment terms according to
# their potential relevance to cancer-associated biological processes.
#
# This function:
# - uses the complete GO BP enrichment table
# - retains only statistically significant GO terms for cancer-related
#   interpretation
# - preserves the original GO identifiers, descriptions, enrichment statistics,
#   adjusted p-values, gene sets, and community assignments
# - distinguishes direct cancer-related terminology from broader biological
#   processes commonly associated with cancer
# - assigns one or more functional categories according to explicit keyword
#   matching in the GO term description
# - does not modify statistical significance or enrichment results
# - does not infer activation, inhibition, causality, or cancer specificity
#
# CANCER RELEVANCE:
# - direct_cancer:
#     GO descriptions explicitly containing cancer, tumor, tumour, neoplasm,
#     carcinoma, metastasis, or metastatic terminology
#
# - cancer_associated_process:
#     significant GO terms matching predefined cancer-associated biological
#     categories
#
# - not_classified:
#     significant GO terms not matching the predefined categories
#
# INPUT:
# - enrichment_results:
#     complete GO Biological Process enrichment table extracted from the
#     multi-ontology results generated by perform_leiden_go_enrichment(),
#     containing at least:
#     - ID
#     - Description
#     - p.adjust
#     - community
#     - significant
#
# OUTPUT:
# - list containing:
#
#   - terms:
#       all significant GO BP terms with additional columns:
#       - cancer_relevance
#       - cancer_category
#
#   - summary:
#       number of significant, directly cancer-related, cancer-associated,
#       and unclassified GO BP terms for every local Leiden community
#
#   - category_summary:
#       number of annotated GO BP terms for each cancer-associated functional
#       category within every local Leiden community
#
# ------------------------------------------------------------------------------

annotate_cancer_related_go_bp <- function(
    enrichment_results
) {
  
  required_columns <- c(
    "ID",
    "Description",
    "p.adjust",
    "community",
    "significant"
  )
  
  if (!all(required_columns %in% colnames(enrichment_results))) {
    stop(
      "enrichment_results must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  significant_terms <- enrichment_results[
    enrichment_results$significant %in% TRUE,
    ,
    drop = FALSE
  ]
  
  if (nrow(significant_terms) == 0L) {
    stop(
      "No significant GO Biological Process terms are available"
    )
  }
  
  descriptions <- tolower(
    significant_terms$Description
  )
  
  direct_cancer_pattern <- paste(
    c(
      "cancer",
      "tumor",
      "tumour",
      "neoplasm",
      "carcinoma",
      "metastasis",
      "metastatic"
    ),
    collapse = "|"
  )
  
  category_patterns <- list(
    cell_proliferation = paste(
      c(
        "cell proliferation",
        "proliferation"
      ),
      collapse = "|"
    ),
    cell_death_apoptosis = paste(
      c(
        "apopt",
        "cell death"
      ),
      collapse = "|"
    ),
    cell_cycle = paste(
      c(
        "cell cycle",
        "mitotic",
        "mitosis"
      ),
      collapse = "|"
    ),
    dna_damage_repair = paste(
      c(
        "dna repair",
        "dna damage",
        "dna recombination",
        "chromosome segregation"
      ),
      collapse = "|"
    ),
    migration_invasion_adhesion = paste(
      c(
        "cell migration",
        "cell motility",
        "cell adhesion",
        "integrin",
        "extracellular matrix"
      ),
      collapse = "|"
    ),
    angiogenesis_vasculature = paste(
      c(
        "angiogenesis",
        "vasculature",
        "blood vessel",
        "vascular"
      ),
      collapse = "|"
    ),
    immune_inflammatory = paste(
      c(
        "immune",
        "inflamm",
        "cytokine",
        "interleukin",
        "leukocyte",
        "lymphocyte"
      ),
      collapse = "|"
    ),
    signaling = paste(
      c(
        "signal transduction",
        "signaling pathway",
        "signalling pathway",
        "receptor signaling",
        "receptor signalling"
      ),
      collapse = "|"
    )
  )
  
  direct_cancer_match <- grepl(
    direct_cancer_pattern,
    descriptions
  )
  
  category_matches <- lapply(
    category_patterns,
    function(pattern) {
      grepl(
        pattern,
        descriptions
      )
    }
  )
  
  category_matrix <- do.call(
    cbind,
    category_matches
  )
  
  colnames(category_matrix) <- names(
    category_patterns
  )
  
  has_associated_category <- apply(
    category_matrix,
    1,
    any
  )
  
  significant_terms$cancer_relevance <- ifelse(
    direct_cancer_match,
    "direct_cancer",
    ifelse(
      has_associated_category,
      "cancer_associated_process",
      "not_classified"
    )
  )
  
  significant_terms$cancer_category <- vapply(
    seq_len(nrow(significant_terms)),
    function(i) {
      
      matched_categories <- names(category_patterns)[
        category_matrix[i, ]
      ]
      
      if (direct_cancer_match[i]) {
        matched_categories <- c(
          "direct_cancer",
          matched_categories
        )
      }
      
      if (length(matched_categories) == 0L) {
        return("not_classified")
      }
      
      paste(
        matched_categories,
        collapse = ";"
      )
    },
    character(1)
  )
  
  relevance_levels <- c(
    "direct_cancer",
    "cancer_associated_process",
    "not_classified"
  )
  
  communities <- sort(
    unique(significant_terms$community)
  )
  
  summary_list <- lapply(
    communities,
    function(current_community) {
      
      current_terms <- significant_terms[
        significant_terms$community == current_community,
        ,
        drop = FALSE
      ]
      
      data.frame(
        community = current_community,
        significant_go_bp_terms = nrow(current_terms),
        direct_cancer_terms = sum(
          current_terms$cancer_relevance == "direct_cancer"
        ),
        cancer_associated_terms = sum(
          current_terms$cancer_relevance == "cancer_associated_process"
        ),
        not_classified_terms = sum(
          current_terms$cancer_relevance == "not_classified"
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  
  summary <- do.call(
    rbind,
    summary_list
  )
  
  rownames(summary) <- NULL
  
  category_summary_list <- list()
  result_index <- 1L
  
  for (current_community in communities) {
    
    community_index <- significant_terms$community == current_community
    
    for (current_category in names(category_patterns)) {
      
      category_summary_list[[result_index]] <- data.frame(
        community = current_community,
        cancer_category = current_category,
        n_go_bp_terms = sum(
          community_index &
            category_matrix[, current_category]
        ),
        stringsAsFactors = FALSE
      )
      
      result_index <- result_index + 1L
    }
    
    category_summary_list[[result_index]] <- data.frame(
      community = current_community,
      cancer_category = "direct_cancer",
      n_go_bp_terms = sum(
        community_index &
          direct_cancer_match
      ),
      stringsAsFactors = FALSE
    )
    
    result_index <- result_index + 1L
  }
  
  category_summary <- do.call(
    rbind,
    category_summary_list
  )
  
  rownames(category_summary) <- NULL
  
  significant_terms$cancer_relevance <- factor(
    significant_terms$cancer_relevance,
    levels = relevance_levels
  )
  
  significant_terms <- significant_terms[
    order(
      significant_terms$community,
      significant_terms$cancer_relevance,
      significant_terms$p.adjust
    ),
    ,
    drop = FALSE
  ]
  
  rownames(significant_terms) <- NULL
  
  return(
    list(
      terms = significant_terms,
      summary = summary,
      category_summary = category_summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: extract_gene_community_cancer_related_go_bp
# ------------------------------------------------------------------------------

# Extract cancer-related GO Biological Process results for the Leiden community
# containing a selected gene.
#
# This function:
# - identifies the Leiden community containing the selected gene
# - extracts all cancer-related GO BP annotation results assigned to that
#   community
# - orders terms according to cancer relevance and adjusted p-value
# - separately retains terms classified as direct cancer-related or
#   cancer-associated biological processes
# - does not modify the original GO enrichment statistics or classifications
#
# INPUT:
# - gene:
#     gene symbol used to identify the corresponding Leiden community
#
# - node_assignment:
#     Leiden node-assignment table containing at least:
#     - gene
#     - leiden_community
#
# - cancer_related_results:
#     GO BP annotation table generated by annotate_cancer_related_go_bp(),
#     containing at least:
#     - ID
#     - Description
#     - GeneRatio
#     - p.adjust
#     - Count
#     - community
#     - cancer_relevance
#     - cancer_category
#
# OUTPUT:
# - list containing:
#
#   - gene:
#       selected gene
#
#   - community:
#       Leiden community containing the selected gene
#
#   - community_name:
#       community identifier in the format used by the GO enrichment results
#
#   - terms:
#       all significant GO BP terms assigned to the selected-gene community
#
#   - cancer_associated_terms:
#       subset containing only direct cancer-related and broader
#       cancer-associated biological processes
#
# ------------------------------------------------------------------------------

extract_gene_community_cancer_related_go_bp <- function(
    gene,
    node_assignment,
    cancer_related_results
) {
  
  required_node_columns <- c(
    "gene",
    "leiden_community"
  )
  
  if (!all(required_node_columns %in% colnames(node_assignment))) {
    stop(
      "node_assignment must contain: ",
      paste(required_node_columns, collapse = ", ")
    )
  }
  
  required_result_columns <- c(
    "ID",
    "Description",
    "GeneRatio",
    "p.adjust",
    "Count",
    "community",
    "cancer_relevance",
    "cancer_category"
  )
  
  if (!all(required_result_columns %in% colnames(cancer_related_results))) {
    stop(
      "cancer_related_results must contain: ",
      paste(required_result_columns, collapse = ", ")
    )
  }
  
  gene_community <- unique(
    node_assignment$leiden_community[
      node_assignment$gene == gene
    ]
  )
  
  if (length(gene_community) == 0L) {
    stop(
      "Gene not found in Leiden node assignments: ",
      gene
    )
  }
  
  if (
    length(gene_community) != 1L ||
    is.na(gene_community)
  ) {
    stop(
      "The selected gene must belong to exactly one Leiden community"
    )
  }
  
  community_name <- paste0(
    "community_",
    gene_community
  )
  
  community_terms <- cancer_related_results[
    cancer_related_results$community == community_name,
    ,
    drop = FALSE
  ]
  
  community_terms <- community_terms[
    order(
      community_terms$cancer_relevance,
      community_terms$p.adjust
    ),
    ,
    drop = FALSE
  ]
  
  rownames(community_terms) <- NULL
  
  cancer_associated_terms <- community_terms[
    as.character(
      community_terms$cancer_relevance
    ) != "not_classified",
    ,
    drop = FALSE
  ]
  
  rownames(cancer_associated_terms) <- NULL
  
  return(
    list(
      gene = gene,
      community = gene_community,
      community_name = community_name,
      terms = community_terms,
      cancer_associated_terms = cancer_associated_terms
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: save_cancer_related_go_bp_outputs
# ------------------------------------------------------------------------------

# Save cancer-related GO Biological Process annotation results.
#
# This function:
# - saves the complete cancer-related GO BP annotation object
# - exports the complete annotated significant-term table
# - exports community-level and functional-category summaries
# - exports the GO BP results corresponding to the community containing the
#   selected gene
# - preserves separate filename prefixes for the complete local analysis and
#   the selected-gene community analysis
#
# INPUT:
# - cancer_related_results:
#     output generated by annotate_cancer_related_go_bp(), containing:
#     - terms
#     - summary
#     - category_summary
#
# - gene_community_results:
#     output generated by extract_gene_community_cancer_related_go_bp()
#
# - output_dir:
#     directory where the RDS object is saved
#
# - tables_dir:
#     directory where TSV tables are saved
#
# - prefix:
#     filename prefix used for the complete local cancer-related GO BP analysis
#
# - gene_prefix:
#     filename prefix used for the selected-gene community table
#
# OUTPUT:
# - invisibly returns TRUE after all outputs have been successfully saved
#
# ------------------------------------------------------------------------------

save_cancer_related_go_bp_outputs <- function(
    cancer_related_results,
    gene_community_results,
    output_dir,
    tables_dir,
    prefix,
    gene_prefix
) {
  
  required_result_elements <- c(
    "terms",
    "summary",
    "category_summary"
  )
  
  if (
    !is.list(cancer_related_results) ||
    !all(required_result_elements %in% names(cancer_related_results))
  ) {
    stop(
      "cancer_related_results must contain: ",
      paste(required_result_elements, collapse = ", ")
    )
  }
  
  required_gene_elements <- c(
    "gene",
    "community",
    "community_name",
    "terms",
    "cancer_associated_terms"
  )
  
  if (
    !is.list(gene_community_results) ||
    !all(required_gene_elements %in% names(gene_community_results))
  ) {
    stop(
      "gene_community_results must contain: ",
      paste(required_gene_elements, collapse = ", ")
    )
  }
  
  if (
    length(prefix) != 1L ||
    is.na(prefix) ||
    !nzchar(prefix)
  ) {
    stop(
      "prefix must be a single non-empty character string"
    )
  }
  
  if (
    length(gene_prefix) != 1L ||
    is.na(gene_prefix) ||
    !nzchar(gene_prefix)
  ) {
    stop(
      "gene_prefix must be a single non-empty character string"
    )
  }
  
  if (!dir.exists(output_dir)) {
    stop(
      "Output directory does not exist: ",
      output_dir
    )
  }
  
  if (!dir.exists(tables_dir)) {
    stop(
      "Tables directory does not exist: ",
      tables_dir
    )
  }
  
  saveRDS(
    cancer_related_results,
    file.path(
      output_dir,
      paste0(
        prefix,
        "_cancer_related_go_bp.rds"
      )
    )
  )
  
  write.table(
    cancer_related_results$terms,
    file = file.path(
      tables_dir,
      paste0(
        prefix,
        "_cancer_related_go_bp_terms.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    cancer_related_results$summary,
    file = file.path(
      tables_dir,
      paste0(
        prefix,
        "_cancer_related_go_bp_summary.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    cancer_related_results$category_summary,
    file = file.path(
      tables_dir,
      paste0(
        prefix,
        "_cancer_related_go_bp_category_summary.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    gene_community_results$terms,
    file = file.path(
      tables_dir,
      paste0(
        gene_prefix,
        "_cancer_related_go_bp.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# FUNCTION: visualize_gene_centered_downstream_subgraph
# ------------------------------------------------------------------------------

# Visualize the structure of a local gene-centered downstream subgraph.
#
# This function:
# - uses the directed local subgraph generated in STEP 17C
# - preserves all retained nodes and edges
# - calculates a reproducible force-directed layout
# - distinguishes nodes according to downstream distance using both color
#   and shape
# - explicitly highlights and labels the selected root gene
# - suppresses labels for all other genes to preserve readability
# - reduces edge opacity to improve visualization of dense networks
# - does not use Leiden communities, hub status, or functional annotations
#
# INPUT:
# - local_subgraph:
#     local gene-centered subgraph object generated by
#     prepare_gene_centered_local_subgraph(), containing:
#     - graph
#     - nodes
#     - edges
#     - summary
#
# - selected_gene:
#     gene symbol of the root gene to highlight
#
# - seed:
#     random seed used for reproducible graph layout
#
# OUTPUT:
# - list containing:
#
#   - plot:
#       ggraph object representing the local downstream network
#
#   - summary:
#       list containing:
#       - network: general network and visualization statistics
#       - distance_counts: number of nodes at each downstream distance
#
# ------------------------------------------------------------------------------

visualize_gene_centered_downstream_subgraph <- function(
    local_subgraph,
    selected_gene,
    seed = 1
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for network visualization"
    )
  }
  
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop(
      "Package 'ggraph' is required for network visualization"
    )
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for network visualization"
    )
  }
  
  if (!all(c("graph", "nodes") %in% names(local_subgraph))) {
    stop(
      "local_subgraph must contain 'graph' and 'nodes'"
    )
  }
  
  graph <- local_subgraph$graph
  nodes <- local_subgraph$nodes
  
  if (!inherits(graph, "igraph")) {
    stop(
      "local_subgraph$graph must be an igraph object"
    )
  }
  
  if (!igraph::is_directed(graph)) {
    stop(
      "The local gene-centered graph must be directed"
    )
  }
  
  required_columns <- c(
    "gene",
    "downstream_distance",
    "is_selected_gene"
  )
  
  if (!all(required_columns %in% colnames(nodes))) {
    stop(
      "local node table must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (
    length(selected_gene) != 1L ||
    is.na(selected_gene) ||
    !nzchar(selected_gene)
  ) {
    stop(
      "selected_gene must be a single non-empty gene symbol"
    )
  }
  
  if (!selected_gene %in% nodes$gene) {
    stop(
      "selected_gene is not present in the local subgraph"
    )
  }
  
  node_index <- match(
    igraph::V(graph)$name,
    nodes$gene
  )
  
  if (anyNA(node_index)) {
    stop(
      "Some graph vertices could not be matched to the local node table"
    )
  }
  
  graph_nodes <- nodes[
    node_index,
    ,
    drop = FALSE
  ]
  
  igraph::V(graph)$downstream_distance <- factor(
    graph_nodes$downstream_distance,
    levels = sort(
      unique(graph_nodes$downstream_distance)
    )
  )
  
  igraph::V(graph)$selected_gene <- (
    igraph::V(graph)$name == selected_gene
  )
  
  set.seed(seed)
  
  plot <- ggraph::ggraph(
    graph,
    layout = "fr"
  ) +
    ggraph::geom_edge_link(
      ggplot2::aes(
        alpha = combined_score
      ),
      edge_width = 0.1,
      show.legend = FALSE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(
        color = downstream_distance,
        shape = downstream_distance,
        size = selected_gene
      ),
      alpha = 0.85
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(
        label = ifelse(
          selected_gene,
          name,
          ""
        )
      ),
      size = 4,
      fontface = "bold",
      vjust = -1
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "0" = "#D55E00",
        "1" = "#0072B2",
        "2" = "#999999"
      ),
      name = "Downstream distance"
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        "0" = 16,
        "1" = 17,
        "2" = 15
      ),
      name = "Downstream distance"
    ) +
    ggplot2::scale_size_manual(
      values = c(
        `FALSE` = 1,
        `TRUE` = 6
      ),
      guide = "none"
    ) +
    ggraph::scale_edge_alpha(
      range = c(
        0.03,
        0.015
      )
    ) +
    ggplot2::labs(
      title = paste0(
        selected_gene,
        " Local Downstream Regulatory Network"
      ),
      subtitle = paste0(
        "Downstream distance <= ",
        max(graph_nodes$downstream_distance)
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      legend.position = "right"
    )
  
  distance_counts <- as.data.frame(
    table(
      graph_nodes$downstream_distance
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(distance_counts) <- c(
    "downstream_distance",
    "n_nodes"
  )
  
  network_summary <- data.frame(
    metric = c(
      "selected_gene",
      "visualized_nodes",
      "visualized_edges",
      "maximum_downstream_distance",
      "layout_seed"
    ),
    value = c(
      selected_gene,
      igraph::vcount(graph),
      igraph::ecount(graph),
      max(
        as.numeric(
          as.character(
            graph_nodes$downstream_distance
          )
        )
      ),
      seed
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      plot = plot,
      summary = list(
        network = network_summary,
        distance_counts = distance_counts
      )
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: visualize_gene_centered_leiden_communities
# ------------------------------------------------------------------------------

# Visualize Leiden community organization within a local gene-centered
# downstream regulatory network.
#
# This function:
# - uses the directed local subgraph generated in STEP 17C
# - preserves all retained nodes and edges
# - integrates the Leiden community assignments obtained in STEP 17C
# - calculates a reproducible force-directed layout
# - colors nodes according to Leiden community
# - explicitly highlights and labels the selected gene
# - preserves regulatory interactions with reduced edge opacity
# - does not apply additional node or edge filtering
#
# INPUT:
# - local_subgraph:
#     local gene-centered subgraph object generated by
#     prepare_gene_centered_local_subgraph()
#
# - node_assignment:
#     node table generated by
#     assign_integrated_grn_leiden_membership(), containing the
#     leiden_community assignment for every retained node
#
# - selected_gene:
#     gene symbol of the root gene to highlight
#
# - seed:
#     random seed used for reproducible graph layout
#
# OUTPUT:
# - list containing:
#
#   - plot:
#       ggraph object representing Leiden community organization
#
#   - summary:
#       data.frame reporting network size, number of communities,
#       selected gene, and selected-gene community
#
# ------------------------------------------------------------------------------

visualize_gene_centered_leiden_communities <- function(
    local_subgraph,
    node_assignment,
    selected_gene,
    seed = 1
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for network visualization"
    )
  }
  
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop(
      "Package 'ggraph' is required for network visualization"
    )
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for network visualization"
    )
  }
  
  if (!all(c("graph", "nodes") %in% names(local_subgraph))) {
    stop(
      "local_subgraph must contain 'graph' and 'nodes'"
    )
  }
  
  graph <- local_subgraph$graph
  
  if (!inherits(graph, "igraph")) {
    stop(
      "local_subgraph$graph must be an igraph object"
    )
  }
  
  if (!igraph::is_directed(graph)) {
    stop(
      "The local gene-centered graph must be directed"
    )
  }
  
  required_columns <- c(
    "gene",
    "leiden_community"
  )
  
  if (!all(required_columns %in% colnames(node_assignment))) {
    stop(
      "node_assignment must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (
    length(selected_gene) != 1L ||
    is.na(selected_gene) ||
    !nzchar(selected_gene)
  ) {
    stop(
      "selected_gene must be a single non-empty gene symbol"
    )
  }
  
  if (!selected_gene %in% node_assignment$gene) {
    stop(
      "selected_gene is not present in node_assignment"
    )
  }
  
  node_index <- match(
    igraph::V(graph)$name,
    node_assignment$gene
  )
  
  if (anyNA(node_index)) {
    stop(
      "Some graph vertices could not be matched to node_assignment"
    )
  }
  
  graph_nodes <- node_assignment[
    node_index,
    ,
    drop = FALSE
  ]
  
  if (anyNA(graph_nodes$leiden_community)) {
    stop(
      "Some graph vertices do not have a Leiden community assignment"
    )
  }
  
  igraph::V(graph)$leiden_community <- factor(
    graph_nodes$leiden_community,
    levels = sort(
      unique(graph_nodes$leiden_community)
    )
  )
  
  igraph::V(graph)$selected_gene <- (
    igraph::V(graph)$name == selected_gene
  )
  
  selected_community <- as.character(
    graph_nodes$leiden_community[
      graph_nodes$gene == selected_gene
    ]
  )
  
  set.seed(seed)
  
  plot <- ggraph::ggraph(
    graph,
    layout = "fr"
  ) +
    ggraph::geom_edge_link(
      ggplot2::aes(
        alpha = combined_score
      ),
      edge_width = 0.1,
      show.legend = FALSE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(
        color = leiden_community,
        size = selected_gene
      ),
      alpha = 0.85
    ) +
    ggraph::geom_node_point(
      data = function(x) {
        x[x$selected_gene, , drop = FALSE]
      },
      shape = 21,
      size = 6,
      stroke = 1.2,
      fill = "white"
    ) +
    ggraph::geom_node_point(
      data = function(x) {
        x[x$selected_gene, , drop = FALSE]
      },
      ggplot2::aes(
        color = leiden_community
      ),
      size = 4
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(
        label = ifelse(
          selected_gene,
          name,
          ""
        )
      ),
      size = 4,
      fontface = "bold",
      vjust = -1.2
    ) +
    ggplot2::scale_size_manual(
      values = c(
        `FALSE` = 1,
        `TRUE` = 4
      ),
      guide = "none"
    ) +
    ggraph::scale_edge_alpha(
      range = c(
        0.03,
        0.15
      )
    ) +
    ggplot2::labs(
      title = paste0(
        selected_gene,
        " Local Leiden Community Structure"
      ),
      subtitle = paste0(
        selected_gene,
        " belongs to Leiden community ",
        selected_community
      ),
      color = "Leiden community"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      legend.position = "right"
    )
  
  community_counts <- as.data.frame(
    table(
      graph_nodes$leiden_community
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(community_counts) <- c(
    "leiden_community",
    "n_nodes"
  )
  
  network_summary <- data.frame(
    metric = c(
      "selected_gene",
      "selected_gene_community",
      "visualized_nodes",
      "visualized_edges",
      "number_of_communities",
      "layout_seed"
    ),
    value = c(
      selected_gene,
      selected_community,
      igraph::vcount(graph),
      igraph::ecount(graph),
      length(
        unique(
          graph_nodes$leiden_community
        )
      ),
      seed
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      plot = plot,
      summary = list(
        network = network_summary,
        community_counts = community_counts
      )
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: visualize_gene_centered_hubs_tf_lncrna
# ------------------------------------------------------------------------------

# Visualize TF, lncRNA, and community-specific hub organization within a local
# gene-centered downstream regulatory network.
#
# This function:
# - uses the directed local subgraph generated in STEP 17C
# - preserves all retained nodes and edges
# - draws PCG nodes first as a light background layer
# - draws TF and lncRNA nodes above the PCG background
# - draws all community-specific hubs above standard nodes
# - explicitly highlights the selected gene as the top node layer
# - labels the selected gene, all regulatory hubs, and global hubs that are TFs
# - keeps other global and target hubs visible through increased node size
#   without labeling them
# - uses a reproducible force-directed layout
# - does not apply additional filtering to nodes or edges
#
# INPUT:
# - local_subgraph:
#     local gene-centered subgraph object generated by
#     prepare_gene_centered_local_subgraph()
#
# - node_assignment:
#     node table containing:
#     - gene
#     - node_type
#     - leiden_community
#
# - community_hubs:
#     hub-analysis object generated by identify_leiden_community_hubs(),
#     containing:
#     - regulatory_hubs
#     - target_hubs
#     - global_hubs
#     - thresholds
#
# - selected_gene:
#     gene symbol of the root gene to highlight
#
# - seed:
#     random seed used for reproducible graph layout
#
# OUTPUT:
# - list containing:
#
#   - plot:
#       ggraph object showing gene types and local hub organization
#
#   - summary:
#       data.frame reporting network size, gene-type composition, hub counts,
#       and number of labeled genes
#
#   - labeled_nodes:
#       table containing only the genes labeled in the visualization together
#       with node type, Leiden community, and hub annotations
#
# ------------------------------------------------------------------------------

visualize_gene_centered_hubs_tf_lncrna <- function(
    local_subgraph,
    node_assignment,
    community_hubs,
    selected_gene,
    seed = 1
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for network visualization"
    )
  }
  
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop(
      "Package 'ggraph' is required for network visualization"
    )
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for network visualization"
    )
  }
  
  if (!all(c("graph", "nodes") %in% names(local_subgraph))) {
    stop(
      "local_subgraph must contain 'graph' and 'nodes'"
    )
  }
  
  graph <- local_subgraph$graph
  
  if (!inherits(graph, "igraph")) {
    stop(
      "local_subgraph$graph must be an igraph object"
    )
  }
  
  if (!igraph::is_directed(graph)) {
    stop(
      "The local gene-centered graph must be directed"
    )
  }
  
  required_node_columns <- c(
    "gene",
    "node_type",
    "leiden_community"
  )
  
  if (!all(required_node_columns %in% colnames(node_assignment))) {
    stop(
      "node_assignment must contain: ",
      paste(required_node_columns, collapse = ", ")
    )
  }
  
  required_hub_objects <- c(
    "regulatory_hubs",
    "target_hubs",
    "global_hubs"
  )
  
  if (!all(required_hub_objects %in% names(community_hubs))) {
    stop(
      "community_hubs must contain: ",
      paste(required_hub_objects, collapse = ", ")
    )
  }
  
  if (
    length(selected_gene) != 1L ||
    is.na(selected_gene) ||
    !nzchar(selected_gene)
  ) {
    stop(
      "selected_gene must be a single non-empty gene symbol"
    )
  }
  
  if (!selected_gene %in% node_assignment$gene) {
    stop(
      "selected_gene is not present in node_assignment"
    )
  }
  
  node_index <- match(
    igraph::V(graph)$name,
    node_assignment$gene
  )
  
  if (anyNA(node_index)) {
    stop(
      "Some graph vertices could not be matched to node_assignment"
    )
  }
  
  graph_nodes <- node_assignment[
    node_index,
    ,
    drop = FALSE
  ]
  
  regulatory_hub_genes <- unique(
    community_hubs$regulatory_hubs$gene
  )
  
  target_hub_genes <- unique(
    community_hubs$target_hubs$gene
  )
  
  global_hub_genes <- unique(
    community_hubs$global_hubs$gene
  )
  
  all_hub_genes <- unique(
    c(
      regulatory_hub_genes,
      target_hub_genes,
      global_hub_genes
    )
  )
  
  graph_nodes$is_regulatory_hub <- (
    graph_nodes$gene %in% regulatory_hub_genes
  )
  
  graph_nodes$is_target_hub <- (
    graph_nodes$gene %in% target_hub_genes
  )
  
  graph_nodes$is_global_hub <- (
    graph_nodes$gene %in% global_hub_genes
  )
  
  graph_nodes$is_any_hub <- (
    graph_nodes$gene %in% all_hub_genes
  )
  
  graph_nodes$is_selected_gene <- (
    graph_nodes$gene == selected_gene
  )
  
  graph_nodes$is_labelled <- (
    graph_nodes$is_selected_gene |
      graph_nodes$is_regulatory_hub |
      (
        graph_nodes$is_global_hub &
          graph_nodes$node_type == "TF"
      )
  )
  
  graph_nodes$node_display_type <- factor(
    graph_nodes$node_type,
    levels = c(
      "PCG",
      "TF",
      "lncRNA"
    )
  )
  
  igraph::V(graph)$node_display_type <- graph_nodes$node_display_type
  igraph::V(graph)$is_any_hub <- graph_nodes$is_any_hub
  igraph::V(graph)$is_selected_gene <- graph_nodes$is_selected_gene
  igraph::V(graph)$is_labelled <- graph_nodes$is_labelled
  
  igraph::V(graph)$label_gene <- ifelse(
    graph_nodes$is_labelled,
    graph_nodes$gene,
    ""
  )
  
  set.seed(seed)
  
  plot <- ggraph::ggraph(
    graph,
    layout = "fr"
  ) +
    
    ggraph::geom_edge_link(
      ggplot2::aes(
        alpha = combined_score
      ),
      edge_width = 0.1,
      show.legend = FALSE
    ) +
    
    ggraph::geom_node_point(
      data = function(x) {
        x[
          x$node_display_type == "PCG" &
            !x$is_any_hub &
            !x$is_selected_gene,
          ,
          drop = FALSE
        ]
      },
      ggplot2::aes(
        color = node_display_type
      ),
      size = 0.8,
      alpha = 0.30
    ) +
    
    ggraph::geom_node_point(
      data = function(x) {
        x[
          x$node_display_type %in% c("TF", "lncRNA") &
            !x$is_any_hub &
            !x$is_selected_gene,
          ,
          drop = FALSE
        ]
      },
      ggplot2::aes(
        color = node_display_type
      ),
      size = 1.8,
      alpha = 0.9
    ) +
    
    ggraph::geom_node_point(
      data = function(x) {
        x[
          x$is_any_hub &
            !x$is_selected_gene,
          ,
          drop = FALSE
        ]
      },
      ggplot2::aes(
        color = node_display_type
      ),
      size = 3.5,
      alpha = 1
    ) +
    
    ggraph::geom_node_point(
      data = function(x) {
        x[
          x$is_selected_gene,
          ,
          drop = FALSE
        ]
      },
      shape = 21,
      size = 7,
      stroke = 1.3,
      fill = "white"
    ) +
    
    ggraph::geom_node_point(
      data = function(x) {
        x[
          x$is_selected_gene,
          ,
          drop = FALSE
        ]
      },
      ggplot2::aes(
        color = node_display_type
      ),
      size = 5,
      alpha = 1
    ) +
    
    ggraph::geom_node_text(
      data = function(x) {
        x[
          x$is_labelled &
            !x$is_selected_gene,
          ,
          drop = FALSE
        ]
      },
      ggplot2::aes(
        label = label_gene
      ),
      size = 3,
      fontface = "bold",
      repel = TRUE
    ) +
    
    ggraph::geom_node_text(
      data = function(x) {
        x[
          x$is_selected_gene,
          ,
          drop = FALSE
        ]
      },
      ggplot2::aes(
        label = label_gene
      ),
      size = 4,
      fontface = "bold",
      vjust = -1.4
    ) +
    
    ggplot2::scale_color_manual(
      values = c(
        "PCG" = "#BDBDBD",
        "TF" = "#0072B2",
        "lncRNA" = "#D55E00"
      ),
      name = "Gene type"
    ) +
    
    ggraph::scale_edge_alpha(
      range = c(
        0.03,
        0.15
      )
    ) +
    
    ggplot2::labs(
      title = paste0(
        selected_gene,
        " Local Regulatory Hubs and Gene Types"
      ),
      subtitle = "TFs, lncRNAs and community-specific hubs"
    ) +
    
    ggplot2::theme_void() +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      legend.position = "right"
    )
  
  labeled_nodes <- graph_nodes[
    graph_nodes$is_labelled,
    c(
      "gene",
      "node_type",
      "leiden_community",
      "is_regulatory_hub",
      "is_target_hub",
      "is_global_hub",
      "is_selected_gene"
    ),
    drop = FALSE
  ]
  
  labeled_nodes <- labeled_nodes[
    order(
      !labeled_nodes$is_selected_gene,
      labeled_nodes$leiden_community,
      labeled_nodes$gene
    ),
    ,
    drop = FALSE
  ]
  
  rownames(labeled_nodes) <- NULL
  
  summary <- data.frame(
    metric = c(
      "selected_gene",
      "visualized_nodes",
      "visualized_edges",
      "TF_nodes",
      "lncRNA_nodes",
      "PCG_nodes",
      "unique_hub_genes",
      "regulatory_hub_genes",
      "target_hub_genes",
      "global_hub_genes",
      "labeled_nodes",
      "layout_seed"
    ),
    value = c(
      selected_gene,
      igraph::vcount(graph),
      igraph::ecount(graph),
      sum(graph_nodes$node_type == "TF"),
      sum(graph_nodes$node_type == "lncRNA"),
      sum(graph_nodes$node_type == "PCG"),
      length(all_hub_genes),
      length(regulatory_hub_genes),
      length(target_hub_genes),
      length(global_hub_genes),
      nrow(labeled_nodes),
      seed
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      plot = plot,
      summary = summary,
      labeled_nodes = labeled_nodes
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: visualize_go_enrichment_by_community
# ------------------------------------------------------------------------------

# Visualize Gene Ontology enrichment results across Leiden communities.
#
# This function:
# - accepts enrichment results for a single GO ontology
# - retains statistically significant terms for visualization
# - selects the most statistically significant terms independently within each
#   Leiden community
# - converts GeneRatio values to numeric ratios for plotting
# - visualizes adjusted p-values and numbers of contributing genes
# - represents lower adjusted p-values using darker colors
# - applies top-N selection only for visualization
# - does not modify the stored enrichment results
#
# INPUT:
# - enrichment_results:
#     complete enrichment table for a single GO ontology extracted from the
#     multi-ontology results generated by perform_leiden_go_enrichment(),
#     containing:
#     - ID
#     - Description
#     - GeneRatio
#     - p.adjust
#     - Count
#     - community
#     - ontology
#     - significant
#
# - ontology:
#     GO ontology to visualize:
#     - "BP" for Biological Process
#     - "MF" for Molecular Function
#     - "CC" for Cellular Component
#
# - top_terms_per_community:
#     maximum number of significant GO terms displayed for each community
#
# OUTPUT:
# - list containing:
#
#   - plot:
#       dot plot showing the top significant GO terms across Leiden communities
#
#   - terms:
#       exact GO terms used in the visualization
#
#   - summary:
#       data.frame reporting numbers of available and visualized terms
#
# ------------------------------------------------------------------------------

visualize_go_enrichment_by_community <- function(
    enrichment_results,
    ontology,
    top_terms_per_community = 10
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for GO enrichment visualization"
    )
  }
  
  valid_ontologies <- c(
    "BP",
    "MF",
    "CC"
  )
  
  if (
    length(ontology) != 1L ||
    is.na(ontology) ||
    !ontology %in% valid_ontologies
  ) {
    stop(
      "ontology must be one of: ",
      paste(valid_ontologies, collapse = ", ")
    )
  }
  
  required_columns <- c(
    "ID",
    "Description",
    "GeneRatio",
    "p.adjust",
    "Count",
    "community",
    "ontology",
    "significant"
  )
  
  if (!all(required_columns %in% colnames(enrichment_results))) {
    stop(
      "enrichment_results must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (
    length(top_terms_per_community) != 1L ||
    is.na(top_terms_per_community) ||
    top_terms_per_community < 1
  ) {
    stop(
      "top_terms_per_community must be a positive integer"
    )
  }
  
  ontology_results <- enrichment_results[
    enrichment_results$ontology == ontology,
    ,
    drop = FALSE
  ]
  
  if (nrow(ontology_results) == 0L) {
    stop(
      "No enrichment results are available for ontology ",
      ontology
    )
  }
  
  significant_results <- ontology_results[
    ontology_results$significant %in% TRUE,
    ,
    drop = FALSE
  ]
  
  if (nrow(significant_results) == 0L) {
    stop(
      "No significant enrichment results are available for ontology ",
      ontology
    )
  }
  
  parse_gene_ratio <- function(x) {
    
    ratio_parts <- strsplit(
      as.character(x),
      "/",
      fixed = TRUE
    )
    
    vapply(
      ratio_parts,
      function(parts) {
        
        if (length(parts) != 2L) {
          return(NA_real_)
        }
        
        numerator <- as.numeric(parts[1])
        denominator <- as.numeric(parts[2])
        
        if (
          is.na(numerator) ||
          is.na(denominator) ||
          denominator == 0
        ) {
          return(NA_real_)
        }
        
        numerator / denominator
      },
      numeric(1)
    )
  }
  
  significant_results$gene_ratio_numeric <- parse_gene_ratio(
    significant_results$GeneRatio
  )
  
  communities <- sort(
    unique(significant_results$community)
  )
  
  community_term_list <- lapply(
    communities,
    function(current_community) {
      
      current_terms <- significant_results[
        significant_results$community == current_community,
        ,
        drop = FALSE
      ]
      
      current_terms <- current_terms[
        order(
          current_terms$p.adjust,
          -current_terms$Count
        ),
        ,
        drop = FALSE
      ]
      
      head(
        current_terms,
        top_terms_per_community
      )
    }
  )
  
  community_terms <- do.call(
    rbind,
    community_term_list
  )
  
  rownames(community_terms) <- NULL
  
  community_terms$community <- factor(
    community_terms$community,
    levels = communities
  )
  
  community_terms$term_label <- paste0(
    community_terms$Description,
    " [",
    community_terms$ID,
    "]"
  )
  
  community_terms$term_label <- factor(
    community_terms$term_label,
    levels = rev(
      unique(
        community_terms$term_label
      )
    )
  )
  
  ontology_name <- switch(
    ontology,
    BP = "Biological Process",
    MF = "Molecular Function",
    CC = "Cellular Component"
  )
  
  enrichment_plot <- ggplot2::ggplot(
    community_terms,
    ggplot2::aes(
      x = gene_ratio_numeric,
      y = term_label,
      size = Count,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.9
    ) +
    ggplot2::scale_color_gradient(
      low = "black",
      high = "lightcyan"
    ) +
    ggplot2::facet_wrap(
      ~ community,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::labs(
      title = paste0(
        "CASC19 Local Community GO ",
        ontology_name,
        " Enrichment"
      ),
      subtitle = paste0(
        "Top ",
        top_terms_per_community,
        " significant terms per community"
      ),
      x = "Gene ratio",
      y = NULL,
      size = "Gene count",
      color = "Adjusted p-value"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 7
      )
    )
  
  summary <- data.frame(
    metric = c(
      "ontology",
      "significant_terms",
      "communities_with_significant_terms",
      "terms_visualized"
    ),
    value = c(
      ontology,
      nrow(significant_results),
      length(communities),
      nrow(community_terms)
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      plot = enrichment_plot,
      terms = community_terms,
      summary = summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: visualize_go_enrichment_selected_community
# ------------------------------------------------------------------------------

# Visualize Gene Ontology enrichment results within a selected Leiden community.
#
# This function:
# - accepts enrichment results for a single GO ontology
# - selects the Leiden community of interest
# - retains statistically significant terms for visualization
# - selects the most statistically significant terms within the community
# - converts GeneRatio values to numeric ratios for plotting
# - visualizes adjusted p-values and numbers of contributing genes
# - represents lower adjusted p-values using darker colors
# - applies top-N selection only for visualization
# - does not modify the stored enrichment results
#
# INPUT:
# - enrichment_results:
#     complete enrichment table for a single GO ontology extracted from the
#     multi-ontology results generated by perform_leiden_go_enrichment(),
#     containing:
#     - ID
#     - Description
#     - GeneRatio
#     - p.adjust
#     - Count
#     - community
#     - ontology
#     - significant
#
# - selected_community:
#     Leiden community to visualize
#
# - ontology:
#     GO ontology to visualize:
#     - "BP" for Biological Process
#     - "MF" for Molecular Function
#     - "CC" for Cellular Component
#
# - top_terms:
#     maximum number of significant GO terms displayed
#
# OUTPUT:
# - list containing:
#
#   - plot:
#       dot plot showing the top significant GO terms in the selected community
#
#   - terms:
#       exact GO terms used in the visualization
#
#   - summary:
#       data.frame reporting numbers of available and visualized terms
#
# ------------------------------------------------------------------------------

visualize_go_enrichment_selected_community <- function(
    enrichment_results,
    selected_community,
    ontology,
    top_terms = 10
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for GO enrichment visualization"
    )
  }
  
  valid_ontologies <- c(
    "BP",
    "MF",
    "CC"
  )
  
  if (
    length(ontology) != 1L ||
    is.na(ontology) ||
    !ontology %in% valid_ontologies
  ) {
    stop(
      "ontology must be one of: ",
      paste(valid_ontologies, collapse = ", ")
    )
  }
  
  required_columns <- c(
    "ID",
    "Description",
    "GeneRatio",
    "p.adjust",
    "Count",
    "community",
    "ontology",
    "significant"
  )
  
  if (!all(required_columns %in% colnames(enrichment_results))) {
    stop(
      "enrichment_results must contain: ",
      paste(required_columns, collapse = ", ")
    )
  }
  
  if (
    length(selected_community) != 1L ||
    is.na(selected_community)
  ) {
    stop(
      "selected_community must contain exactly one Leiden community"
    )
  }
  
  if (
    length(top_terms) != 1L ||
    is.na(top_terms) ||
    top_terms < 1
  ) {
    stop(
      "top_terms must be a positive integer"
    )
  }
  
  selected_community_name <- paste0(
    "community_",
    selected_community
  )
  
  community_results <- enrichment_results[
    enrichment_results$ontology == ontology &
      enrichment_results$community == selected_community_name,
    ,
    drop = FALSE
  ]
  
  if (nrow(community_results) == 0L) {
    stop(
      "No enrichment results are available for ",
      ontology,
      " in ",
      selected_community_name
    )
  }
  
  significant_results <- community_results[
    community_results$significant %in% TRUE,
    ,
    drop = FALSE
  ]
  
  if (nrow(significant_results) == 0L) {
    stop(
      "No significant ",
      ontology,
      " enrichment results are available for ",
      selected_community_name
    )
  }
  
  ratio_parts <- strsplit(
    as.character(significant_results$GeneRatio),
    "/",
    fixed = TRUE
  )
  
  significant_results$gene_ratio_numeric <- vapply(
    ratio_parts,
    function(parts) {
      
      if (length(parts) != 2L) {
        return(NA_real_)
      }
      
      numerator <- as.numeric(parts[1])
      denominator <- as.numeric(parts[2])
      
      if (
        is.na(numerator) ||
        is.na(denominator) ||
        denominator == 0
      ) {
        return(NA_real_)
      }
      
      numerator / denominator
    },
    numeric(1)
  )
  
  significant_results <- significant_results[
    order(
      significant_results$p.adjust,
      -significant_results$Count
    ),
    ,
    drop = FALSE
  ]
  
  community_terms <- head(
    significant_results,
    top_terms
  )
  
  community_terms$term_label <- paste0(
    community_terms$Description,
    " [",
    community_terms$ID,
    "]"
  )
  
  community_terms$term_label <- factor(
    community_terms$term_label,
    levels = rev(
      community_terms$term_label
    )
  )
  
  ontology_name <- switch(
    ontology,
    BP = "Biological Process",
    MF = "Molecular Function",
    CC = "Cellular Component"
  )
  
  enrichment_plot <- ggplot2::ggplot(
    community_terms,
    ggplot2::aes(
      x = gene_ratio_numeric,
      y = term_label,
      size = Count,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.9
    ) +
    ggplot2::scale_color_gradient(
      low = "black",
      high = "lightcyan"
    ) +
    ggplot2::labs(
      title = paste0(
        "GO ",
        ontology_name,
        " Enrichment in the CASC19 Community"
      ),
      subtitle = paste0(
        "Leiden community ",
        selected_community
      ),
      x = "Gene ratio",
      y = NULL,
      size = "Gene count",
      color = "Adjusted p-value"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(
        size = 8
      )
    )
  
  summary <- data.frame(
    metric = c(
      "ontology",
      "selected_community",
      "significant_terms",
      "terms_visualized"
    ),
    value = c(
      ontology,
      selected_community,
      nrow(significant_results),
      nrow(community_terms)
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      plot = enrichment_plot,
      terms = community_terms,
      summary = summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: visualize_gene_centered_functional_enrichment
# ------------------------------------------------------------------------------

# Visualize GO Biological Process enrichment and cancer-associated functional
# results within a local gene-centered regulatory network.
#
# This function:
# - uses the GO Biological Process enrichment results generated in STEP 17E
#   and extracted from the multi-ontology enrichment analysis
# - retains statistically significant terms for visualization
# - selects the most statistically significant terms independently within each
#   Leiden community
# - converts GeneRatio values to numeric ratios for plotting
# - visualizes adjusted p-values and numbers of contributing genes
# - represents lower adjusted p-values using darker colors
# - separately visualizes cancer-associated GO BP terms from the community
#   containing the selected gene
# - applies top-N selection only for visualization
# - does not modify the stored enrichment results
#
# INPUT:
# - enrichment_results:
#     complete GO Biological Process enrichment table extracted from the
#     multi-ontology results generated by perform_leiden_go_enrichment(),
#     containing:
#     - ID
#     - Description
#     - GeneRatio
#     - p.adjust
#     - Count
#     - community
#     - significant
#
# - cancer_related_results:
#     GO BP table generated by annotate_cancer_related_go_bp(), containing:
#     - ID
#     - Description
#     - GeneRatio
#     - p.adjust
#     - Count
#     - community
#     - cancer_relevance
#     - cancer_category
#
# - selected_community:
#     Leiden community containing the selected gene
#
# - top_terms_per_community:
#     maximum number of significant GO BP terms displayed for each community
#
# - top_cancer_terms:
#     maximum number of cancer-associated GO BP terms displayed for the
#     selected-gene community
#
# OUTPUT:
# - list containing:
#
#   - community_plot:
#       dot plot showing the top significant GO BP terms across local
#       Leiden communities
#
#   - cancer_plot:
#       dot plot showing cancer-associated GO BP terms in the community
#       containing the selected gene
#
#   - community_terms:
#       exact GO BP terms used in the community-level visualization
#
#   - cancer_terms:
#       exact cancer-associated GO BP terms used in the selected-community
#       visualization
#
#   - summary:
#       data.frame reporting numbers of available and visualized terms
#
# ------------------------------------------------------------------------------

visualize_gene_centered_functional_enrichment <- function(
    enrichment_results,
    cancer_related_results,
    selected_community,
    top_terms_per_community = 10,
    top_cancer_terms = 20
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for functional enrichment visualization"
    )
  }
  
  required_enrichment_columns <- c(
    "ID",
    "Description",
    "GeneRatio",
    "p.adjust",
    "Count",
    "community",
    "significant"
  )
  
  if (!all(required_enrichment_columns %in% colnames(enrichment_results))) {
    stop(
      "enrichment_results must contain: ",
      paste(required_enrichment_columns, collapse = ", ")
    )
  }
  
  required_cancer_columns <- c(
    "ID",
    "Description",
    "GeneRatio",
    "p.adjust",
    "Count",
    "community",
    "cancer_relevance",
    "cancer_category"
  )
  
  if (!all(required_cancer_columns %in% colnames(cancer_related_results))) {
    stop(
      "cancer_related_results must contain: ",
      paste(required_cancer_columns, collapse = ", ")
    )
  }
  
  if (
    length(top_terms_per_community) != 1L ||
    is.na(top_terms_per_community) ||
    top_terms_per_community < 1
  ) {
    stop(
      "top_terms_per_community must be a positive integer"
    )
  }
  
  if (
    length(top_cancer_terms) != 1L ||
    is.na(top_cancer_terms) ||
    top_cancer_terms < 1
  ) {
    stop(
      "top_cancer_terms must be a positive integer"
    )
  }
  
  parse_gene_ratio <- function(x) {
    
    ratio_parts <- strsplit(
      as.character(x),
      "/",
      fixed = TRUE
    )
    
    vapply(
      ratio_parts,
      function(parts) {
        
        if (length(parts) != 2L) {
          return(NA_real_)
        }
        
        numerator <- as.numeric(parts[1])
        denominator <- as.numeric(parts[2])
        
        if (
          is.na(numerator) ||
          is.na(denominator) ||
          denominator == 0
        ) {
          return(NA_real_)
        }
        
        numerator / denominator
      },
      numeric(1)
    )
  }
  
  significant_results <- enrichment_results[
    enrichment_results$significant %in% TRUE,
    ,
    drop = FALSE
  ]
  
  if (nrow(significant_results) == 0L) {
    stop(
      "No significant GO BP enrichment results are available"
    )
  }
  
  significant_results$gene_ratio_numeric <- parse_gene_ratio(
    significant_results$GeneRatio
  )
  
  communities <- sort(
    unique(significant_results$community)
  )
  
  community_term_list <- lapply(
    communities,
    function(current_community) {
      
      current_terms <- significant_results[
        significant_results$community == current_community,
        ,
        drop = FALSE
      ]
      
      current_terms <- current_terms[
        order(
          current_terms$p.adjust,
          -current_terms$Count
        ),
        ,
        drop = FALSE
      ]
      
      head(
        current_terms,
        top_terms_per_community
      )
    }
  )
  
  community_terms <- do.call(
    rbind,
    community_term_list
  )
  
  rownames(community_terms) <- NULL
  
  community_terms$community <- factor(
    community_terms$community,
    levels = communities
  )
  
  community_terms$term_label <- paste0(
    community_terms$Description,
    " [",
    community_terms$ID,
    "]"
  )
  
  community_terms$term_label <- factor(
    community_terms$term_label,
    levels = rev(
      unique(
        community_terms$term_label
      )
    )
  )
  
  community_plot <- ggplot2::ggplot(
    community_terms,
    ggplot2::aes(
      x = gene_ratio_numeric,
      y = term_label,
      size = Count,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.9
    ) +
    ggplot2::scale_color_gradient(
      low = "black",
      high = "lightcyan"
    ) +
    ggplot2::facet_wrap(
      ~ community,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::labs(
      title = "GO Biological Process Enrichment Across CASC19 Local Communities",
      subtitle = paste0(
        "Top ",
        top_terms_per_community,
        " significant terms per community"
      ),
      x = "Gene ratio",
      y = NULL,
      size = "Gene count",
      color = "Adjusted p-value"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 7
      )
    )
  
  selected_community_name <- paste0(
    "community_",
    selected_community
  )
  
  cancer_terms <- cancer_related_results[
    cancer_related_results$community == selected_community_name &
      cancer_related_results$cancer_relevance != "not_classified",
    ,
    drop = FALSE
  ]
  
  if (nrow(cancer_terms) == 0L) {
    stop(
      "No cancer-associated GO BP terms were found for the selected community"
    )
  }
  
  cancer_terms$gene_ratio_numeric <- parse_gene_ratio(
    cancer_terms$GeneRatio
  )
  
  cancer_terms <- cancer_terms[
    order(
      cancer_terms$p.adjust,
      -cancer_terms$Count
    ),
    ,
    drop = FALSE
  ]
  
  cancer_terms <- head(
    cancer_terms,
    top_cancer_terms
  )
  
  cancer_terms$term_label <- paste0(
    cancer_terms$Description,
    " [",
    cancer_terms$ID,
    "]"
  )
  
  cancer_terms$term_label <- factor(
    cancer_terms$term_label,
    levels = rev(
      cancer_terms$term_label
    )
  )
  
  cancer_plot <- ggplot2::ggplot(
    cancer_terms,
    ggplot2::aes(
      x = gene_ratio_numeric,
      y = term_label,
      size = Count,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point(
      alpha = 0.9
    ) +
    ggplot2::scale_color_gradient(
      low = "black",
      high = "lightcyan"
    ) +
    ggplot2::labs(
      title = "Cancer-Associated GO Biological Processes in the CASC19 Community",
      subtitle = paste0(
        "Leiden community ",
        selected_community
      ),
      x = "Gene ratio",
      y = NULL,
      size = "Gene count",
      color = "Adjusted p-value"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(
        size = 8
      )
    )
  
  summary <- data.frame(
    metric = c(
      "significant_go_bp_terms",
      "communities_with_significant_terms",
      "community_terms_visualized",
      "selected_community",
      "selected_community_cancer_associated_terms",
      "cancer_terms_visualized"
    ),
    value = c(
      nrow(significant_results),
      length(communities),
      nrow(community_terms),
      selected_community,
      sum(
        cancer_related_results$community == selected_community_name &
          cancer_related_results$cancer_relevance != "not_classified"
      ),
      nrow(cancer_terms)
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      community_plot = community_plot,
      cancer_plot = cancer_plot,
      community_terms = community_terms,
      cancer_terms = cancer_terms,
      summary = summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: characterize_selected_community_regulatory_relationships
# ------------------------------------------------------------------------------

# Characterize directed regulatory relationships between the Leiden community
# containing a selected gene and all other communities in a local regulatory
# network.
#
# This function:
# - identifies the Leiden community containing the selected gene
# - maps regulator and target genes to their Leiden communities
# - retains only inter-community edges involving the selected-gene community
# - preserves the original direction of regulatory interactions
# - classifies relationships as outgoing or incoming relative to the
#   selected-gene community
# - summarizes interactions by connected community and regulatory edge class
# - preserves original GRN edge attributes
# - does not apply additional edge filtering or network inference
#
# INPUT:
# - local_subgraph:
#     local gene-centered regulatory subgraph generated by
#     prepare_gene_centered_local_subgraph(), containing:
#     - graph
#     - nodes
#     - edges
#     - summary
#
# - node_assignment:
#     Leiden node-assignment table containing:
#     - gene
#     - leiden_community
#
# - selected_gene:
#     gene whose Leiden community is used as the focal community
#
# OUTPUT:
# - list containing:
#
#   - selected_gene:
#       selected gene
#
#   - selected_community:
#       Leiden community containing the selected gene
#
#   - edges:
#       complete inter-community edge table involving the selected community
#
#   - community_summary:
#       summary of incoming and outgoing regulatory relationships for every
#       connected community
#
#   - edge_type_summary:
#       summary stratified by regulatory interaction class
#
#   - sign_summary:
#       summary stratified by Spearman correlation sign
#
#   - summary:
#       general analysis summary
#
# ------------------------------------------------------------------------------

characterize_selected_community_regulatory_relationships <- function(
    local_subgraph,
    node_assignment,
    selected_gene
) {
  
  if (!is.list(local_subgraph) || !"edges" %in% names(local_subgraph)) {
    stop(
      "local_subgraph must contain an 'edges' table"
    )
  }
  
  required_node_columns <- c(
    "gene",
    "leiden_community"
  )
  
  if (!all(required_node_columns %in% colnames(node_assignment))) {
    stop(
      "node_assignment must contain: ",
      paste(required_node_columns, collapse = ", ")
    )
  }
  
  edges <- local_subgraph$edges
  
  required_edge_columns <- c(
    "regulator",
    "target"
  )
  
  if (!all(required_edge_columns %in% colnames(edges))) {
    stop(
      "local_subgraph$edges must contain: ",
      paste(required_edge_columns, collapse = ", ")
    )
  }
  
  if (anyDuplicated(node_assignment$gene) > 0L) {
    stop(
      "node_assignment contains duplicated gene identifiers"
    )
  }
  
  selected_gene_rows <- node_assignment[
    node_assignment$gene == selected_gene,
    ,
    drop = FALSE
  ]
  
  if (nrow(selected_gene_rows) != 1L) {
    stop(
      "selected_gene must occur exactly once in node_assignment"
    )
  }
  
  selected_community <- selected_gene_rows$leiden_community[1]
  
  regulator_index <- match(
    edges$regulator,
    node_assignment$gene
  )
  
  target_index <- match(
    edges$target,
    node_assignment$gene
  )
  
  if (anyNA(regulator_index) || anyNA(target_index)) {
    stop(
      "Some local GRN edge genes are missing from the Leiden node assignment"
    )
  }
  
  edges$regulator_community <- node_assignment$leiden_community[
    regulator_index
  ]
  
  edges$target_community <- node_assignment$leiden_community[
    target_index
  ]
  
  inter_edges <- edges[
    edges$regulator_community != edges$target_community &
      (
        edges$regulator_community == selected_community |
          edges$target_community == selected_community
      ),
    ,
    drop = FALSE
  ]
  
  if (nrow(inter_edges) == 0L) {
    stop(
      "No inter-community regulatory relationships involve the selected community"
    )
  }
  
  inter_edges$direction <- ifelse(
    inter_edges$regulator_community == selected_community,
    "outgoing",
    "incoming"
  )
  
  inter_edges$connected_community <- ifelse(
    inter_edges$direction == "outgoing",
    inter_edges$target_community,
    inter_edges$regulator_community
  )
  
  inter_edges$selected_community <- selected_community
  
  if ("spearman_rho" %in% colnames(inter_edges)) {
    
    inter_edges$correlation_sign <- ifelse(
      inter_edges$spearman_rho > 0,
      "positive",
      ifelse(
        inter_edges$spearman_rho < 0,
        "negative",
        "zero"
      )
    )
    
  } else {
    
    inter_edges$correlation_sign <- NA_character_
  }
  
  connected_communities <- sort(
    unique(
      inter_edges$connected_community
    )
  )
  
  directions <- c(
    "outgoing",
    "incoming"
  )
  
  community_summary_list <- list()
  
  summary_index <- 1L
  
  for (current_community in connected_communities) {
    
    for (current_direction in directions) {
      
      current_edges <- inter_edges[
        inter_edges$connected_community == current_community &
          inter_edges$direction == current_direction,
        ,
        drop = FALSE
      ]
      
      if (nrow(current_edges) == 0L) {
        next
      }
      
      mean_score <- if (
        "combined_score" %in% colnames(current_edges)
      ) {
        mean(
          current_edges$combined_score,
          na.rm = TRUE
        )
      } else {
        NA_real_
      }
      
      median_score <- if (
        "combined_score" %in% colnames(current_edges)
      ) {
        stats::median(
          current_edges$combined_score,
          na.rm = TRUE
        )
      } else {
        NA_real_
      }
      
      mean_abs_rho <- if (
        "spearman_rho" %in% colnames(current_edges)
      ) {
        mean(
          abs(current_edges$spearman_rho),
          na.rm = TRUE
        )
      } else {
        NA_real_
      }
      
      lnctard_supported_edges <- if (
        "lnctard_supported" %in% colnames(current_edges)
      ) {
        sum(
          current_edges$lnctard_supported %in% TRUE,
          na.rm = TRUE
        )
      } else {
        NA_integer_
      }
      
      community_summary_list[[summary_index]] <- data.frame(
        selected_community = selected_community,
        connected_community = current_community,
        direction = current_direction,
        n_edges = nrow(current_edges),
        mean_combined_score = mean_score,
        median_combined_score = median_score,
        mean_absolute_spearman_rho = mean_abs_rho,
        lnctard_supported_edges = lnctard_supported_edges,
        stringsAsFactors = FALSE
      )
      
      summary_index <- summary_index + 1L
    }
  }
  
  community_summary <- do.call(
    rbind,
    community_summary_list
  )
  
  rownames(community_summary) <- NULL
  
  if ("edge_type" %in% colnames(inter_edges)) {
    
    edge_type_summary <- as.data.frame(
      table(
        connected_community = inter_edges$connected_community,
        direction = inter_edges$direction,
        edge_type = inter_edges$edge_type
      ),
      stringsAsFactors = FALSE
    )
    
    edge_type_summary <- edge_type_summary[
      edge_type_summary$Freq > 0,
      ,
      drop = FALSE
    ]
    
    colnames(edge_type_summary)[
      colnames(edge_type_summary) == "Freq"
    ] <- "n_edges"
    
  } else {
    
    edge_type_summary <- data.frame()
  }
  
  if ("spearman_rho" %in% colnames(inter_edges)) {
    
    sign_summary <- as.data.frame(
      table(
        connected_community = inter_edges$connected_community,
        direction = inter_edges$direction,
        correlation_sign = inter_edges$correlation_sign
      ),
      stringsAsFactors = FALSE
    )
    
    sign_summary <- sign_summary[
      sign_summary$Freq > 0,
      ,
      drop = FALSE
    ]
    
    colnames(sign_summary)[
      colnames(sign_summary) == "Freq"
    ] <- "n_edges"
    
  } else {
    
    sign_summary <- data.frame()
  }
  
  analysis_summary <- data.frame(
    metric = c(
      "selected_gene",
      "selected_community",
      "connected_communities",
      "intercommunity_edges",
      "outgoing_edges",
      "incoming_edges"
    ),
    value = c(
      selected_gene,
      selected_community,
      length(connected_communities),
      nrow(inter_edges),
      sum(inter_edges$direction == "outgoing"),
      sum(inter_edges$direction == "incoming")
    ),
    stringsAsFactors = FALSE
  )
  
  return(
    list(
      selected_gene = selected_gene,
      selected_community = selected_community,
      edges = inter_edges,
      community_summary = community_summary,
      edge_type_summary = edge_type_summary,
      sign_summary = sign_summary,
      summary = analysis_summary
    )
  )
}

# ------------------------------------------------------------------------------
# FUNCTION: visualize_selected_community_regulatory_relationships
# ------------------------------------------------------------------------------

# Visualize directed regulatory relationships between the Leiden community
# containing a selected gene and the other communities in a local regulatory
# network.
#
# This function generates:
#
# - a directed community-level network summarizing regulatory connections
#   involving the selected-gene community
#
# - a regulatory-class bar plot showing the numbers of incoming and outgoing
#   interactions for each connected community and interaction class
#
# No additional GRN filtering is performed.
#
# INPUT:
# - relationship_results:
#     output generated by
#     characterize_selected_community_regulatory_relationships()
#
# OUTPUT:
# - list containing:
#     - community_network_plot
#     - edge_type_plot
#     - community_network_nodes
#     - community_network_edges
#
# ------------------------------------------------------------------------------

visualize_selected_community_regulatory_relationships <- function(
    relationship_results
) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package 'igraph' is required for inter-community visualization"
    )
  }
  
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop(
      "Package 'ggraph' is required for inter-community visualization"
    )
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for inter-community visualization"
    )
  }
  
  required_objects <- c(
    "selected_gene",
    "selected_community",
    "edges",
    "community_summary",
    "edge_type_summary"
  )
  
  if (!all(required_objects %in% names(relationship_results))) {
    stop(
      "relationship_results must contain: ",
      paste(required_objects, collapse = ", ")
    )
  }
  
  selected_gene <- relationship_results$selected_gene
  selected_community <- relationship_results$selected_community
  inter_edges <- relationship_results$edges
  
  pair_summary <- stats::aggregate(
    rep(1L, nrow(inter_edges)),
    by = list(
      regulator_community = inter_edges$regulator_community,
      target_community = inter_edges$target_community
    ),
    FUN = sum
  )
  
  colnames(pair_summary)[
    colnames(pair_summary) == "x"
  ] <- "n_edges"
  
  pair_summary$from <- paste0(
    "community_",
    pair_summary$regulator_community
  )
  
  pair_summary$to <- paste0(
    "community_",
    pair_summary$target_community
  )
  
  community_ids <- sort(
    unique(
      c(
        inter_edges$regulator_community,
        inter_edges$target_community
      )
    )
  )
  
  community_nodes <- data.frame(
    name = paste0(
      "community_",
      community_ids
    ),
    community = community_ids,
    is_selected_community = community_ids == selected_community,
    stringsAsFactors = FALSE
  )
  
  community_nodes$node_label <- ifelse(
    community_nodes$is_selected_community,
    paste0(
      "Community ",
      community_nodes$community,
      "\n",
      selected_gene
    ),
    paste0(
      "Community ",
      community_nodes$community
    )
  )
  
  community_graph <- igraph::graph_from_data_frame(
    d = pair_summary[
      ,
      c(
        "from",
        "to",
        "n_edges"
      )
    ],
    directed = TRUE,
    vertices = community_nodes
  )
  
  selected_node_name <- paste0(
    "community_",
    selected_community
  )
  
  selected_node_index <- match(
    selected_node_name,
    igraph::V(community_graph)$name
  )
  
  if (is.na(selected_node_index)) {
    stop(
      "Selected community was not found in the community-level graph"
    )
  }
  
  community_layout <- igraph::layout_as_star(
    community_graph,
    center = selected_node_index
  )
  
  community_network_plot <- ggraph::ggraph(
    community_graph,
    layout = "manual",
    x = community_layout[, 1],
    y = community_layout[, 2]
  ) +
    ggraph::geom_edge_fan(
      ggplot2::aes(
        width = n_edges,
        alpha = n_edges
      ),
      arrow = grid::arrow(
        length = grid::unit(
          3,
          "mm"
        ),
        type = "closed"
      ),
      start_cap = ggraph::circle(
        7,
        "mm"
      ),
      end_cap = ggraph::circle(
        7,
        "mm"
      ),
      show.legend = TRUE
    ) +
    ggraph::scale_edge_width(
      range = c(
        0.6,
        4
      ),
      name = "Regulatory edges"
    ) +
    ggraph::scale_edge_alpha(
      range = c(
        0.35,
        0.85
      ),
      guide = "none"
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(
        size = is_selected_community,
        shape = is_selected_community
      ),
      fill = "white",
      color = "black",
      stroke = 1.2
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(
        label = node_label
      ),
      fontface = "bold",
      vjust = -1.5,
      size = 4
    ) +
    ggplot2::scale_size_manual(
      values = c(
        `FALSE` = 7,
        `TRUE` = 12
      ),
      guide = "none"
    ) +
    ggplot2::scale_shape_manual(
      values = c(
        `FALSE` = 21,
        `TRUE` = 22
      ),
      guide = "none"
    ) +
    ggplot2::labs(
      title = paste0(
        selected_gene,
        " Community Regulatory Relationships"
      ),
      subtitle = paste0(
        "Directed inter-community relationships involving Leiden community ",
        selected_community
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5
      ),
      legend.position = "right"
    )
  
  edge_type_summary <- relationship_results$edge_type_summary
  
  if (nrow(edge_type_summary) > 0L) {
    
    edge_type_summary$connected_community_label <- paste0(
      "Community ",
      edge_type_summary$connected_community
    )
    
    edge_type_summary$direction <- factor(
      edge_type_summary$direction,
      levels = c(
        "outgoing",
        "incoming"
      )
    )
    
    edge_type_plot <- ggplot2::ggplot(
      edge_type_summary,
      ggplot2::aes(
        x = connected_community_label,
        y = n_edges,
        fill = direction
      )
    ) +
      ggplot2::geom_col(
        position = "dodge"
      ) +
      ggplot2::facet_wrap(
        ~ edge_type,
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = paste0(
          selected_gene,
          " Community Regulatory Interaction Classes"
        ),
        subtitle = paste0(
          "Incoming and outgoing relationships involving Leiden community ",
          selected_community
        ),
        x = "Connected Leiden community",
        y = "Number of regulatory edges",
        fill = "Direction"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5
        ),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
        legend.position = "right"
      )
    
  } else {
    
    edge_type_plot <- NULL
  }
  
  return(
    list(
      community_network_plot = community_network_plot,
      edge_type_plot = edge_type_plot,
      community_network_nodes = community_nodes,
      community_network_edges = pair_summary
    )
  )
}
