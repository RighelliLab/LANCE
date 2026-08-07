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

  obj <- Seurat::CreateSeuratObject(
    counts = dge,
    min.cells = min_cells,
    min.features = min_features
  )

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
# Load a single-cell dataset stored in AnnData (.h5ad) format and convert it
# into a Seurat object.
#
# The function:
# - imports gene expression data from an h5ad file
# - converts the expression matrix to a sparse format
# - loads and attaches cell metadata
# - performs a consistency check between cells and metadata
# - calculates mitochondrial gene percentages
# - applies mitochondrial quality-control filtering
# - optionally saves the resulting Seurat object to disk
#
# INPUT:
# - h5ad_file: path to the AnnData (.h5ad) file
# - meta_file: path to the metadata annotation file
# - mt_pattern: regular expression identifying mitochondrial genes
# - mt_threshold: maximum allowed mitochondrial percentage
# - save_path: optional path for saving the Seurat object
#
# OUTPUT:
# - filtered Seurat object ready for downstream analysis
# ------------------------------------------------------------------------------
load_h5ad_to_seurat <- function(
    h5ad_file,
    meta_file,
    mt_pattern = "^MT-",
    mt_threshold = 5,
    save_path = NULL
) {

  anndata <- reticulate::import("anndata")
  adata <- anndata$read_h5ad(h5ad_file)

  # AnnData stores X as cells x genes; Seurat expects genes x cells.
  X <- adata$X
  dge <- Matrix::Matrix(X, sparse = TRUE)
  obs_names <- unlist(adata$obs_names$to_list(), use.names = FALSE)
  var_names <- unlist(adata$var_names$to_list(), use.names = FALSE)

  if (!all(dim(dge) == c(length(obs_names), length(var_names)))) {
    stop("Dimensioni di X, obs_names e var_names non coerenti")
  }

  dge <- t(dge)
  rownames(dge) <- var_names
  colnames(dge) <- obs_names

  if (anyDuplicated(rownames(dge))) {
    stop("I nomi dei geni non sono univoci")
  }
  if (anyDuplicated(colnames(dge))) {
    stop("I barcode delle cellule non sono univoci")
  }

  meta <- read.table(
    meta_file,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
  )

  if (ncol(dge) != nrow(meta)) {
    stop("Mismatch tra cellule e metadata")
  }
  if (anyDuplicated(rownames(meta))) {
    stop("Cell barcodes in metadata are not unique")
  }
  if (!setequal(colnames(dge), rownames(meta))) {
    stop("Cell barcodes and metadata row names do not match")
  }
  meta <- meta[colnames(dge), , drop = FALSE]

  obj <- Seurat::CreateSeuratObject(counts = dge)

  obj <- Seurat::AddMetaData(obj, meta)

  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    obj,
    pattern = mt_pattern
  )

  obj <- subset(obj, subset = percent.mt < mt_threshold)

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
  DefaultAssay(obj) <- "RNA"
  obj <- JoinLayers(obj)
  obj <- NormalizeData(obj)
  Idents(obj) <- "Cell_type"
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
# - excludes ambiguous genes from both returned biotype sets
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

# ------------------------------------------------------------------------------
# FUNCTION: plot_gini_results
# ------------------------------------------------------------------------------
# Visualize Gini coefficient distributions for gene biotypes.
#
# The function:
# - creates boxplot comparing lncRNA vs protein-coding genes
# - adds statistical comparison
# - annotates gene counts
# - applies fixed color palette for reproducibility
#
# INPUT:
# - GiniCoef_clean: cleaned Gini results table
# - lnc_n: number of lncRNA genes
# - prot_n: number of protein-coding genes
#
# OUTPUT:
# - ggplot object
# ------------------------------------------------------------------------------
plot_gini_results <- function(GiniCoef_clean, lnc_n, prot_n) {

  ggplot(GiniCoef_clean, aes(x = type, y = gini, fill = type)) +
    geom_boxplot() +
    ggpubr::stat_compare_means(hjust = 0) +
    theme_bw() +
    scale_fill_manual(values = kelly()[c(12,9)]) +
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