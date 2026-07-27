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

  X <- adata$X
  dge <- Matrix::Matrix(X, sparse = TRUE)
  
  rownames(dge) <- adata$obs_names$to_list()
  colnames(dge) <- adata$var_names$to_list()
  
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
  
  # Merge all objects
  merged_obj <- merge(
    x = obj_list[[1]],
    y = obj_list[-1]
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
  obj <- JoinLayers(obj)
  
  logcounts <- GetAssayData(obj, layer = "data")
  
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
  
  result <- scMRMA(
    input = logcounts,
    species = "Hs",
    normalizedData = TRUE
  )
  
  obj$scMRMA_annotation <- result$multiR$annotationResult[colnames(obj), 1]
  obj$cell_type1 <- result$multiR$annotationResult[colnames(obj), 1]
  
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
# - does NOT remove ambiguous genes (kept for downstream logic consistency)
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
  
  lnc  <- gencode[gencode$type == "lncRNA", ]
  prot <- gencode[gencode$type == "protein_coding", ]
  
  genes_obj <- rownames(obj)
  
  lnc_raw  <- intersect(genes_obj, lnc$SYMBOL)
  prot_raw <- intersect(genes_obj, prot$SYMBOL)
  
  ambiguous_genes <- intersect(lnc_raw, prot_raw)
  
  lnc_features <- unique(lnc_raw)
  prot_features <- unique(prot_raw)
  
  cat("\n================ DEBUG ================\n")
  cat("lnc raw:", length(lnc_raw), "\n")
  cat("prot raw:", length(prot_raw), "\n")
  cat("ambiguous:", length(ambiguous_genes), "\n")
  cat("lnc after filter:", length(lnc_features), "\n")
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
compute_gini_dotplot <- function(obj, features_chunk) {
  
  dp <- DotPlot(obj, features = features_chunk)
  df <- dp$data
  
  df <- df[, c("features.plot", "pct.exp")]
  
  out <- aggregate(pct.exp ~ features.plot, df, function(x) Gini(x))
  
  colnames(out) <- c("gene", "gini")
  return(out)
}

compute_gini_pipeline <- function(obj, features, type, chunk_size = 1000) {
  
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
    stat_compare_means(hjust = 0) +
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
  
  avg_expression <- AverageExpression(obj, features = gene)$RNA
  
  dotplot_data <- DotPlot(obj, features = gene)$data
  
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
  
  colnames(out) <- c("celltype", "gini")
  
  return(out)
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_epithelial_lnc_analysis
# ------------------------------------------------------------------------------
# Extract lncRNA genes from GENCODE annotation that are present in the dataset.
#
# This function:
# - filters GENCODE for lncRNA genes
# - keeps only genes present in the Seurat object
# - returns clean gene vector for downstream analysis
#
# INPUT:
# - obj: Seurat object
# - gencode: annotation data.frame
#
# OUTPUT:
# - character vector of lncRNA genes
# ------------------------------------------------------------------------------
prepare_epithelial_lnc_analysis <- function(obj, gencode) {
  
  lnc <- gencode[gencode$type == "lncRNA", ]
  lnc <- lnc[lnc$SYMBOL %in% rownames(obj), ]
  
  unique(lnc$SYMBOL)
}

# ------------------------------------------------------------------------------
# FUNCTION: run_epithelial_markers
# ------------------------------------------------------------------------------
# Perform differential expression analysis for epithelial cells.
#
# This function:
# - ensures correct Seurat layer structure
# - sets cell identity to Cell_type
# - runs FindMarkers for epithelial cells
# - filters significant genes (adjusted p-value < 0.01)
#
# INPUT:
# - obj: Seurat object
# - genes: vector of genes to test
# - epithelial_label: name of epithelial population
#
# OUTPUT:
# - filtered FindMarkers result (significant genes only)
# ------------------------------------------------------------------------------
run_epithelial_markers <- function(obj, genes, epithelial_label = "Epithelial cells") {
  
  obj <- JoinLayers(obj)
  Idents(obj) <- "Cell_type"
  
  M <- FindMarkers(
    object = obj,
    ident.1 = epithelial_label,
    only.pos = TRUE,
    min.pct = 0,
    logfc.threshold = 0,
    features = genes
  )
  
  M <- subset(M, p_val_adj < 0.01)
  
  return(M)
}

# ------------------------------------------------------------------------------
# FUNCTION: filter_dotplot_by_markers
# ------------------------------------------------------------------------------
# Filter DotPlot data using significant marker genes.
#
# This function:
# - subsets DotPlot output for selected genes
# - optionally orders cell types for consistent plotting
#
# INPUT:
# - dotplot_data: Seurat DotPlot output
# - markers: FindMarkers result
# - cell_order: optional factor ordering for cell types
#
# OUTPUT:
# - filtered DotPlot data.frame
# ------------------------------------------------------------------------------
filter_dotplot_by_markers <- function(dotplot_data, markers, cell_order = NULL) {
  
  df <- subset(dotplot_data, features.plot %in% rownames(markers))
  
  if (!is.null(cell_order)) {
    df$id <- factor(df$id, levels = cell_order)
  }
  
  return(df)
}

# ------------------------------------------------------------------------------
# FUNCTION: select_example_gene
# ------------------------------------------------------------------------------
# Select a representative gene for visualization.
#
# Logic:
# - uses a predefined gene if present in dataset
# - otherwise falls back to first gene from marker list
#
# INPUT:
# - obj: Seurat object
# - gene: preferred gene name
# - fallback_genes: vector of candidate genes
#
# OUTPUT:
# - selected gene name
# ------------------------------------------------------------------------------
select_example_gene <- function(obj, gene = "LINC00523", fallback_genes = NULL) {
  
  if (gene %in% rownames(obj[["RNA"]])) {
    return(gene)
  } else {
    return(fallback_genes[1])
  }
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_lnc_boxplot
# ------------------------------------------------------------------------------
# Boxplot of scaled expression across cell types.
#
# INPUT:
# - df: DotPlot-derived dataframe
# - palette: color vector
#
# OUTPUT:
# - ggplot object
# ------------------------------------------------------------------------------
plot_lnc_boxplot <- function(df, palette) {
  
  ggplot(df, aes(x = reorder(id, avg.exp.scaled),
                 y = avg.exp.scaled,
                 fill = id)) +
    geom_boxplot() +
    coord_flip() +
    scale_fill_manual(values = palette) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab("Expression (Scaled)") +
    xlab("Cell Type")
}

# ------------------------------------------------------------------------------
# FUNCTION: prepare_tcga_overlap
# ------------------------------------------------------------------------------
# Prepare TCGA dataset for overlap analysis.
#
# INPUT:
# - tcga_file: path to TCGA table
# - gencode: annotation table
# - marker_genes: epithelial marker genes
#
# OUTPUT:
# - list with filtered TCGA table and gene sets
# ------------------------------------------------------------------------------
prepare_tcga_overlap <- function(tcga_file, gencode, marker_genes) {
  
  COAD <- read.table(tcga_file, header = TRUE)
  
  lnc <- gencode[gencode$type == "lncRNA", ]
  
  COAD <- COAD[COAD$Gene %in% lnc$SYMBOL, ]
  
  x <- list(
    COAD$Gene[COAD$log2FC > 0],
    COAD$Gene[COAD$log2FC < 0],
    marker_genes
  )
  
  names(x) <- c("TCGA-Up", "TCGA-Down", "Epithelial")
  
  list(COAD = COAD, x = x)
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_tcga_barplot
# ------------------------------------------------------------------------------
# Barplot of TCGA differential expression results.
#
# INPUT:
# - COAD: filtered TCGA dataframe
#
# OUTPUT:
# - ggplot object
# ------------------------------------------------------------------------------
plot_tcga_barplot <- function(COAD) {
  
  ggplot(COAD, aes(x = reorder(Gene, log2FC),
                   y = log2FC,
                   fill = log2FC)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient2(
      midpoint = 0,
      low = viridis::viridis(3)[1],
      high = viridis::viridis(3)[3]
    ) +
    theme_bw() +
    xlab("LncRNAs") +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 6)
    )
}

# ------------------------------------------------------------------------------
# FUNCTION: plot_tcga_scatter
# ------------------------------------------------------------------------------
# Scatter plot Tumor vs Normal expression (TCGA).
#
# INPUT:
# - COAD: TCGA dataframe
# - color: label color
#
# OUTPUT:
# - ggplot object
# ------------------------------------------------------------------------------
plot_tcga_scatter <- function(COAD, color) {
  
  ggplot(COAD, aes(x = log(Tumor + 0.001),
                   y = log(Normal + 0.001))) +
    geom_point() +
    theme_bw() +
    ylab("log(Normal Expression)") +
    xlab("log(Tumor Expression)") +
    geom_label_repel(
      nudge_y = 3,
      nudge_x = 0.5,
      color = color,
      fill = "white",
      max.overlaps = 100000,
      size = 2,
      label = ifelse(COAD$Normal == 0, COAD$Gene, "")
    )
}

