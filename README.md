# LANCE

`LANCE` refactors the original `new_pipeline/` colorectal cancer single-cell analysis into an installable R/Bioconductor-style package.

## Scope

The package provides utilities for:

- importing tab-delimited counts or `.h5ad` files into Seurat objects;
- merging and preprocessing single-cell datasets;
- quantifying cell-type specificity with Gini-based summaries; and
- building lncRNA-focused gene regulatory networks from pseudobulk expression profiles.

## Installation

```r
remotes::install_local(".")
```

## Main workflow areas

1. **Data import** with `load_txt_counts_to_seurat()` or `load_h5ad_to_seurat()`.
2. **Seurat preprocessing** with `merge_seurat_objects()`, `prepare_scrmra_input()`, `run_scrmra_annotation()`, and `seurat_preprocess()`.
3. **Gene specificity analysis** with `gene_biotype_matching()`, `compute_gini_pipeline()`, and related plotting helpers.
4. **GRN analysis** from pseudobulk aggregation through GENIE3/Spearman integration.
5. **Integrated network analysis** for lncRNA, TF, and protein-coding interaction classes.

## Documentation

See `vignettes/LANCE.Rmd` for a package-oriented walkthrough derived from `new_pipeline/main.Rmd`.
