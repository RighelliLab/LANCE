test_that("prepare_grn_gene_sets keeps disjoint lnc and protein-coding sets", {
  expression_matrix <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8),
    nrow = 4,
    dimnames = list(c("LNC1", "GENE1", "AMBIG", "GENE2"), c("S1", "S2"))
  )

  gencode <- data.frame(
    SYMBOL = c("LNC1", "GENE1", "AMBIG", "AMBIG", "GENE2"),
    type = c("lncRNA", "protein_coding", "lncRNA", "protein_coding", "protein_coding"),
    stringsAsFactors = FALSE
  )

  gene_sets <- prepare_grn_gene_sets(expression_matrix, gencode)

  expect_equal(gene_sets$lnc, c("LNC1", "AMBIG"))
  expect_equal(gene_sets$prot, c("GENE1", "GENE2"))
  expect_equal(gene_sets$ambiguous, "AMBIG")
  expect_length(intersect(gene_sets$lnc, gene_sets$prot), 0)
})

test_that("calculate_grn_expression_statistics summarizes expression and detection", {
  expression_matrix <- matrix(
    c(1, 2, 3, 4, 2, 2, 2, 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("LNC1", "GENE1"), c("S1", "S2", "S3", "S4"))
  )
  count_matrix <- matrix(
    c(0, 1, 3, 4, 1, 1, 0, 0),
    nrow = 2,
    byrow = TRUE,
    dimnames = dimnames(expression_matrix)
  )

  stats_tbl <- calculate_grn_expression_statistics(expression_matrix, count_matrix)

  expect_equal(stats_tbl$gene, c("LNC1", "GENE1"))
  expect_equal(stats_tbl$mean_logcpm, c(2.5, 2))
  expect_equal(stats_tbl$n_detected_pseudobulk, c(3, 2))
})

test_that("summarize_grn_correlations counts thresholded pairs", {
  correlation_matrix <- matrix(
    c(1, 0.8, -0.4,
      0.8, 1, 0.2,
      -0.4, 0.2, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  summary_obj <- summarize_grn_correlations(correlation_matrix, thresholds = c(0.3, 0.7))

  expect_equal(summary_obj$summary$value[summary_obj$summary$statistic == "total_unique_pairs"], 3)
  expect_equal(summary_obj$threshold_counts$n_pairs, c(2, 1))
})

test_that("integrate_genie3_spearman drops self edges and preserves scores", {
  weights <- matrix(
    c(0.9, 0.4,
      0.2, 0.8),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("TF1", "TF2"), c("TF1", "GENE1"))
  )

  spearman <- matrix(
    c(1, 0.5, -0.3,
      0.5, 1, 0.7,
      -0.3, 0.7, 1),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("TF1", "TF2", "GENE1"), c("TF1", "TF2", "GENE1"))
  )

  edges <- integrate_genie3_spearman(weights, spearman)

  expect_false(any(edges$regulator == edges$target))
  expect_true(all(c("regulator", "target", "genie3_weight", "spearman_rho") %in% colnames(edges)))
})

test_that("construct_final_grn annotates correlation sign", {
  integrated_edges <- data.frame(
    regulator = c("TF1", "TF1", "TF2"),
    target = c("GENE1", "GENE2", "GENE3"),
    genie3_weight = c(0.9, 0.2, 0.8),
    spearman_rho = c(0.7, -0.8, -0.2),
    stringsAsFactors = FALSE
  )

  final_edges <- construct_final_grn(integrated_edges, genie3_cutoff = 0.5, spearman_threshold = 0.5)

  expect_equal(nrow(final_edges), 1)
  expect_equal(final_edges$correlation_sign, "positive")
})

test_that("prepare_grn_node_table and connectivity calculations stay consistent", {
  grn_edges <- data.frame(
    regulator = c("TF1", "TF1", "TF2"),
    target = c("LNC1", "GENE1", "LNC1"),
    genie3_weight = c(1, 1, 1),
    spearman_rho = c(0.6, 0.5, -0.5),
    correlation_sign = c("positive", "positive", "negative"),
    stringsAsFactors = FALSE
  )
  gene_sets <- list(lnc = "LNC1", prot = c("TF1", "TF2", "GENE1"))

  nodes <- prepare_grn_node_table(grn_edges, gene_sets, tf_regulators = c("TF1", "TF2"))
  connectivity <- calculate_grn_node_connectivity(grn_edges, nodes)

  expect_true(all(c("TF1", "TF2", "LNC1", "GENE1") %in% connectivity$gene))
  expect_equal(connectivity$out_degree[connectivity$gene == "TF1"], 2)
  expect_equal(connectivity$in_degree[connectivity$gene == "LNC1"], 2)
})
