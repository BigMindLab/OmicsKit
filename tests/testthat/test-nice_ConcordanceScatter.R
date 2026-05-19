# tests/testthat/test-nice_ConcordanceScatter.R
test_that("nice_ConcordanceScatter returns a ggplot object", {
  de_x <- data.frame(
    gene = c("ESR1", "PGR", "EGFR", "ERBB2", "TP53"),
    logFC = c(2.5, 1.8, -1.6, 0.4, -0.2),
    padj = c(0.001, 0.01, 0.02, 0.50, 0.80)
  )
  
  de_y <- data.frame(
    gene = c("ESR1", "PGR", "EGFR", "ERBB2", "TP53"),
    logFC = c(0.7, 0.4, -0.5, 0.3, -0.1),
    padj = c(0.002, 0.03, 0.04, 0.20, 0.90)
  )
  
  conc <- concordanceDE(
    de_x = de_x,
    de_y = de_y,
    gene_col = "gene",
    logfc_col = "logFC",
    padj_col = "padj",
    padj_threshold = 0.05,
    logfc_threshold = c(1, 0.20)
  )
  
  p <- nice_ConcordanceScatter(
    concordance_result = conc,
    genes_label = c("ESR1", "EGFR"),
    method = "spearman"
  )
  
  expect_s3_class(p, "ggplot")
})
test_that("nice_ConcordanceScatter works without gene labels", {
  de_x <- data.frame(
    gene = c("ESR1", "PGR", "EGFR", "ERBB2", "TP53"),
    logFC = c(2.5, 1.8, -1.6, 0.4, -0.2),
    padj = c(0.001, 0.01, 0.02, 0.50, 0.80)
  )
  
  de_y <- data.frame(
    gene = c("ESR1", "PGR", "EGFR", "ERBB2", "TP53"),
    logFC = c(0.7, 0.4, -0.5, 0.3, -0.1),
    padj = c(0.002, 0.03, 0.04, 0.20, 0.90)
  )
  
  conc <- concordanceDE(
    de_x = de_x,
    de_y = de_y,
    gene_col = "gene",
    logfc_col = "logFC",
    padj_col = "padj",
    padj_threshold = 0.05,
    logfc_threshold = c(1, 0.20)
  )
  
  p <- nice_ConcordanceScatter(conc)
  
  expect_s3_class(p, "ggplot")
})