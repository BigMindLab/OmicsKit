test_that("concordanceDE classifies categories correctly", {
  de_x <- data.frame(
    gene = c("A", "B", "C", "D", "E"),
    logFC = c(2, 2, 2, 0.2, 0.1),
    padj = c(0.01, 0.01, 0.01, 0.80, 0.80)
  )

  de_y <- data.frame(
    gene = c("A", "B", "C", "D", "E"),
    logFC = c(0.5, -0.5, 0.1, -0.5, 0.1),
    padj = c(0.01, 0.01, 0.80, 0.01, 0.80)
  )

  res <- concordanceDE(de_x, de_y, logfc_threshold = c(1, 0.2))
  cats <- setNames(as.character(res$table$category), res$table$gene)

  expect_equal(cats[["A"]], "concordant")
  expect_equal(cats[["B"]], "discordant")
  expect_equal(cats[["C"]], "only_x")
  expect_equal(cats[["D"]], "only_y")
  expect_equal(cats[["E"]], "not_significant")
  expect_equal(res$concordance_score, 0.5)
})
