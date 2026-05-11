testthat::test_that(".resolve_ensembl_host handles current", {
  expect_equal(.resolve_ensembl_host("current"), "https://www.ensembl.org")
})

testthat::test_that(".get_symbol_attribute returns HGNC for human", {
  res <- .get_symbol_attribute(NULL, "hsapiens_gene_ensembl")
  expect_equal(res$attr, "hgnc_symbol")
  expect_equal(res$source, "HGNC")
})

testthat::test_that(".resolve_symbol_column sets symbol for human", {
  df <- data.frame(
    hgnc_symbol = c("", "TP53"),
    external_gene_name = c("TP53", "TP53"),
    stringsAsFactors = FALSE
  )

  out <- .resolve_symbol_column(df, is_human = TRUE, symbol_attr = "hgnc_symbol")
  expect_equal(out$symbol[1], "TP53") # fallback to external_gene_name
  expect_equal(out$symbol[2], "TP53") # hgnc_symbol preferred
})

testthat::test_that(".resolve_symbol_column sets symbol for non-human", {
  df <- data.frame(
    mgi_symbol = c("Actb"),
    stringsAsFactors = FALSE
  )

  out <- .resolve_symbol_column(df, is_human = FALSE, symbol_attr = "mgi_symbol")
  expect_equal(out$symbol[1], "Actb")
})
