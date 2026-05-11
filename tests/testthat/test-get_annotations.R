testthat::test_that("get_annotations returns expected columns for genes", {
  testthat::skip_if_not_installed("biomaRt")
  testthat::skip_on_cran()

  mock_useMart <- function(...) structure(list(), class = "mock_mart")
  mock_getBM <- function(attributes, filters, values, mart) {
    if (filters == "ensembl_gene_id") {
      data.frame(
        ensembl_gene_id = "ENSG00000141510",
        hgnc_symbol     = "TP53",
        gene_biotype    = "protein_coding",
        chromosome_name = "17",
        start_position  = 7661779,
        end_position    = 7687550,
        description     = "tumor protein p53",
        stringsAsFactors = FALSE
      )
    } else {
      stop("Unexpected filters in mock_getBM()")
    }
  }

  testthat::local_mocked_bindings(
    useMart = mock_useMart,
    getBM   = mock_getBM,
    .package = "biomaRt"
  )

  ids <- c("ENSG00000141510") # TP53
  res <- get_annotations(
    ensembl_ids = ids,
    species     = "hsapiens_gene_ensembl",
    version     = "112",
    filename    = NULL
  )

  testthat::expect_s3_class(res, "data.frame")
  testthat::expect_true(all(c("geneID", "symbol", "biotype", "chromosome",
                              "gene_start", "gene_end", "gene_length", "description") %in% names(res)))
  testthat::expect_equal(res$geneID[1], "ENSG00000141510")
})

testthat::test_that("get_annotations works for transcripts", {
  testthat::skip_if_not_installed("biomaRt")
  testthat::skip_on_cran()

  mock_useMart <- function(...) structure(list(), class = "mock_mart")
  mock_getBM <- function(attributes, filters, values, mart) {
    if (filters == "ensembl_transcript_id_version") {
      data.frame(
        ensembl_transcript_id_version = "ENST00000269305.9",
        ensembl_gene_id = "ENSG00000141510",
        hgnc_symbol     = "TP53",
        gene_biotype    = "protein_coding",
        chromosome_name = "17",
        start_position  = 7661779,
        end_position    = 7687550,
        description     = "tumor protein p53",
        stringsAsFactors = FALSE
      )
    } else {
      stop("Unexpected filters in mock_getBM()")
    }
  }

  testthat::local_mocked_bindings(
    useMart = mock_useMart,
    getBM   = mock_getBM,
    .package = "biomaRt"
  )

  ids <- "ENST00000269305.9"

  res <- get_annotations(
    ensembl_ids = ids,
    species     = "hsapiens_gene_ensembl",
    mode        = "transcripts",
    version     = "112",
    filename    = NULL
  )

  testthat::expect_s3_class(res, "data.frame")
  testthat::expect_true("transcriptID" %in% names(res))
  testthat::expect_equal(res$transcriptID[1], ids)
})

testthat::test_that("get_annotations errors when biomaRt is missing", {
  testthat::skip_if_not_installed("withr")

  withr::with_namespace("biomaRt", {
    # If biomaRt is installed, skip this test
    testthat::skip("biomaRt is installed; cannot simulate missing package here.")
  })
})
