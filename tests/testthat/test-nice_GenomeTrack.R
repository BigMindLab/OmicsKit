test_that("nice_GenomeTrack validates inputs and returns tracks", {
  skip_if_not_installed("Gviz")
  skip_if_not_installed("biomaRt")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomeInfoDb")

  host <- .resolve_ensembl_host("113")
  testthat::skip_if_offline(gsub("^https?://", "", host))


  # ---- Minimal annotations (no BioMart queries) -----------------------------
  annotations <- data.frame(
    geneID     = c("ENSG000001", "ENSG000002"),
    symbol     = c("GENE1", "GENE2"),
    chromosome = c("1", "1"),
    gene_start = c(150, 400),
    gene_end   = c(200, 450),
    strand     = c("+", "-"),
    stringsAsFactors = FALSE
  )

  region <- c(chr = "chr1", start = 100, end = 1000)

  # ---- Error when both annotations and gene_ids are provided ----------------
  expect_error(
    nice_GenomeTrack(
      region = region,
      annotations = annotations,
      gene_ids = c("GENE1")
    ),
    "Provide only one"
  )

  # ---- Error on invalid region ---------------------------------------------
  expect_error(
    nice_GenomeTrack(
      region = c(chr = "chr1", start = 200, end = 100),
      annotations = annotations
    ),
    "start must be less than or equal to end"
  )

  # ---- Error when annotations missing required columns ----------------------
  bad_annotations <- data.frame(
    geneID = "ENSG000001",
    symbol = "GENE1"
  )
  expect_error(
    nice_GenomeTrack(
      region = region,
      annotations = bad_annotations
    ),
    "must include columns"
  )

  # ---- Successful run with annotations --------------------------------------
  out_pdf <- tempfile(fileext = ".pdf")
  tracks <- nice_GenomeTrack(
    region = region,
    annotations = annotations,
    show_transcripts = FALSE,
    export_pdf = out_pdf
  )

  expect_type(tracks, "list")
  expect_length(tracks, 2) # GenomeAxis + Genes

  # ---- Error when track_sizes length does not match track_list --------------
  expect_error(
    nice_GenomeTrack(
      region = region,
      annotations = annotations,
      show_transcripts = FALSE,
      track_sizes = c(1) # wrong length
    ),
    "track_sizes"
  )
})
