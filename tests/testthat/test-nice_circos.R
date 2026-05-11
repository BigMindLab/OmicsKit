# Tests for nice_circos ----------------------------------------------------
# These tests do not require internet access.

# Helpers -----------------------------------------------------------------

get_circos_internal <- function(name) {
  get(name, envir = environment(nice_circos), inherits = FALSE)
}

toy_cytoband <- function(chromosomes = c("chrA", "chrB"),
                         lengths = c(chrA = 1000, chrB = 1500)) {
  data.frame(
    chr = chromosomes,
    start = 0,
    end = as.numeric(lengths[chromosomes]),
    name = chromosomes,
    gieStain = "gneg",
    stringsAsFactors = FALSE
  )
}

toy_tracks <- function() {
  list(
    RNA = data.frame(
      chr = c("chrA", "chrA", "chrB"),
      start = c(1, 300, 100),
      end = c(100, 450, 250),
      value = c(2, 4, 3),
      name = c("gene1", "gene2", "gene3"),
      stringsAsFactors = FALSE
    ),
    CNV = data.frame(
      chr = c("chrA", "chrB"),
      start = c(1, 1),
      end = c(900, 1200),
      value = c(0.3, -0.2),
      name = c("seg1", "seg2"),
      stringsAsFactors = FALSE
    ),
    Mutation = data.frame(
      chr = c("chrA", "chrB"),
      start = c(500, 700),
      end = c(501, 701),
      value = c(0.45, 0.60),
      name = c("TP53_p.X", "BRCA1_p.Y"),
      stringsAsFactors = FALSE
    )
  )
}

toy_tracks_config <- function() {
  x <- toy_tracks()

  list(
    RNA = list(
      data = x$RNA,
      type = "histogram",
      color = "#7F77DD",
      height = 0.08,
      gradient = FALSE
    ),
    CNV = list(
      data = x$CNV,
      type = "line",
      color = "#1D9E75",
      height = 0.08
    ),
    Mutation = list(
      data = x$Mutation,
      type = "mutation",
      color = "#60A5FA",
      height = 0.06
    )
  )
}


# resolve_genome -----------------------------------------------------------

test_that("resolve_genome validates genome_build input", {
  expect_error(resolve_genome(NULL), "`genome_build`")
  expect_error(resolve_genome(character()), "`genome_build`")
  expect_error(resolve_genome(NA_character_), "`genome_build`")
  expect_error(resolve_genome(""), "`genome_build`")
  expect_error(resolve_genome(c("ce10", "dm6")), "`genome_build`")
})

test_that("resolve_genome resolves ce10 without internet", {
  g <- resolve_genome("ce10")

  expect_true(is_circos_genome(g))
  expect_s3_class(g, "circos_genome")
  expect_equal(g$genome_build, "ce10")
  expect_equal(g$source, "built-in")
  expect_equal(
    g$chromosomes,
    c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX")
  )
  expect_equal(g$plot_type, c("axis", "labels"))

  expect_true(is.data.frame(g$cytoband))
  expect_true(all(c("chr", "start", "end", "name", "gieStain") %in% names(g$cytoband)))
  expect_true(all(g$cytoband$chr %in% g$chromosomes))
  expect_true(all(g$cytoband$end > g$cytoband$start))

  expect_type(g$chr_lengths, "double")
  expect_equal(names(g$chr_lengths), g$chromosomes)
})

test_that("resolve_genome respects chromosome_index order for ce10", {
  g <- resolve_genome("ce10", chromosome_index = c("chrX", "chrI"))

  expect_equal(g$chromosomes, c("chrX", "chrI"))
  expect_equal(unique(g$cytoband$chr), c("chrX", "chrI"))
  expect_equal(names(g$chr_lengths), c("chrX", "chrI"))
})

test_that("resolve_genome resolves dm6 without internet", {
  g <- resolve_genome("dm6", chromosome_index = c("chrX", "chr2L"))

  expect_true(is_circos_genome(g))
  expect_equal(g$genome_build, "dm6")
  expect_equal(g$source, "built-in")
  expect_equal(g$chromosomes, c("chrX", "chr2L"))
})

test_that("resolve_genome rejects bad chromosome_index", {
  expect_error(
    resolve_genome("ce10", chromosome_index = c("chrI", "chrBad")),
    "chrBad"
  )

  expect_error(
    resolve_genome("ce10", chromosome_index = c("chrI", "chrI")),
    "duplicated"
  )
})

test_that("is_circos_genome works", {
  g <- resolve_genome("dm6")

  expect_true(is_circos_genome(g))
  expect_false(is_circos_genome(list()))
  expect_false(is_circos_genome(data.frame()))
})

test_that("print.circos_genome returns object invisibly", {
  g <- resolve_genome("ce10", chromosome_index = c("chrI", "chrX"))

  out <- capture.output(ret <- print(g))

  expect_true(any(grepl("Genome: ce10", out)))
  expect_true(any(grepl("Chromosomes:", out)))
  expect_identical(ret, g)
})


# Cytoband internals -------------------------------------------------------

test_that(".validate_circos_cytoband accepts standard cytoband", {
  validate <- get_circos_internal(".validate_circos_cytoband")

  cb <- toy_cytoband()
  out <- validate(cb)

  expect_true(is.data.frame(out))
  expect_equal(names(out), c("chr", "start", "end", "name", "gieStain"))
  expect_equal(out$chr, cb$chr)
})

test_that(".validate_circos_cytoband accepts unnamed UCSC-style cytoband", {
  validate <- get_circos_internal(".validate_circos_cytoband")

  cb <- toy_cytoband()
  names(cb) <- paste0("V", seq_len(ncol(cb)))

  out <- validate(cb)

  expect_equal(names(out), c("chr", "start", "end", "name", "gieStain"))
  expect_equal(out$chr, c("chrA", "chrB"))
  expect_equal(out$gieStain, c("gneg", "gneg"))
})

test_that(".validate_circos_cytoband accepts common aliases", {
  validate <- get_circos_internal(".validate_circos_cytoband")

  cb <- data.frame(
    chrom = c("chrA", "chrB"),
    chromStart = c(0, 0),
    chromEnd = c(1000, 2000),
    band = c("p", "q"),
    stain = c("gneg", "gpos50"),
    stringsAsFactors = FALSE
  )

  out <- validate(cb)

  expect_equal(names(out), c("chr", "start", "end", "name", "gieStain"))
  expect_equal(out$chr, c("chrA", "chrB"))
  expect_equal(out$start, c(0, 0))
  expect_equal(out$end, c(1000, 2000))
})

test_that(".validate_circos_cytoband rejects malformed data", {
  validate <- get_circos_internal(".validate_circos_cytoband")

  expect_error(validate(list()), "data frame")

  bad_cols <- data.frame(chr = "chrA", start = 0)
  expect_error(validate(bad_cols), "columns")

  bad_coords <- toy_cytoband()
  bad_coords$start[1] <- 100
  bad_coords$end[1] <- 10
  expect_error(validate(bad_coords), "invalid")
})

test_that(".filter_circos_cytoband filters and orders chromosomes", {
  filter_cb <- get_circos_internal(".filter_circos_cytoband")

  cb <- toy_cytoband(chromosomes = c("chrA", "chrB"))
  out <- filter_cb(cb, chromosome_index = c("chrB", "chrA"))

  expect_equal(out$chromosomes, c("chrB", "chrA"))
  expect_equal(unique(out$cytoband$chr), c("chrB", "chrA"))
})

test_that(".synthetic_cytoband creates valid cytoband", {
  synthetic <- get_circos_internal(".synthetic_cytoband")

  out <- synthetic(c(chrA = 100, chrB = 200))

  expect_equal(names(out), c("chr", "start", "end", "name", "gieStain"))
  expect_equal(out$chr, c("chrA", "chrB"))
  expect_equal(out$start, c(0, 0))
  expect_equal(out$end, c(100, 200))
  expect_equal(out$gieStain, c("gneg", "gneg"))

  expect_error(synthetic(c(100, 200)), "named")
  expect_error(synthetic(c(chrA = -1)), "positive")
})


# Track normalization internals ------------------------------------------

test_that(".normalize_circos_track standardizes data frame columns", {
  normalize <- get_circos_internal(".normalize_circos_track")

  df <- data.frame(
    chrom = c("chrA", "chrA"),
    chromStart = c(1, 100),
    chromEnd = c(50, 200),
    score = c(1.2, 2.4),
    id = c("a", "b"),
    stringsAsFactors = FALSE
  )

  out <- normalize(df, track_name = "RNA", track_type = "histogram")

  expect_equal(out$type, "histogram")
  expect_equal(names(out$data), c("chr", "start", "end", "value", "name"))
  expect_equal(out$data$chr, c("chrA", "chrA"))
  expect_equal(out$data$value, c(1.2, 2.4))
  expect_equal(out$data$name, c("a", "b"))
})

test_that(".normalize_circos_track infers track types from names", {
  normalize <- get_circos_internal(".normalize_circos_track")

  base_df <- data.frame(
    chr = "chrA",
    start = 1,
    end = 10,
    value = 1,
    stringsAsFactors = FALSE
  )

  expect_equal(normalize(base_df, "CNV", NULL)$type, "line")
  expect_equal(normalize(base_df, "methylation", NULL)$type, "points")
  expect_equal(normalize(base_df, "somatic_mutations", NULL)$type, "mutation")
  expect_equal(normalize(base_df, "RNA", NULL)$type, "histogram")
})

test_that(".normalize_circos_track handles mutation tracks without value", {
  normalize <- get_circos_internal(".normalize_circos_track")

  df <- data.frame(
    chr = "chrA",
    start = 10,
    end = 11,
    name = "TP53_p.X",
    stringsAsFactors = FALSE
  )

  out <- normalize(df, track_name = "Mutation", track_type = "mutation")

  expect_equal(out$type, "mutation")
  expect_equal(out$data$value, 1)
  expect_equal(out$data$name, "TP53_p.X")
})

test_that(".normalize_circos_track warns for non-mutation tracks without value", {
  normalize <- get_circos_internal(".normalize_circos_track")

  df <- data.frame(
    chr = "chrA",
    start = 10,
    end = 20,
    stringsAsFactors = FALSE
  )

  expect_warning(
    out <- normalize(df, track_name = "RNA", track_type = "histogram"),
    "using value = 1"
  )

  expect_equal(out$data$value, 1)
})

test_that(".normalize_circos_track supports config list values", {
  normalize <- get_circos_internal(".normalize_circos_track")

  df <- data.frame(
    chr = "chrA",
    start = 1,
    end = 10,
    value = 2,
    stringsAsFactors = FALSE
  )

  out <- normalize(
    list(
      data = df,
      type = "histogram",
      color = "#123456",
      height = 0.07,
      gradient = FALSE
    ),
    track_name = "RNA",
    track_type = NULL
  )

  expect_equal(out$type, "histogram")
  expect_equal(out$color, "#123456")
  expect_equal(out$height, 0.07)

  if ("gradient" %in% names(out)) {
    expect_false(out$gradient)
  }
})

test_that(".normalize_circos_track rejects invalid tracks", {
  normalize <- get_circos_internal(".normalize_circos_track")

  bad_coord <- data.frame(
    chr = "chrA",
    start = 20,
    end = 10,
    value = 1,
    stringsAsFactors = FALSE
  )

  expect_error(
    normalize(bad_coord, track_name = "bad", track_type = "histogram"),
    "invalid"
  )

  valid_coord <- data.frame(
    chr = "chrA",
    start = 1,
    end = 10,
    value = 1,
    stringsAsFactors = FALSE
  )

  expect_error(
    normalize(valid_coord, track_name = "bad", track_type = "unknown"),
    "must be one of"
  )
})


# Metadata and links -------------------------------------------------------

test_that(".normalize_circos_metadata accepts named vector metadata", {
  normalize_metadata <- get_circos_internal(".normalize_circos_metadata")

  out <- normalize_metadata(c(Subtype = "Basal", Sample = "Tumor"))

  expect_equal(names(out), c("variable", "group"))
  expect_equal(out$variable, c("Subtype", "Sample"))
  expect_equal(out$group, c("Basal", "Tumor"))
})

test_that(".normalize_circos_metadata accepts genomic metadata", {
  normalize_metadata <- get_circos_internal(".normalize_circos_metadata")

  md <- data.frame(
    chrom = "chrA",
    chromStart = 1,
    chromEnd = 100,
    group = "gain",
    stringsAsFactors = FALSE
  )

  out <- normalize_metadata(md)

  expect_true(all(c("chr", "start", "end", "group", "variable") %in% names(out)))
  expect_equal(out$variable, "metadata")
  expect_equal(out$chr, "chrA")
  expect_equal(out$group, "gain")
})

test_that(".normalize_circos_metadata rejects unsupported multirow metadata", {
  normalize_metadata <- get_circos_internal(".normalize_circos_metadata")

  md <- data.frame(
    subtype = c("A", "B"),
    stage = c("I", "II"),
    stringsAsFactors = FALSE
  )

  expect_error(normalize_metadata(md), "must be")
})

test_that(".circos_validate_link_data validates and filters links", {
  validate_links <- get_circos_internal(".circos_validate_link_data")

  links <- data.frame(
    chr1 = c("chrA", "chrZ"),
    start1 = c(1, 1),
    end1 = c(10, 10),
    chr2 = c("chrB", "chrB"),
    start2 = c(20, 20),
    end2 = c(30, 30),
    stringsAsFactors = FALSE
  )

  expect_warning(
    out <- validate_links(links, chromosomes = c("chrA", "chrB")),
    "Dropping"
  )

  expect_equal(nrow(out), 1)
  expect_equal(out$chr1, "chrA")
})

test_that(".circos_validate_link_data rejects malformed links", {
  validate_links <- get_circos_internal(".circos_validate_link_data")

  bad_cols <- data.frame(
    chr1 = "chrA",
    start1 = 1,
    end1 = 10
  )

  expect_error(validate_links(bad_cols), "must have columns")

  bad_coords <- data.frame(
    chr1 = "chrA",
    start1 = 20,
    end1 = 10,
    chr2 = "chrB",
    start2 = 1,
    end2 = 10,
    stringsAsFactors = FALSE
  )

  expect_error(validate_links(bad_coords), "invalid")
})


# nice_circos argument validation -----------------------------------------

test_that("nice_circos validates top-level arguments", {
  skip_if_not_installed("circlize")

  expect_error(
    nice_circos(data_tracks = data.frame()),
    "`data_tracks` must be a list"
  )

  expect_error(
    nice_circos(
      genome_build = "ce10",
      data_tracks = list(
        A = data.frame(
          chr = "chrI",
          start = 1,
          end = 10,
          value = 1,
          stringsAsFactors = FALSE
        )
      ),
      track_types = c("line", "points")
    ),
    "one value per data track"
  )

  expect_error(
    nice_circos(
      genome_build = "ce10",
      data_tracks = list(
        A = data.frame(
          chr = "chrI",
          start = 1,
          end = 10,
          value = 1,
          stringsAsFactors = FALSE
        ),
        B = data.frame(
          chr = "chrI",
          start = 20,
          end = 30,
          value = 2,
          stringsAsFactors = FALSE
        )
      ),
      track_types = c(A = "line")
    ),
    "one value per data track"
  )

  expect_error(
    nice_circos(
      genome_build = "ce10",
      data_tracks = list(),
      ideogram = FALSE
    ),
    "data_tracks"
  )
})

test_that("nice_circos errors when export_pdf directory does not exist", {
  skip_if_not_installed("circlize")

  bad_file <- file.path(tempdir(), "missing_directory_for_test", "plot.pdf")

  expect_error(
    nice_circos(
      data_tracks = toy_tracks_config(),
      custom_cytoband = toy_cytoband(),
      chromosome_index = c("chrA", "chrB"),
      plot_type = c("axis", "labels"),
      export_pdf = bad_file
    ),
    "does not exist"
  )
})


# Plotting with custom cytoband -------------------------------------------

test_that("nice_circos draws with custom cytoband without internet", {
  skip_if_not_installed("circlize")

  out_pdf <- tempfile(fileext = ".pdf")
  res <- NULL

  expect_error(
    res <- suppressWarnings(suppressMessages(
      nice_circos(
        data_tracks = toy_tracks_config(),
        custom_cytoband = toy_cytoband(),
        chromosome_index = c("chrA", "chrB"),
        plot_type = c("axis", "labels"),
        show_labels = TRUE,
        show_legend = FALSE,
        export_pdf = out_pdf
      )
    )),
    NA
  )

  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 0)

  expect_type(res, "list")
  expect_true(all(c(
    "tracks",
    "metadata",
    "track_settings",
    "link_data",
    "genome_build",
    "genome",
    "chromosomes",
    "plot_type",
    "cytoband",
    "pdf_opened"
  ) %in% names(res)))

  expect_equal(res$chromosomes, c("chrA", "chrB"))
  expect_equal(res$plot_type, c("axis", "labels"))
  expect_true(res$pdf_opened)

  expect_true(all(c("RNA", "CNV", "Mutation") %in% names(res$track_settings)))
  expect_equal(res$track_settings$RNA$type, "histogram")
  expect_equal(res$track_settings$CNV$type, "line")
  expect_equal(res$track_settings$Mutation$type, "mutation")
})

test_that("nice_circos filters tracks to selected chromosomes", {
  skip_if_not_installed("circlize")

  out_pdf <- tempfile(fileext = ".pdf")
  res <- NULL

  expect_error(
    res <- suppressWarnings(suppressMessages(
      nice_circos(
        data_tracks = toy_tracks_config(),
        custom_cytoband = toy_cytoband(),
        chromosome_index = "chrA",
        plot_type = c("axis", "labels"),
        export_pdf = out_pdf
      )
    )),
    NA
  )

  expect_equal(res$chromosomes, "chrA")
  expect_true(all(res$tracks$RNA$data$chr == "chrA"))
  expect_true(all(res$tracks$CNV$data$chr == "chrA"))
  expect_true(all(res$tracks$Mutation$data$chr == "chrA"))
})

test_that("nice_circos warns when a track becomes empty after chromosome filtering", {
  skip_if_not_installed("circlize")

  tracks <- list(
    OffTarget = data.frame(
      chr = "chrZ",
      start = 1,
      end = 10,
      value = 1,
      stringsAsFactors = FALSE
    )
  )

  out_pdf <- tempfile(fileext = ".pdf")

  expect_warning(
    suppressMessages(
      nice_circos(
        data_tracks = tracks,
        custom_cytoband = toy_cytoband(),
        chromosome_index = "chrA",
        plot_type = c("axis", "labels"),
        export_pdf = out_pdf
      )
    ),
    "no rows"
  )

  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 0)
})


# Plotting with built-in genome -------------------------------------------

test_that("nice_circos draws with built-in ce10 genome without internet", {
  skip_if_not_installed("circlize")

  cnv <- data.frame(
    chr = c("chrI", "chrII"),
    start = c(1, 1),
    end = c(1e5, 2e5),
    value = c(0.3, -0.2),
    name = c("seg1", "seg2"),
    stringsAsFactors = FALSE
  )

  out_pdf <- tempfile(fileext = ".pdf")
  res <- NULL

  expect_error(
    res <- suppressWarnings(suppressMessages(
      nice_circos(
        genome_build = "ce10",
        data_tracks = list(CNV = cnv),
        track_types = "line",
        chromosome_index = c("chrI", "chrII"),
        plot_type = c("axis", "labels"),
        show_legend = FALSE,
        export_pdf = out_pdf
      )
    )),
    NA
  )

  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 0)

  expect_equal(res$chromosomes, c("chrI", "chrII"))
  expect_true(is_circos_genome(res$genome))
  expect_equal(res$genome$source, "built-in")
})


# Metadata and links in plotting ------------------------------------------

test_that("nice_circos supports metadata rings", {
  skip_if_not_installed("circlize")

  out_pdf <- tempfile(fileext = ".pdf")

  metadata <- c(
    Sample = "TCGA-test",
    Subtype = "Basal"
  )

  res <- NULL

  expect_error(
    res <- suppressWarnings(suppressMessages(
      nice_circos(
        data_tracks = toy_tracks_config(),
        metadata = metadata,
        custom_cytoband = toy_cytoband(),
        chromosome_index = c("chrA", "chrB"),
        plot_type = c("axis", "labels"),
        export_pdf = out_pdf
      )
    )),
    NA
  )

  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 0)

  expect_equal(res$metadata$variable, c("Sample", "Subtype"))
  expect_equal(res$metadata$group, c("TCGA-test", "Basal"))
})

test_that("nice_circos supports genomic metadata rings", {
  skip_if_not_installed("circlize")

  metadata <- data.frame(
    chr = c("chrA", "chrB"),
    start = c(1, 1),
    end = c(500, 700),
    group = c("gain", "loss"),
    variable = c("CN_state", "CN_state"),
    stringsAsFactors = FALSE
  )

  out_pdf <- tempfile(fileext = ".pdf")
  res <- NULL

  expect_error(
    res <- suppressWarnings(suppressMessages(
      nice_circos(
        data_tracks = toy_tracks_config(),
        metadata = metadata,
        custom_cytoband = toy_cytoband(),
        chromosome_index = c("chrA", "chrB"),
        plot_type = c("axis", "labels"),
        export_pdf = out_pdf
      )
    )),
    NA
  )

  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 0)

  expect_equal(nrow(res$metadata), 2)
  expect_equal(unique(res$metadata$variable), "CN_state")
})

test_that("nice_circos supports genomic links", {
  skip_if_not_installed("circlize")

  links <- data.frame(
    chr1 = "chrA",
    start1 = 100,
    end1 = 200,
    chr2 = "chrB",
    start2 = 300,
    end2 = 400,
    stringsAsFactors = FALSE
  )

  out_pdf <- tempfile(fileext = ".pdf")
  res <- NULL

  expect_error(
    res <- suppressWarnings(suppressMessages(
      nice_circos(
        data_tracks = toy_tracks_config(),
        link_data = links,
        custom_cytoband = toy_cytoband(),
        chromosome_index = c("chrA", "chrB"),
        plot_type = c("axis", "labels"),
        export_pdf = out_pdf
      )
    )),
    NA
  )

  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 0)

  expect_equal(nrow(res$link_data), 1)
  expect_equal(res$link_data$chr1, "chrA")
  expect_equal(res$link_data$chr2, "chrB")
})


# Legend ------------------------------------------------------------------

test_that("nice_circos can draw a legend when ComplexHeatmap is installed", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("ComplexHeatmap")

  out_pdf <- tempfile(fileext = ".pdf")
  res <- NULL

  expect_error(
    res <- suppressWarnings(suppressMessages(
      nice_circos(
        data_tracks = toy_tracks_config(),
        custom_cytoband = toy_cytoband(),
        chromosome_index = c("chrA", "chrB"),
        plot_type = c("axis", "labels"),
        show_legend = TRUE,
        export_pdf = out_pdf
      )
    )),
    NA
  )

  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 0)
  expect_true(length(res$track_settings) > 0)
})


# File import --------------------------------------------------------------

test_that("nice_circos imports BED tracks without internet", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("rtracklayer")

  bed_file <- tempfile(fileext = ".bed")

  writeLines(
    c(
      "chrA\t0\t100\tmut1\t1",
      "chrB\t50\t150\tmut2\t1"
    ),
    con = bed_file
  )

  tracks <- list(
    Mutations = list(
      path = bed_file,
      type = "mutation",
      color = "#60A5FA",
      height = 0.06
    )
  )

  out_pdf <- tempfile(fileext = ".pdf")
  res <- NULL

  expect_error(
    res <- suppressWarnings(suppressMessages(
      nice_circos(
        data_tracks = tracks,
        custom_cytoband = toy_cytoband(),
        chromosome_index = c("chrA", "chrB"),
        plot_type = c("axis", "labels"),
        show_labels = TRUE,
        export_pdf = out_pdf
      )
    )),
    NA
  )

  expect_true(file.exists(out_pdf))
  expect_gt(file.info(out_pdf)$size, 0)

  expect_equal(res$track_settings$Mutations$type, "mutation")
  expect_equal(res$track_settings$Mutations$n, 2)
})

test_that("nice_circos errors for missing track files", {
  skip_if_not_installed("circlize")

  missing_file <- tempfile(fileext = ".bed")

  expect_error(
    nice_circos(
      data_tracks = list(
        Missing = list(
          path = missing_file,
          type = "mutation"
        )
      ),
      custom_cytoband = toy_cytoband(),
      chromosome_index = c("chrA", "chrB"),
      plot_type = c("axis", "labels"),
      export_pdf = tempfile(fileext = ".pdf")
    ),
    "Track file not found"
  )
})

test_that("nice_circos errors for unsupported track file extensions", {
  skip_if_not_installed("circlize")

  txt_file <- tempfile(fileext = ".txt")
  writeLines("chrA\t1\t10", txt_file)

  expect_error(
    nice_circos(
      data_tracks = list(
        BadFile = list(
          path = txt_file,
          type = "histogram"
        )
      ),
      custom_cytoband = toy_cytoband(),
      chromosome_index = c("chrA", "chrB"),
      plot_type = c("axis", "labels"),
      export_pdf = tempfile(fileext = ".pdf")
    ),
    "Unsupported file format"
  )
})
