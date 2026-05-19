test_that("trend_filter keeps genes with patient-level trend consistency", {
  fx <- make_trend_fixture()

  out <- trend_filter(
    expr = fx$expr,
    sampledata = fx$sampledata,
    results = list(Tumor_vs_Normal = fx$res),
    baseline = "normal",
    conditions = c(Tumor_vs_Normal = "tumor"),
    sample.col = "sample_id",
    patient.col = "patient_id",
    group.col = "sample_type",
    ratio = 1.1,
    scale = "linear"
  )

  expect_type(out, "list")
  expect_named(
    out,
    c("Tumor_vs_Normal", "TrendGenes", "Diagnostics", "Summary", "RemovedGenes")
  )

  expect_setequal(out$TrendGenes, c("g_up_pass", "g_down_pass"))
  expect_setequal(out$Tumor_vs_Normal$ensembl, c("g_up_pass", "g_down_pass"))

  expect_s3_class(out$Diagnostics, "data.frame")
  expect_s3_class(out$Summary, "data.frame")
})

test_that("trend_filter reports diagnostic reasons for removed genes", {
  fx <- make_trend_fixture()

  out <- trend_filter(
    expr = fx$expr,
    sampledata = fx$sampledata,
    results = list(Tumor_vs_Normal = fx$res),
    baseline = "normal",
    conditions = c(Tumor_vs_Normal = "tumor"),
    ratio = 1.1,
    scale = "linear"
  )

  diagnostics <- out$Diagnostics

  expect_equal(
    diagnostics$reason[diagnostics$gene == "g_up_fail"],
    "inconsistent_trend"
  )

  expect_equal(
    diagnostics$reason[diagnostics$gene == "g_down_fail"],
    "inconsistent_trend"
  )

  expect_equal(
    diagnostics$reason[diagnostics$gene == "g_no_direction"],
    "no_direction"
  )

  expect_equal(
    diagnostics$reason[diagnostics$gene == "g_missing_expr"],
    "missing_in_expr"
  )

  expect_true("g_up_fail" %in% out$RemovedGenes)
  expect_true("g_down_fail" %in% out$RemovedGenes)
})

test_that("trend_filter supports one data frame in results", {
  fx <- make_trend_fixture()

  out <- trend_filter(
    expr = fx$expr,
    sampledata = fx$sampledata,
    results = fx$res,
    baseline = "normal",
    conditions = "tumor",
    ratio = 1.1,
    scale = "linear"
  )

  expect_named(out, c("Comparison1", "TrendGenes", "Diagnostics", "Summary", "RemovedGenes"))
  expect_setequal(out$TrendGenes, c("g_up_pass", "g_down_pass"))
})

test_that("trend_filter handles log2-scale expression values", {
  fx <- make_trend_fixture(log_scale = TRUE)

  out <- trend_filter(
    expr = fx$expr,
    sampledata = fx$sampledata,
    results = list(Tumor_vs_Normal = fx$res),
    baseline = "normal",
    conditions = c(Tumor_vs_Normal = "tumor"),
    ratio = 1.1,
    scale = "log2"
  )

  expect_setequal(out$TrendGenes, c("g_up_pass", "g_down_pass"))
})

test_that("trend_filter warns or errors on unpaired patients depending on option", {
  fx <- make_trend_unpaired_fixture()

  expect_warning(
    trend_filter(
      expr = fx$expr,
      sampledata = fx$sampledata,
      results = list(Tumor_vs_Normal = fx$res),
      baseline = "normal",
      conditions = c(Tumor_vs_Normal = "tumor"),
      require_complete_pairs = FALSE
    ),
    "Using only paired patients",
    fixed = TRUE
  )

  expect_error(
    trend_filter(
      expr = fx$expr,
      sampledata = fx$sampledata,
      results = list(Tumor_vs_Normal = fx$res),
      baseline = "normal",
      conditions = c(Tumor_vs_Normal = "tumor"),
      require_complete_pairs = TRUE
    ),
    "Unpaired patients found",
    fixed = TRUE
  )
})

test_that("trend_filter validates required inputs", {
  fx <- make_trend_fixture()

  sampledata_missing_col <- fx$sampledata
  sampledata_missing_col$patient_id <- NULL

  expect_error(
    trend_filter(
      expr = fx$expr,
      sampledata = sampledata_missing_col,
      results = list(Tumor_vs_Normal = fx$res),
      baseline = "normal",
      conditions = c(Tumor_vs_Normal = "tumor")
    ),
    "columns are missing from sampledata",
    fixed = TRUE
  )

  expect_error(
    trend_filter(
      expr = fx$expr,
      sampledata = fx$sampledata,
      results = list(Tumor_vs_Normal = fx$res),
      baseline = "normal",
      conditions = c(Tumor_vs_Normal = "tumor"),
      ratio = 1
    ),
    "ratio must be a single numeric value greater than 1",
    fixed = TRUE
  )

  expect_error(
    trend_filter(
      expr = fx$expr,
      sampledata = fx$sampledata,
      results = list(Tumor_vs_Normal = fx$res),
      baseline = "normal",
      conditions = c(Tumor_vs_Normal = "tumor"),
      scale = "raw"
    )
  )
})

test_that("trend_filter supports multiple comparisons", {
  fx <- make_trend_fixture()

  res2 <- fx$res
  res2$log2FoldChange <- -res2$log2FoldChange

  out <- trend_filter(
    expr = fx$expr,
    sampledata = fx$sampledata,
    results = list(
      Tumor_vs_Normal = fx$res,
      Tumor_reversed = res2
    ),
    baseline = "normal",
    conditions = c(
      Tumor_vs_Normal = "tumor",
      Tumor_reversed = "tumor"
    ),
    ratio = 1.1,
    scale = "linear"
  )

  expect_true("Tumor_vs_Normal" %in% names(out))
  expect_true("Tumor_reversed" %in% names(out))
  expect_s3_class(out$Summary, "data.frame")
  expect_equal(nrow(out$Summary), 2)
})
