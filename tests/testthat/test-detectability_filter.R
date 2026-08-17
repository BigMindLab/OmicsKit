test_that("detectability_filter keeps detectable UP and DOWN genes in one comparison", {
  fx <- make_detect_fixture()

  out <- detectability_filter(
    norm.counts = fx$norm_counts,
    df.BvsA = fx$df_b,
    samples.baseline = fx$samples_baseline,
    samples.condition1 = fx$samples_condition1,
    cutoffs = c(50, 50, 0)
  )

  expect_type(out, "list")
  expect_named(out, c("Comparison1", "DetectGenes"))

  expect_setequal(
    out$DetectGenes,
    c("g_up_detect", "g_down_detect")
  )

  expect_setequal(
    out$Comparison1$ensembl,
    c("g_up_detect", "g_down_detect")
  )

  expect_false("g_low_baseMean" %in% out$DetectGenes)
  expect_false("g_up_low_expr" %in% out$DetectGenes)
  expect_false("g_down_low_baseline" %in% out$DetectGenes)
  expect_false("g_no_direction" %in% out$DetectGenes)
})

test_that("detectability_filter supports optional second and third comparisons", {
  fx <- make_detect_fixture()

  out <- detectability_filter(
    norm.counts = fx$norm_counts,
    df.BvsA = fx$df_b,
    df.CvsA = fx$df_c,
    df.DvsA = fx$df_d,
    samples.baseline = fx$samples_baseline,
    samples.condition1 = fx$samples_condition1,
    samples.condition2 = fx$samples_condition2,
    samples.condition3 = fx$samples_condition3,
    cutoffs = c(50, 50, 0)
  )

  expect_named(
    out,
    c("Comparison1", "DetectGenes", "Comparison2", "Comparison3")
  )

  expect_true("g_c_detect" %in% out$DetectGenes)
  expect_true("g_c_detect" %in% out$Comparison2$ensembl)

  expect_setequal(
    out$Comparison1$ensembl,
    c("g_up_detect", "g_down_detect")
  )
})

test_that("detectability_filter removes duplicated genes from DetectGenes", {
  fx <- make_detect_fixture()

  df_b_dup <- rbind(fx$df_b, fx$df_b[fx$df_b$ensembl == "g_up_detect", ])

  out <- detectability_filter(
    norm.counts = fx$norm_counts,
    df.BvsA = df_b_dup,
    samples.baseline = fx$samples_baseline,
    samples.condition1 = fx$samples_condition1,
    cutoffs = c(50, 50, 0)
  )

  expect_equal(length(out$DetectGenes), length(unique(out$DetectGenes)))
})

test_that("detectability_filter validates cutoff length", {
  fx <- make_detect_fixture()

  expect_error(
    detectability_filter(
      norm.counts = fx$norm_counts,
      df.BvsA = fx$df_b,
      samples.baseline = fx$samples_baseline,
      samples.condition1 = fx$samples_condition1,
      cutoffs = c(50, 50)
    ),
    "Cutoffs vector must contain three values",
    fixed = TRUE
  )
})
