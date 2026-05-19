# Tests for get_cox().
# Place this file in tests/testthat/.

test_that("get_cox validates required inputs", {
  skip_if_not_installed("survival")
  skip_if_not_installed("broom")
  skip_if_not_installed("dplyr")

  dat <- make_cox_test_data()

  expect_error(
    get_cox(list(), vars = "gene_status", verbose = FALSE),
    "'data' should be a data.frame",
    fixed = TRUE
  )

  expect_error(
    get_cox(dat, time_col = "missing_time", vars = "gene_status", verbose = FALSE),
    "Column 'missing_time' not found"
  )

  expect_error(
    get_cox(dat, event_col = "missing_event", vars = "gene_status", verbose = FALSE),
    "Column 'missing_event' not found"
  )

  expect_error(
    get_cox(dat, vars = "does_not_exist", verbose = FALSE),
    "variables were not found"
  )

  expect_error(
    get_cox(dat, vars = "gene_status", model = "adjusted", verbose = FALSE),
    "'adjust_vars' is required",
    fixed = TRUE
  )
})

test_that("get_cox fits univariable categorical Cox models", {
  skip_if_not_installed("survival")
  skip_if_not_installed("broom")
  skip_if_not_installed("dplyr")

  dat <- make_cox_test_data()

  out <- get_cox(
    data = dat,
    time_col = "PFI.time",
    event_col = "PFI",
    vars = "gene_status",
    model = "univariable",
    ref_levels = list(gene_status = "WT"),
    min_events_per_level = 1,
    verbose = FALSE
  )

  expect_s3_class(out, "data.frame")
  expect_true(all(c(
    "model", "model_id", "variable", "term", "term_clean", "reference",
    "HR", "CI_low", "CI_high", "p.value", "n_used", "n_events"
  ) %in% names(out)))
  expect_equal(unique(out$model), "univariable")
  expect_equal(unique(out$variable), "gene_status")
  expect_equal(unique(out$reference), "WT")
  expect_true(all(is.finite(out$HR)))
  expect_true(all(out$CI_low > 0 & out$CI_high > 0))
  expect_true(all(out$p.value >= 0 & out$p.value <= 1))
})

test_that("get_cox fits multivariable models with categorical and continuous terms", {
  skip_if_not_installed("survival")
  skip_if_not_installed("broom")
  skip_if_not_installed("dplyr")

  dat <- make_cox_test_data()

  out <- get_cox(
    data = dat,
    time_col = "PFI.time",
    event_col = "PFI",
    vars = c("gene_status", "age"),
    model = "multivariable",
    ref_levels = list(gene_status = "WT"),
    min_events_per_level = 1,
    verbose = FALSE
  )

  expect_equal(unique(out$model), "multivariable")
  expect_equal(unique(out$model_id), "multivariable")
  expect_setequal(unique(out$variable), c("gene_status", "age"))
  expect_true(any(out$reference == "continuous"))
  expect_true(all(is.finite(out$HR)))
})

test_that("get_cox adjusted models drop adjustment terms unless requested", {
  skip_if_not_installed("survival")
  skip_if_not_installed("broom")
  skip_if_not_installed("dplyr")

  dat <- make_cox_test_data()

  out_drop <- get_cox(
    data = dat,
    time_col = "PFI.time",
    event_col = "PFI",
    vars = "gene_status",
    model = "adjusted",
    adjust_vars = "age",
    keep_adjust_terms = FALSE,
    ref_levels = list(gene_status = "WT"),
    min_events_per_level = 1,
    verbose = FALSE
  )

  expect_equal(unique(out_drop$model), "adjusted")
  expect_equal(unique(out_drop$variable), "gene_status")
  expect_equal(unique(out_drop$adjusted_for), "age")

  out_keep <- get_cox(
    data = dat,
    time_col = "PFI.time",
    event_col = "PFI",
    vars = "gene_status",
    model = "adjusted",
    adjust_vars = "age",
    keep_adjust_terms = TRUE,
    ref_levels = list(gene_status = "WT"),
    min_events_per_level = 1,
    verbose = FALSE
  )

  expect_setequal(unique(out_keep$variable), c("gene_status", "age"))
})

test_that("get_cox skips sparse categorical variables but keeps valid models", {
  skip_if_not_installed("survival")
  skip_if_not_installed("broom")
  skip_if_not_installed("dplyr")

  dat <- make_cox_test_data()

  out <- get_cox(
    data = dat,
    time_col = "PFI.time",
    event_col = "PFI",
    vars = c("gene_status", "sparse_var"),
    model = "univariable",
    ref_levels = list(gene_status = "WT"),
    min_n_per_level = 10,
    min_events_per_level = 1,
    verbose = FALSE
  )

  expect_equal(unique(out$variable), "gene_status")
  expect_false("sparse_var" %in% out$variable)
})

test_that("get_cox automatic variable selection excludes ID-like, survival-like, and underscored columns", {
  skip_if_not_installed("survival")
  skip_if_not_installed("broom")
  skip_if_not_installed("dplyr")

  dat <- make_cox_test_data()[, c(
    "sample", "PFI.time", "PFI", "gene_status", "days_to_death", "_hidden_omic"
  )]

  out <- get_cox(
    data = dat,
    time_col = "PFI.time",
    event_col = "PFI",
    vars = NULL,
    model = "univariable",
    ref_levels = list(gene_status = "WT"),
    min_events_per_level = 1,
    verbose = FALSE
  )

  expect_equal(unique(out$variable), "gene_status")
})

test_that("get_cox errors when no Cox model can be fitted", {
  skip_if_not_installed("survival")
  skip_if_not_installed("broom")
  skip_if_not_installed("dplyr")

  dat <- make_cox_test_data()
  dat$PFI <- 2L

  expect_error(
    get_cox(
      data = dat,
      time_col = "PFI.time",
      event_col = "PFI",
      vars = "gene_status",
      model = "univariable",
      min_events_per_level = 1,
      verbose = FALSE
    ),
    "No Cox model could be fitted"
  )
})
