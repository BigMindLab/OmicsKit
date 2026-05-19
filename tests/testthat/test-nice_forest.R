# Tests for nice_forest().
# Place this file in tests/testthat/.

test_that("nice_forest returns a ggplot object", {
  skip_if_not_installed("ggplot2")

  tab <- make_forest_test_data()

  p <- nice_forest(tab, title = "PFI adjusted Cox forest plot")

  expect_s3_class(p, "ggplot")
})

test_that("nice_forest can return the filtered plot table", {
  skip_if_not_installed("ggplot2")

  tab <- make_forest_test_data()

  out <- nice_forest(
    tab,
    p_display = 0.02,
    sort_by = "input",
    return_table = TRUE
  )

  expect_named(out, c("plot", "table"))
  expect_s3_class(out$plot, "ggplot")
  expect_s3_class(out$table, "data.frame")
  expect_equal(nrow(out$table), 2)
  expect_true(all(out$table$.p_plot <= 0.02))
  expect_true(all(out$table$.significance_plot == "p < 0.05"))
})

test_that("nice_forest sorts rows by estimate or p-value", {
  skip_if_not_installed("ggplot2")

  tab <- make_forest_test_data()

  by_estimate <- nice_forest(tab, sort_by = "estimate", return_table = TRUE)$table
  expect_equal(by_estimate$.estimate_plot, sort(tab$HR))

  by_p <- nice_forest(tab, sort_by = "p.value", return_table = TRUE)$table
  expect_equal(by_p$.p_plot, sort(tab$p.value))
})

test_that("nice_forest supports user-provided label columns", {
  skip_if_not_installed("ggplot2")

  tab <- make_forest_test_data()
  tab$pretty_label <- c("PAM50 HER2", "DNA methylation cluster 3", "RPPA Reactive")

  out <- nice_forest(
    tab,
    label_col = "pretty_label",
    sort_by = "input",
    return_table = TRUE
  )

  expect_equal(as.character(out$table$.label_plot), tab$pretty_label)
})

test_that("nice_forest validates input columns and scalar arguments", {
  skip_if_not_installed("ggplot2")

  tab <- make_forest_test_data()

  expect_error(
    nice_forest(tab[, setdiff(names(tab), "CI_high")]),
    "required columns were not found"
  )

  expect_error(
    nice_forest(tab, p_display = 2),
    "'p_display' should be a single numeric value between 0 and 1",
    fixed = TRUE
  )

  expect_error(
    nice_forest(tab, vline = c(1, 2)),
    "'vline' should be a single numeric value",
    fixed = TRUE
  )
})

test_that("nice_forest errors when filtering removes all rows", {
  skip_if_not_installed("ggplot2")

  tab <- make_forest_test_data()

  expect_error(
    nice_forest(tab, p_display = 0.001),
    "No rows remained after filtering"
  )
})

test_that("nice_forest requires positive estimates on a log scale", {
  skip_if_not_installed("ggplot2")

  tab <- make_forest_test_data()
  tab$HR[1] <- 0

  expect_error(
    nice_forest(tab, log_scale = TRUE),
    "must be > 0"
  )

  expect_s3_class(nice_forest(tab, log_scale = FALSE), "ggplot")
})
