testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("ggplot2")
testthat::skip_if_not_installed("pROC")
testthat::skip_if_not_installed("tibble")

test_that("nice_ROC returns plot, AUC table, DeLong test, and ROC objects", {
  dat <- make_binary_glm_data(n = 180, seed = 321)
  fit_x1 <- stats::glm(y ~ x1, data = dat, family = stats::binomial())
  fit_x1_x2 <- stats::glm(y ~ x1 + x2, data = dat, family = stats::binomial())

  out <- nice_ROC(
    models = list("x1 only" = fit_x1, "x1 + x2" = fit_x1_x2),
    data = dat,
    outcome = "y",
    return_data = TRUE
  )

  expect_type(out, "list")
  expect_named(out, c("plot", "auc_table", "delong_test", "roc_objects"))
  expect_s3_class(out$plot, "ggplot")

  expect_s3_class(out$auc_table, "tbl_df")
  expect_equal(nrow(out$auc_table), 2L)
  expect_named(
    out$auc_table,
    c("model", "auc", "ci_lower", "ci_upper", "n_cases", "n_controls")
  )
  expect_equal(out$auc_table$model, c("x1 only", "x1 + x2"))
  expect_true(all(out$auc_table$auc >= 0 & out$auc_table$auc <= 1))
  expect_true(all(out$auc_table$ci_lower >= 0 & out$auc_table$ci_lower <= 1))
  expect_true(all(out$auc_table$ci_upper >= 0 & out$auc_table$ci_upper <= 1))
  expect_equal(unique(out$auc_table$n_cases), sum(dat$y == 1))
  expect_equal(unique(out$auc_table$n_controls), sum(dat$y == 0))

  expect_s3_class(out$delong_test, "htest")
  expect_named(out$roc_objects, c("x1 only", "x1 + x2"))
  expect_true(all(vapply(out$roc_objects, inherits, logical(1), "roc")))
})

test_that("nice_ROC returns a ggplot object when return_data is FALSE", {
  dat <- make_binary_glm_data(n = 150, seed = 654)
  fit <- stats::glm(y ~ x1 + x2, data = dat, family = stats::binomial())

  p <- nice_ROC(
    models = list("full model" = fit),
    data = dat,
    outcome = "y",
    plot_title = "Test ROC",
    plot_subtitle = "Held-out sample"
  )

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Test ROC")
  expect_equal(p$labels$subtitle, "Held-out sample")
})

test_that("nice_ROC accepts numeric probability vectors", {
  dat <- make_binary_glm_data(n = 150, seed = 987)
  fit <- stats::glm(y ~ x1 + x2, data = dat, family = stats::binomial())
  probs <- stats::predict(fit, newdata = dat, type = "response")

  out <- nice_ROC(
    models = list("probabilities" = as.numeric(probs)),
    data = dat,
    outcome = "y",
    show_ci = FALSE,
    return_data = TRUE
  )

  expect_equal(out$auc_table$model, "probabilities")
  expect_null(out$delong_test)
  expect_named(out$roc_objects, "probabilities")
  expect_s3_class(out$plot, "ggplot")
})

test_that("nice_ROC can compare two models without adding DeLong output", {
  dat <- make_binary_glm_data(n = 160, seed = 246)
  fit_x1 <- stats::glm(y ~ x1, data = dat, family = stats::binomial())
  fit_x1_x2 <- stats::glm(y ~ x1 + x2, data = dat, family = stats::binomial())

  out <- nice_ROC(
    models = list("x1 only" = fit_x1, "x1 + x2" = fit_x1_x2),
    data = dat,
    outcome = "y",
    show_delong = FALSE,
    return_data = TRUE
  )

  expect_null(out$delong_test)
  expect_equal(nrow(out$auc_table), 2L)
})

test_that("nice_ROC validates inputs", {
  dat <- make_binary_glm_data(n = 120, seed = 135)
  fit <- stats::glm(y ~ x1 + x2, data = dat, family = stats::binomial())

  expect_error(
    nice_ROC(models = list(fit), data = dat, outcome = "y"),
    "named"
  )
  expect_error(
    nice_ROC(models = list("model" = fit), data = dat, outcome = c("y", "z")),
    "single character string"
  )
  expect_error(
    nice_ROC(models = list("model" = fit), data = dat, outcome = "missing_y"),
    "not found"
  )

  non_binary <- dat
  non_binary$y[1] <- 2L
  expect_error(
    nice_ROC(models = list("model" = fit), data = non_binary, outcome = "y"),
    "binary 0/1"
  )
  expect_error(
    nice_ROC(models = list("prob" = rep(0.5, nrow(dat))), data = NULL, outcome = "y"),
    "data` must be provided"
  )
  expect_error(
    nice_ROC(models = list("prob" = rep(0.5, 3)), data = dat, outcome = "y"),
    "has length"
  )
  expect_error(
    nice_ROC(models = list("bad" = "not a model"), data = dat, outcome = "y"),
    "glm object or numeric probability vector"
  )
})
