testthat::skip_if_not_installed("broom")
testthat::skip_if_not_installed("dplyr")
testthat::skip_if_not_installed("tibble")

test_that("get_glm returns tidy GLM results and model metadata", {
  dat <- make_binary_glm_data()

  res <- suppressMessages(get_glm(
    data = dat,
    outcome = "y",
    predictors = c("x1", "x2"),
    family = "binomial",
    adjust_method = "BH"
  ))

  expect_s3_class(res, "get_glm_result")
  expect_s3_class(res, "tbl_df")
  expect_equal(
    names(res),
    c(
      "term", "estimate", "ci_lower", "ci_upper", "std_error",
      "statistic", "p_value", "q_value", "significance"
    )
  )

  expect_false("(Intercept)" %in% res$term)
  expect_true(all(c("x1", "x2B") %in% res$term))
  expect_true(all(is.finite(res$estimate)))
  expect_true(all(res$estimate > 0))

  expect_s3_class(attr(res, "model"), "glm")
  expect_equal(attr(res, "formula"), "y ~ x1 + x2")
  expect_equal(attr(res, "family"), "binomial")
  expect_equal(attr(res, "link"), "logit")
  expect_equal(attr(res, "n_obs"), nrow(dat))
  expect_true(is.finite(attr(res, "AIC")))
  expect_true(attr(res, "exponentiate"))
  expect_equal(attr(res, "adjust_method"), "BH")

  expect_equal(
    res$q_value,
    stats::p.adjust(res$p_value, method = "BH"),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})

test_that("get_glm can keep the intercept and adjust all reported p-values", {
  dat <- make_binary_glm_data()

  res <- suppressMessages(get_glm(
    data = dat,
    outcome = "y",
    predictors = c("x1", "x2"),
    family = stats::binomial(),
    adjust_method = "holm",
    remove_intercept = FALSE
  ))

  expect_true("(Intercept)" %in% res$term)
  expect_equal(attr(res, "adjust_method"), "holm")
  expect_equal(
    res$q_value,
    stats::p.adjust(res$p_value, method = "holm"),
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
})

test_that("get_glm auto exponentiation follows the link and can be overridden", {
  bin_dat <- make_binary_glm_data()
  gauss_dat <- make_gaussian_glm_data()
  pois_dat <- make_count_glm_data()

  res_gauss <- suppressMessages(get_glm(
    data = gauss_dat,
    outcome = "y",
    predictors = c("x1", "x2"),
    family = stats::gaussian()
  ))
  expect_false(attr(res_gauss, "exponentiate"))
  expect_equal(attr(res_gauss, "link"), "identity")

  res_pois <- suppressMessages(get_glm(
    data = pois_dat,
    outcome = "y",
    predictors = "x1",
    family = stats::poisson
  ))
  expect_true(attr(res_pois, "exponentiate"))
  expect_equal(attr(res_pois, "link"), "log")
  expect_true(all(res_pois$estimate > 0))

  res_no_exp <- suppressMessages(get_glm(
    data = bin_dat,
    outcome = "y",
    predictors = c("x1", "x2"),
    family = "binomial",
    exponentiate = FALSE
  ))
  expect_false(attr(res_no_exp, "exponentiate"))
})

test_that("get_glm safely quotes non-syntactic and reserved variable names", {
  dat <- make_binary_glm_data()
  weird_dat <- data.frame(
    "case status" = dat$y,
    "age years" = dat$x1,
    "if" = as.numeric(dat$x2 == "B"),
    check.names = FALSE
  )

  res <- suppressMessages(get_glm(
    data = weird_dat,
    outcome = "case status",
    predictors = c("age years", "if"),
    family = "binomial"
  ))

  expect_equal(attr(res, "formula"), "`case status` ~ `age years` + `if`")
  expect_equal(nrow(res), 2L)
  expect_s3_class(attr(res, "model"), "glm")
})

test_that("get_glm validates inputs before fitting", {
  dat <- make_binary_glm_data()

  expect_error(
    get_glm(as.list(dat), "y", c("x1", "x2")),
    "`data` must be"
  )
  expect_error(
    get_glm(dat, c("y", "other"), c("x1", "x2")),
    "`outcome` must be"
  )
  expect_error(
    get_glm(dat, "missing_y", c("x1", "x2")),
    "not found"
  )
  expect_error(
    get_glm(dat, "y", character()),
    "`predictors` must be"
  )
  expect_error(
    get_glm(dat, "y", c("x1", "missing_x")),
    "not in `data`"
  )
  expect_error(
    get_glm(dat, "y", c("x1", "x2"), conf_level = 1),
    "`conf_level` must be"
  )
  expect_error(
    get_glm(dat, "y", c("x1", "x2"), exponentiate = NA),
    "`exponentiate` must be"
  )
  expect_error(
    get_glm(dat, "y", c("x1", "x2"), remove_intercept = NA),
    "`remove_intercept` must be"
  )
  expect_error(
    get_glm(dat, "y", c("x1", "x2"), verbose = NA),
    "`verbose` must be"
  )
  expect_error(
    get_glm(dat, "y", c("x1", "x2"), adjust_method = "bad_method"),
    "should be one of"
  )
  expect_error(
    get_glm(dat, "y", c("x1", "x2"), family = "not_a_family"),
    "Could not find"
  )
  expect_error(
    get_glm(
      dat, "y", c("x1", "x2"),
      family = function(link) stats::binomial(link = link)
    ),
    "could not be evaluated"
  )
})

test_that("print.get_glm_result prints a compact summary and returns input invisibly", {
  dat <- make_binary_glm_data()
  res <- suppressMessages(get_glm(dat, "y", c("x1", "x2"), family = "binomial"))

  printed <- capture.output(returned <- print(res))

  expect_identical(returned, res)
  expect_true(any(grepl("get_glm result", printed, fixed = TRUE)))
  expect_true(any(grepl("q-value correction: BH", printed, fixed = TRUE)))

  bad <- res
  bad$p_value <- NULL
  expect_error(print(bad), "missing columns")
})
