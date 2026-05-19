test_that("nice_KM returns a ggplot object for two strata", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")
  skip_if_not_installed("survminer")

  dat <- make_km_test_data()

  p <- nice_KM(
    data = dat,
    gene = "GENE_muts",
    time_var = "PFI.time",
    event_var = "PFI",
    conf_int = FALSE,
    pval_size = 3
  )

  expect_s3_class(p, "ggplot")
})

test_that("nice_KM returnData returns the survfit object and plot", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")
  skip_if_not_installed("survminer")

  dat <- make_km_test_data()

  out <- nice_KM(
    data = dat,
    gene = "GENE_muts",
    time_var = "PFI.time",
    event_var = "PFI",
    conf_int = FALSE,
    returnData = TRUE,
    pval_size = 3
  )

  expect_true(identical(names(out), c("km_fit", "plot")))
  expect_s3_class(out$km_fit, "survfit")
  expect_s3_class(out$plot, "ggplot")
})

test_that("nice_KM supports three strata when a matching palette is supplied", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")
  skip_if_not_installed("survminer")

  dat <- make_km_test_data()

  p <- nice_KM(
    data = dat,
    gene = "PAM50_group",
    time_var = "PFI.time",
    event_var = "PFI",
    title_prefix = "Subtype ",
    colors = c("#1F8FFF", "#ED4D4D", "#2ECC71"),
    conf_int = FALSE,
    pval_size = 3
  )

  expect_s3_class(p, "ggplot")
})

test_that("nice_KM handles variables with only one observed category", {
  skip_if_not_installed("ggplot2")

  dat <- make_km_test_data()
  dat$GENE_muts <- factor("No", levels = c("No", "Yes"))

  p <- NULL
  expect_warning(
    p <- nice_KM(
      data = dat,
      gene = "GENE_muts",
      time_var = "PFI.time",
      event_var = "PFI"
    ),
    "Only one category"
  )

  expect_s3_class(p, "ggplot")
})

test_that("nice_KM singleton-category branch also works with returnData", {
  skip_if_not_installed("ggplot2")

  dat <- make_km_test_data()
  dat$GENE_muts <- "No"

  out <- NULL
  expect_warning(
    out <- nice_KM(
      data = dat,
      gene = "GENE_muts",
      time_var = "PFI.time",
      event_var = "PFI",
      returnData = TRUE
    ),
    "Only one category"
  )

  expect_true(identical(names(out), c("km_fit", "plot")))
  expect_null(out$km_fit)
  expect_s3_class(out$plot, "ggplot")
})
