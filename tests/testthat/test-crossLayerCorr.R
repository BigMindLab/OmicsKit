test_that("crossLayerCorr matches shared features and samples", {
  set.seed(1)
  mat_x <- matrix(rnorm(30), nrow = 5)
  mat_y <- mat_x + matrix(rnorm(30, sd = 0.1), nrow = 5)

  rownames(mat_x) <- rownames(mat_y) <- paste0("G", 1:5)
  colnames(mat_x) <- colnames(mat_y) <- paste0("S", 1:6)

  res <- crossLayerCorr(mat_x, mat_y, top_n = 4, plot = FALSE)

  expect_equal(res$n_shared_samples, 6)
  expect_equal(res$n_shared_features, 5)
  expect_equal(res$n_features_used, 4)
  expect_true(all(c("sample", "correlation") %in% names(res$correlations)))
})
