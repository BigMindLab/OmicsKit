# Helper data for get_glm() and nice_ROC() tests.
# testthat automatically sources files named helper*.R before tests.

make_binary_glm_data <- function(n = 160, seed = 123) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- factor(sample(c("A", "B"), n, replace = TRUE), levels = c("A", "B"))
  eta <- -0.25 + 0.75 * x1 + 0.45 * (x2 == "B")
  p <- stats::plogis(eta)
  y <- stats::rbinom(n, size = 1, prob = p)

  data.frame(
    y = as.integer(y),
    x1 = x1,
    x2 = x2
  )
}

make_gaussian_glm_data <- function(n = 120, seed = 456) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rnorm(n)
  y <- 1.5 + 0.8 * x1 - 0.35 * x2 + stats::rnorm(n, sd = 0.7)

  data.frame(y = y, x1 = x1, x2 = x2)
}

make_count_glm_data <- function(n = 140, seed = 789) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  exposure <- stats::runif(n, min = 0.8, max = 1.5)
  lambda <- exp(0.2 + 0.35 * x1 + log(exposure))
  y <- stats::rpois(n, lambda = lambda)

  data.frame(y = y, x1 = x1, exposure = exposure)
}
