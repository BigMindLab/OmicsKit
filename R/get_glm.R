#####################
# Function get_glm  #
#####################

# Helper: %||% (NULL coalescing) - defined locally to avoid external dependency
`%||%` <- function(a, b) if (!is.null(a)) a else b


# Helper: quote variable names safely for formula construction
.quote_name <- function(x) {
  needs_quotes <- !grepl("^[A-Za-z.][A-Za-z0-9._]*$", x) ||
    x %in% c(
      "if", "else", "repeat", "while", "function", "for", "in", "next",
      "break", "TRUE", "FALSE", "NULL", "Inf", "NaN", "NA", "NA_integer_",
      "NA_real_", "NA_complex_", "NA_character_"
    )

  if (needs_quotes) {
    paste0("`", gsub("`", "\\\\`", x), "`")
  } else {
    x
  }
}


# Helper: normalize family input into a proper family object
.normalize_glm_family <- function(family, envir = parent.frame()) {

  if (inherits(family, "family")) {
    return(family)
  }

  if (is.character(family)) {
    if (length(family) != 1 || is.na(family)) {
      stop("`family` must be a single character string, family function, or family object.",
           call. = FALSE)
    }

    fam_fun <- NULL

    # First look in the calling environment / attached packages
    if (exists(family, mode = "function", envir = envir, inherits = TRUE)) {
      fam_fun <- get(family, mode = "function", envir = envir, inherits = TRUE)
    }

    # Then look explicitly in stats
    if (is.null(fam_fun) &&
        exists(family, mode = "function", envir = asNamespace("stats"),
               inherits = FALSE)) {
      fam_fun <- get(family, mode = "function", envir = asNamespace("stats"),
                     inherits = FALSE)
    }

    if (is.null(fam_fun)) {
      stop(sprintf("Could not find a GLM family function called '%s'.", family),
           call. = FALSE)
    }

    family <- fam_fun
  }

  if (is.function(family)) {
    family <- tryCatch(
      family(),
      error = function(e) {
        stop(sprintf(
          "`family` function could not be evaluated without arguments: %s",
          conditionMessage(e)
        ), call. = FALSE)
      }
    )
  }

  if (!inherits(family, "family")) {
    stop("`family` must be a character string, family function, or valid family object.",
         call. = FALSE)
  }

  family
}


#' Function to Fit a GLM and Return Tidy Results with Effect Sizes and FDR Correction
#'
#' @description
#' Fits a Generalized Linear Model (GLM) for binary, continuous, count,
#' proportion, rate, or positive skewed outcomes, depending on the specified
#' \code{family}. The function returns a tidy \code{tibble} with coefficients,
#' optional exponentiated effect sizes, 95% confidence intervals, Wald test
#' statistics, raw p-values, and multiple-testing-adjusted q-values.
#'
#' The fitted \code{glm} object is attached as an attribute so that it can be
#' passed directly to downstream functions such as \code{\link{nice_ROC}}
#' without re-fitting.
#'
#' @param data A \code{data.frame} or \code{tibble} containing all outcome and
#'   predictor variables.
#'
#' @param outcome Character string of length 1. Name of the outcome column.
#'   For logistic regression this is typically coded as \code{0}/\code{1},
#'   although other binomial formats accepted by \code{\link[stats]{glm}} can
#'   also be used.
#'
#' @param predictors Character vector. Names of the predictor columns to include.
#'   Categorical variables should be \code{factor} or \code{character}; they are
#'   handled by \code{\link[stats]{glm}} contrast coding automatically.
#'
#' @param family Character string, family function, or \code{\link[stats]{family}}
#'   object passed to \code{\link[stats]{glm}}. It specifies both the assumed
#'   distribution of the outcome and the link function used in the GLM.
#'   Common options include:
#'   \itemize{
#'     \item \code{"binomial"} or \code{binomial(link = "logit")} for binary,
#'       proportion, or success/failure outcomes; this corresponds to logistic
#'       regression when the link is \code{"logit"}.
#'     \item \code{"gaussian"} or \code{gaussian(link = "identity")} for
#'       approximately normally distributed continuous outcomes.
#'     \item \code{"poisson"} or \code{poisson(link = "log")} for count or rate
#'       outcomes.
#'     \item \code{"quasibinomial"} and \code{"quasipoisson"} for binomial or
#'       Poisson-type outcomes with overdispersion.
#'     \item \code{Gamma(link = "log")} for positive, right-skewed continuous
#'       outcomes.
#'     \item \code{inverse.gaussian()} for positive continuous outcomes with
#'       variance increasing strongly with the mean.
#'   }
#'   Custom family objects can also be supplied. Default: \code{"binomial"}.
#'
#' @param adjust_method Character string. Method for p-value adjustment passed
#'   to \code{\link[stats]{p.adjust}}. Options: \code{"BH"},
#'   \code{"bonferroni"}, \code{"holm"}, \code{"BY"}, \code{"fdr"},
#'   \code{"none"}. Default: \code{"BH"}.
#'
#' @param conf_level Numeric value in \code{(0, 1)}. Confidence level for the
#'   interval. Default: \code{0.95}.
#'
#' @param exponentiate Logical or \code{NULL}. If \code{TRUE}, coefficients and
#'   confidence intervals are exponentiated. If \code{NULL}, exponentiation is
#'   applied automatically when the model link is \code{"logit"} or \code{"log"}.
#'   This gives odds ratios for logistic models with logit link, rate ratios for
#'   Poisson-type models with log link, and multiplicative effects for other
#'   log-link models. Default: \code{NULL}.
#'
#' @param remove_intercept Logical. Whether to drop the \code{(Intercept)} row
#'   from the returned table. Default: \code{TRUE}.
#'
#' @param verbose Logical. If \code{TRUE}, prints the model summary to console.
#'   Default: \code{FALSE}.
#'
#' @return
#' A \code{tibble} of class
#' \code{c("get_glm_result", "tbl_df", "tbl", "data.frame")} with the following
#' columns:
#' \describe{
#'   \item{\code{term}}{Predictor name; for factors, includes the level
#'         according to the contrast coding used by \code{glm}.}
#'   \item{\code{estimate}}{Exponentiated coefficient if
#'         \code{exponentiate = TRUE}; otherwise the raw coefficient on the
#'         model linear predictor scale.}
#'   \item{\code{ci_lower}}{Lower bound of the confidence interval on the same
#'         scale as \code{estimate}.}
#'   \item{\code{ci_upper}}{Upper bound of the confidence interval on the same
#'         scale as \code{estimate}.}
#'   \item{\code{std_error}}{Standard error of the coefficient on the linear
#'         predictor scale.}
#'   \item{\code{statistic}}{Wald test statistic. For fixed-dispersion families
#'         such as binomial and Poisson this is typically a z-statistic; for
#'         families with estimated dispersion it may be a t-statistic.}
#'   \item{\code{p_value}}{Two-sided p-value from the Wald test.}
#'   \item{\code{q_value}}{P-value adjusted for multiple comparisons using
#'         the method specified in \code{adjust_method}.}
#'   \item{\code{significance}}{Star annotation based on raw p-value:
#'         \code{"***"} p<0.001, \code{"**"} p<0.01, \code{"*"} p<0.05,
#'         \code{"."} p<0.10, \code{" "} otherwise.}
#' }
#'
#' The following attributes are attached to the returned object:
#' \describe{
#'   \item{\code{model}}{The fitted \code{glm} object.}
#'   \item{\code{formula}}{Character string of the model formula.}
#'   \item{\code{family}}{Family name used.}
#'   \item{\code{link}}{Link function used.}
#'   \item{\code{n_obs}}{Number of observations used in model fitting.}
#'   \item{\code{AIC}}{Akaike Information Criterion. May be \code{NA} for
#'         quasi-likelihood families.}
#'   \item{\code{exponentiate}}{Logical indicating whether coefficients were
#'         exponentiated.}
#'   \item{\code{adjust_method}}{P-value adjustment method used.}
#' }
#'
#' @details
#' Q-values are computed across all reported terms after removing the intercept
#' if \code{remove_intercept = TRUE}. Therefore, the multiple-testing correction
#' pool matches the number of hypotheses shown in the returned table.
#'
#' Automatic exponentiation is based on the link function. Models with
#' \code{logit} link are reported as odds ratios; models with \code{log} link
#' are reported as multiplicative effects such as rate ratios or mean ratios.
#' Models with identity, probit, cloglog, inverse, or other links are not
#' exponentiated automatically.
#'
#' @seealso \code{\link{nice_ROC}} for ROC curve visualisation of logistic
#'   models produced by \code{get_glm}.
#'
#' @examples
#' \dontrun{
#' # Binary outcome: logistic regression
#' res_logit <- get_glm(
#'   data       = train_data,
#'   outcome    = "stage_advanced",
#'   predictors = c("age", "ER", "PR", "HER2", "histology", "menopause"),
#'   family     = "binomial"
#' )
#'
#' # Continuous outcome: linear model through glm
#' res_gaussian <- get_glm(
#'   data       = train_data,
#'   outcome    = "tumor_size",
#'   predictors = c("age", "ER", "PR", "HER2"),
#'   family     = "gaussian"
#' )
#'
#' # Count outcome: Poisson regression
#' res_pois <- get_glm(
#'   data       = train_data,
#'   outcome    = "n_mutations",
#'   predictors = c("age", "stage", "histology"),
#'   family     = "poisson"
#' )
#'
#' # Positive skewed continuous outcome
#' res_gamma <- get_glm(
#'   data       = train_data,
#'   outcome    = "cost",
#'   predictors = c("age", "stage", "treatment"),
#'   family     = Gamma(link = "log")
#' )
#'
#' # Retrieve fitted model
#' fit <- attr(res_logit, "model")
#' }
#'
#' @importFrom stats glm p.adjust nobs AIC
#' @importFrom broom tidy
#' @importFrom dplyr rename mutate select filter case_when
#' @importFrom rlang .data
#'
#' @export
get_glm <- function(
    data,
    outcome,
    predictors,
    family           = "binomial",
    adjust_method    = "BH",
    conf_level       = 0.95,
    exponentiate     = NULL,
    remove_intercept = TRUE,
    verbose          = FALSE
) {

  # -- 1. Input validation ----------------------------------------------------
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or tibble.", call. = FALSE)
  }

  if (!is.character(outcome) || length(outcome) != 1 || is.na(outcome)) {
    stop("`outcome` must be a single non-missing character string.",
         call. = FALSE)
  }

  if (!outcome %in% names(data)) {
    stop(sprintf("`outcome` column '%s' not found in `data`.", outcome),
         call. = FALSE)
  }

  if (!is.character(predictors) || length(predictors) < 1 ||
      any(is.na(predictors))) {
    stop("`predictors` must be a non-empty character vector with no missing values.",
         call. = FALSE)
  }

  missing_preds <- setdiff(predictors, names(data))

  if (length(missing_preds) > 0) {
    stop(sprintf(
      "The following predictors are not in `data`: %s",
      paste(missing_preds, collapse = ", ")
    ), call. = FALSE)
  }

  if (!is.numeric(conf_level) || length(conf_level) != 1 ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("`conf_level` must be a single numeric value between 0 and 1.",
         call. = FALSE)
  }

  if (!is.logical(remove_intercept) || length(remove_intercept) != 1 ||
      is.na(remove_intercept)) {
    stop("`remove_intercept` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1 || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(exponentiate) &&
      (!is.logical(exponentiate) || length(exponentiate) != 1 ||
       is.na(exponentiate))) {
    stop("`exponentiate` must be TRUE, FALSE, or NULL.", call. = FALSE)
  }

  adjust_method <- match.arg(
    adjust_method,
    choices = c("BH", "bonferroni", "holm", "BY", "fdr", "none")
  )

  # -- 2. Normalize family object ---------------------------------------------
  family_obj <- .normalize_glm_family(family, envir = parent.frame())

  fam_name  <- family_obj$family
  link_name <- family_obj$link

  # -- 3. Auto-set exponentiate based on link ---------------------------------
  if (is.null(exponentiate)) {
    exponentiate <- link_name %in% c("logit", "log")
  }

  # -- 4. Build formula and fit model -----------------------------------------
  outcome_quoted <- .quote_name(outcome)
  predictors_quoted <- vapply(predictors, .quote_name, character(1))

  formula_str <- paste(
    outcome_quoted,
    "~",
    paste(predictors_quoted, collapse = " + ")
  )

  model_formula <- stats::as.formula(formula_str)

  fit <- tryCatch(
    stats::glm(
      formula = model_formula,
      data    = data,
      family  = family_obj
    ),
    error = function(e) {
      stop(sprintf("GLM fitting failed: %s", conditionMessage(e)),
           call. = FALSE)
    }
  )

  if (verbose) {
    message("-- Model summary ------------------------------")
    print(summary(fit))
  }

  # -- 5. Tidy extraction via broom -------------------------------------------
  tidy_res <- tryCatch(
    broom::tidy(
      fit,
      conf.int     = TRUE,
      conf.level   = conf_level,
      exponentiate = exponentiate
    ),
    error = function(e) {
      stop(sprintf("Could not extract tidy model results: %s",
                   conditionMessage(e)), call. = FALSE)
    }
  )

  # Standardise column names
  tidy_res <- dplyr::rename(
    tidy_res,
    dplyr::all_of(c(
      std_error = "std.error",
      p_value   = "p.value",
      ci_lower  = "conf.low",
      ci_upper  = "conf.high"
    ))
  )

  # -- 6. Remove intercept before FDR -----------------------------------------
  if (remove_intercept) {
    tidy_res <- dplyr::filter(tidy_res, .data$term != "(Intercept)")
  }

  # -- 7. Multiple testing correction and significance stars ------------------
  tidy_res <- dplyr::mutate(
    tidy_res,
    q_value = stats::p.adjust(.data$p_value, method = adjust_method),
    significance = dplyr::case_when(
      .data$p_value < 0.001 ~ "***",
      .data$p_value < 0.01  ~ "**",
      .data$p_value < 0.05  ~ "*",
      .data$p_value < 0.10  ~ ".",
      TRUE                  ~ " "
    )
  )

  tidy_res <- dplyr::select(
    tidy_res,
    dplyr::all_of(c(
      "term", "estimate", "ci_lower", "ci_upper",
      "std_error", "statistic", "p_value", "q_value", "significance"
    ))
  )

  # -- 8. Attach model and metadata as attributes -----------------------------
  aic_value <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)

  attr(tidy_res, "model")         <- fit
  attr(tidy_res, "formula")       <- formula_str
  attr(tidy_res, "family")        <- fam_name
  attr(tidy_res, "link")          <- link_name
  attr(tidy_res, "n_obs")         <- stats::nobs(fit)
  attr(tidy_res, "AIC")           <- aic_value
  attr(tidy_res, "exponentiate")  <- exponentiate
  attr(tidy_res, "adjust_method") <- adjust_method

  class(tidy_res) <- c("get_glm_result", class(tidy_res))

  return(tidy_res)
}

# -- S3 print method ----------------------------------------------------------

#' @export
#' @noRd
print.get_glm_result <- function(x, digits = 3, ...) {

  fam    <- attr(x, "family")        %||% "?"
  link   <- attr(x, "link")          %||% "?"
  n      <- attr(x, "n_obs")         %||% "?"
  aic    <- attr(x, "AIC")           %||% NA_real_
  adj    <- attr(x, "adjust_method") %||% "BH"
  exp_tf <- attr(x, "exponentiate")

  required_cols <- c(
    "term", "estimate", "ci_lower", "ci_upper",
    "p_value", "q_value", "significance"
  )

  missing_cols <- setdiff(required_cols, names(x))

  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Cannot print `get_glm_result`: missing columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  label <- if (isTRUE(exp_tf)) {
    if (identical(link, "logit")) {
      "OR"
    } else if (identical(link, "log")) {
      "Ratio"
    } else {
      "exp(Estimate)"
    }
  } else {
    "Estimate"
  }

  aic_txt <- if (is.na(aic)) "NA" else sprintf("%.2f", aic)

  cat(sprintf(
    "\n-- get_glm result (%s, link = %s) -----------------------------\n",
    fam, link
  ))

  cat(sprintf(
    "   n = %s  |  AIC = %s  |  q-value correction: %s\n\n",
    n, aic_txt, adj
  ))

  # Important:
  # Drop custom S3 class before formatting/printing to avoid recursive dispatch.
  print_tbl <- as.data.frame(x, stringsAsFactors = FALSE)

  print_tbl$estimate <- sprintf("%.3f", print_tbl$estimate)
  print_tbl$ci <- sprintf(
    "[%.3f, %.3f]",
    print_tbl$ci_lower,
    print_tbl$ci_upper
  )

  print_tbl$p_value <- format.pval(
    print_tbl$p_value,
    digits = digits,
    eps = 0.001
  )

  print_tbl$q_value <- format.pval(
    print_tbl$q_value,
    digits = digits,
    eps = 0.001
  )

  out <- print_tbl[, c(
    "term", "estimate", "ci", "p_value", "q_value", "significance"
  ), drop = FALSE]

  names(out)[c(2, 3)] <- c(label, "95% CI")

  base::print.data.frame(out, row.names = FALSE, ...)

  cat("\n")

  invisible(x)
}
