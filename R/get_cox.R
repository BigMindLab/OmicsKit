#' Fit Cox Proportional Hazards Models and Return a Tidy Table
#'
#' Fits Cox proportional hazards models against a selected survival endpoint and
#' returns a tidy table with hazard ratios, confidence intervals, p-values, model
#' metadata, and labels ready to be plotted with \code{nice_forest()}.
#'
#' @param data A \code{data.frame} containing the survival columns
#'   (\code{time_col}, \code{event_col}) and the variables to test.
#' @param time_col Character. Name of the column with follow-up time.
#'   The column should be numeric or coercible to numeric. Default:
#'   \code{"PFI.time"}.
#' @param event_col Character. Name of the event indicator column.
#'   It should be coded as 0 = censored and 1 = event, or be coercible
#'   to numeric 0/1. Default: \code{"PFI"}.
#' @param vars Character vector. Variables of interest to test. If
#'   \code{NULL}, variables are selected automatically after excluding
#'   survival columns, ID-like columns, columns starting with \code{"_"},
#'   and optionally survival/outcome-like columns.
#' @param model Character. Type of Cox model to fit. One of
#'   \code{"univariable"}, \code{"multivariable"}, or \code{"adjusted"}.
#'   \code{"univariable"} fits one model per variable. \code{"multivariable"}
#'   fits one model containing all variables in \code{vars}. \code{"adjusted"}
#'   fits one model per variable of interest, adjusted by \code{adjust_vars}.
#' @param adjust_vars Character vector. Covariates used when
#'   \code{model = "adjusted"}. Ignored for \code{"univariable"} and
#'   \code{"multivariable"} models. Default: \code{NULL}.
#' @param keep_adjust_terms Logical. If \code{TRUE} and
#'   \code{model = "adjusted"}, adjustment covariate terms are kept in the
#'   returned table. If \code{FALSE}, only terms from the variable of interest
#'   are returned. Default: \code{FALSE}.
#' @param ref_levels A named list to explicitly set reference levels for
#'   categorical variables. Names must be column names and values must be the
#'   desired reference categories. Example:
#'   \code{list(ER_Status_nature2012 = "Negative",
#'              AJCC_Stage_nature2012 = "Stage I")}.
#' @param min_n Integer. Minimum total number of complete observations required
#'   for each Cox model. Default: \code{10}.
#' @param min_n_per_level Integer. For categorical predictors, minimum number
#'   of complete observations required in each category. Default: \code{5}.
#'   Set to \code{0} to disable.
#' @param min_events Integer. Minimum total number of events required for each
#'   Cox model. Default: \code{5}. Set to \code{0} to disable.
#' @param min_events_per_level Integer. For categorical predictors, minimum
#'   number of events required in each category. This is the main safeguard
#'   against sparse categories and infinite Cox coefficients. Default:
#'   \code{5}. Set to \code{0} to disable.
#' @param exclude_survival_like Logical. If \code{TRUE} and \code{vars = NULL},
#'   automatic variable selection excludes columns whose names suggest survival
#'   endpoints or outcome leakage, such as \code{OS}, \code{DSS}, \code{DFI},
#'   \code{PFI}, \code{days_to}, \code{death}, \code{vital_status},
#'   \code{followup}, or \code{survival}. Default: \code{TRUE}.
#' @param verbose Logical. If \code{TRUE}, prints informative messages when
#'   models or variables are skipped. Default: \code{TRUE}.
#'
#' @return A \code{data.frame} with one row per model term and columns including
#'   \code{model}, \code{model_id}, \code{variable}, \code{term},
#'   \code{term_clean}, \code{reference}, \code{HR}, \code{CI_low},
#'   \code{CI_high}, \code{p.value}, \code{n_used}, \code{n_events}, and
#'   \code{adjusted_for}.
#'
#' @details
#' The function performs several checks before fitting each model:
#' \enumerate{
#'   \item Removes empty strings and incomplete rows.
#'   \item Converts character and logical predictors to factors.
#'   \item Drops unused factor levels.
#'   \item Applies explicit reference levels via \code{ref_levels}.
#'   \item Skips sparse categorical variables using \code{min_n_per_level} and
#'         \code{min_events_per_level}.
#'   \item Skips models with possible perfect separation or infinite
#'         coefficients.
#'   \item Returns exponentiated Cox coefficients as hazard ratios.
#' }
#'
#' In this function, \code{"multivariable"} means a Cox model with multiple
#' predictors for one survival endpoint. This is different from
#' \code{"multivariate"}, which usually refers to multiple outcomes.
#'
#' @importFrom survival Surv coxph
#' @importFrom broom tidy
#' @importFrom dplyr bind_rows
#' @importFrom stats as.formula complete.cases relevel
#'
#' @examples
#' \dontrun{
#' # Univariable Cox models
#' cox_uni <- get_cox(
#'   data = df,
#'   time_col = "PFI.time",
#'   event_col = "PFI",
#'   vars = c("ER_Status_nature2012", "PAM50Call_RNAseq"),
#'   model = "univariable"
#' )
#'
#' # Multivariable Cox model
#' cox_multi <- get_cox(
#'   data = df,
#'   time_col = "PFI.time",
#'   event_col = "PFI",
#'   vars = c("age_at_initial_pathologic_diagnosis",
#'            "AJCC_Stage_nature2012",
#'            "PAM50Call_RNAseq"),
#'   model = "multivariable"
#' )
#'
#' # Adjusted Cox models
#' cox_adj <- get_cox(
#'   data = df,
#'   time_col = "PFI.time",
#'   event_col = "PFI",
#'   vars = c("ER_Status_nature2012", "HER2_Final_Status_nature2012"),
#'   adjust_vars = c("age_at_initial_pathologic_diagnosis",
#'                   "AJCC_Stage_nature2012"),
#'   model = "adjusted"
#' )
#'
#' nice_forest(cox_adj)
#' }

#' @seealso
#' \code{\link{nice_forest}} for plotting the tidy Cox model results returned
#' by \code{get_cox()}.
#'
#' \code{\link{nice_KM}} for Kaplan-Meier survival curve visualization.
#'
#' \code{\link[survival]{coxph}} and \code{\link[survival]{Surv}} for the
#' underlying Cox proportional hazards model and survival object.
#'
#' \code{\link[broom]{tidy}} for tidying model outputs.
#'
#' @export
get_cox <- function(data,
                    time_col              = "PFI.time",
                    event_col             = "PFI",
                    vars                  = NULL,
                    model                 = c("univariable", "multivariable", "adjusted"),
                    adjust_vars           = NULL,
                    keep_adjust_terms     = FALSE,
                    ref_levels            = NULL,
                    min_n                 = 10,
                    min_n_per_level       = 5,
                    min_events            = 5,
                    min_events_per_level  = 5,
                    exclude_survival_like = TRUE,
                    verbose               = TRUE) {

  model <- match.arg(model)

  # --- 1. Validations ---
  if (!is.data.frame(data)) {
    stop("'data' should be a data.frame.")
  }

  if (!time_col %in% colnames(data)) {
    stop("Column '", time_col, "' not found in 'data'.")
  }

  if (!event_col %in% colnames(data)) {
    stop("Column '", event_col, "' not found in 'data'.")
  }

  if (!is.null(vars) && !is.character(vars)) {
    stop("'vars' should be NULL or a character vector of column names.")
  }

  if (!is.null(adjust_vars) && !is.character(adjust_vars)) {
    stop("'adjust_vars' should be NULL or a character vector of column names.")
  }

  if (!is.null(ref_levels)) {
    if (!is.list(ref_levels) || is.null(names(ref_levels))) {
      stop("'ref_levels' should be NULL or a named list.")
    }
  }

  is_single_number <- function(x) {
    is.numeric(x) && length(x) == 1 && !is.na(x)
  }

  if (!is_single_number(min_n) || min_n < 1) {
    stop("'min_n' should be a single numeric value >= 1.")
  }

  if (!is_single_number(min_n_per_level) || min_n_per_level < 0) {
    stop("'min_n_per_level' should be a single numeric value >= 0.")
  }

  if (!is_single_number(min_events) || min_events < 0) {
    stop("'min_events' should be a single numeric value >= 0.")
  }

  if (!is_single_number(min_events_per_level) || min_events_per_level < 0) {
    stop("'min_events_per_level' should be a single numeric value >= 0.")
  }

  if (model == "adjusted" && (is.null(adjust_vars) || length(adjust_vars) == 0)) {
    stop("'adjust_vars' is required when model = 'adjusted'.")
  }

  # --- 2. Variable selection ---
  if (is.null(vars)) {
    all_cols <- colnames(data)

    id_like <- grep(
      "(^sample$|^sampleID$|sample_id|patient|barcode|^ID$|_ID$)",
      all_cols,
      ignore.case = TRUE,
      value = TRUE
    )

    survival_like <- character(0)

    if (isTRUE(exclude_survival_like)) {
      survival_like <- grep(
        "^(OS|DSS|DFI|PFI)(\\.time)?$|days_to|death|vital_status|last_follow|followup|survival",
        all_cols,
        ignore.case = TRUE,
        value = TRUE
      )
    }

    exclude_auto <- unique(c(
      time_col,
      event_col,
      all_cols[grepl("^_", all_cols)],
      id_like,
      survival_like,
      adjust_vars
    ))

    vars <- setdiff(all_cols, exclude_auto)
  } else {
    missing_vars <- setdiff(vars, colnames(data))

    if (length(missing_vars) > 0) {
      stop(
        "The following variables were not found in 'data': ",
        paste(missing_vars, collapse = ", ")
      )
    }

    vars <- unique(vars)
  }

  if (!is.null(adjust_vars)) {
    missing_adjust <- setdiff(adjust_vars, colnames(data))

    if (length(missing_adjust) > 0) {
      stop(
        "The following adjustment variables were not found in 'data': ",
        paste(missing_adjust, collapse = ", ")
      )
    }

    adjust_vars <- unique(adjust_vars)
  }

  if (length(vars) == 0) {
    stop("No variables available to test.")
  }

  # --- 3. Helper functions ---
  bt <- function(x) {
    paste0("`", x, "`")
  }

  skip_model <- function(model_id, reason) {
    if (isTRUE(verbose)) {
      message("Model '", model_id, "' skipped: ", reason)
    }
    return(NULL)
  }

  clean_term_label <- function(term, variable, reference) {
    out <- as.character(term)

    bt_var <- bt(variable)

    if (startsWith(out, bt_var)) {
      out <- substring(out, nchar(bt_var) + 1)
    } else if (startsWith(out, variable)) {
      out <- substring(out, nchar(variable) + 1)
    }

    out <- trimws(out)

    if (out %in% c("", "``") && identical(reference, "continuous")) {
      out <- "per 1-unit increase"
    } else if (out %in% c("", "``")) {
      out <- "level vs reference"
    }

    out
  }

  guess_variable <- function(term, candidate_vars) {
    candidate_vars <- candidate_vars[order(nchar(candidate_vars), decreasing = TRUE)]

    for (v in candidate_vars) {
      bt_v <- bt(v)

      if (identical(term, bt_v) ||
          startsWith(term, bt_v) ||
          identical(term, v) ||
          startsWith(term, v)) {
        return(v)
      }
    }

    NA_character_
  }

  make_formula <- function(model_vars) {
    rhs <- paste(bt(model_vars), collapse = " + ")

    paste0(
      "survival::Surv(", bt(time_col), ", ", bt(event_col), ") ~ ",
      rhs
    )
  }

  prepare_model_data <- function(model_vars, model_id) {
    model_vars <- unique(model_vars)
    needed_cols <- unique(c(time_col, event_col, model_vars))

    tmp <- data[, needed_cols, drop = FALSE]

    # Convert survival time to numeric
    tmp[[time_col]] <- suppressWarnings(as.numeric(as.character(tmp[[time_col]])))

    # Convert event indicator to numeric 0/1
    if (is.logical(tmp[[event_col]])) {
      tmp[[event_col]] <- as.integer(tmp[[event_col]])
    } else {
      tmp[[event_col]] <- suppressWarnings(as.numeric(as.character(tmp[[event_col]])))
    }

    # Convert empty strings to NA in character/factor predictors
    for (v in model_vars) {
      if (is.character(tmp[[v]]) || is.factor(tmp[[v]])) {
        empty_idx <- !is.na(tmp[[v]]) & tmp[[v]] == ""
        tmp[[v]][empty_idx] <- NA
      }
    }

    # Remove incomplete rows
    tmp <- tmp[stats::complete.cases(tmp), , drop = FALSE]

    # Remove non-finite survival values
    finite_idx <- is.finite(tmp[[time_col]]) &
      is.finite(tmp[[event_col]]) &
      tmp[[time_col]] >= 0

    # Remove non-finite numeric predictors
    for (v in model_vars) {
      if (is.numeric(tmp[[v]]) || is.integer(tmp[[v]])) {
        finite_idx <- finite_idx & is.finite(tmp[[v]])
      }
    }

    tmp <- tmp[finite_idx, , drop = FALSE]

    if (nrow(tmp) < min_n) {
      return(skip_model(
        model_id,
        paste0("fewer than ", min_n, " complete observations.")
      ))
    }

    event_values <- sort(unique(tmp[[event_col]]))

    if (!all(event_values %in% c(0, 1))) {
      return(skip_model(
        model_id,
        "event column should be coded as 0/1 after removing missing values."
      ))
    }

    n_events <- sum(tmp[[event_col]] == 1, na.rm = TRUE)

    if (min_events > 0 && n_events < min_events) {
      return(skip_model(
        model_id,
        paste0("fewer than ", min_events, " total events.")
      ))
    }

    for (v in model_vars) {
      # Convert character and logical predictors to factors
      if (is.character(tmp[[v]]) || is.logical(tmp[[v]])) {
        tmp[[v]] <- as.factor(tmp[[v]])
      }

      # Drop unused factor levels
      if (is.factor(tmp[[v]])) {
        tmp[[v]] <- droplevels(tmp[[v]])
      }

      # Skip categorical predictors with fewer than two levels
      if (is.factor(tmp[[v]]) && nlevels(tmp[[v]]) < 2) {
        return(skip_model(
          model_id,
          paste0("variable '", v, "' has fewer than two non-empty factor levels.")
        ))
      }

      # Skip non-categorical predictors with no variation
      if (!is.factor(tmp[[v]]) && length(unique(tmp[[v]])) < 2) {
        return(skip_model(
          model_id,
          paste0("variable '", v, "' has no variation.")
        ))
      }

      # Apply explicit reference level
      if (is.factor(tmp[[v]]) && !is.null(ref_levels) && !is.null(ref_levels[[v]])) {
        ref <- ref_levels[[v]]

        if (ref %in% levels(tmp[[v]])) {
          tmp[[v]] <- stats::relevel(tmp[[v]], ref = ref)
        } else {
          warning(
            "Reference '", ref, "' not found in variable '", v,
            "'. Available levels: ",
            paste(levels(tmp[[v]]), collapse = ", ")
          )
        }
      }

      # Check sample size per factor level
      if (is.factor(tmp[[v]])) {
        level_n <- table(tmp[[v]])

        if (min_n_per_level > 0 && any(level_n < min_n_per_level)) {
          small_levels <- names(level_n)[level_n < min_n_per_level]

          return(skip_model(
            model_id,
            paste0(
              "variable '", v, "' has level(s) with fewer than ",
              min_n_per_level, " observations: ",
              paste(small_levels, collapse = ", "),
              "."
            )
          ))
        }

        event_n <- tapply(
          tmp[[event_col]] == 1,
          tmp[[v]],
          sum,
          na.rm = TRUE
        )

        event_n <- event_n[levels(tmp[[v]])]
        event_n[is.na(event_n)] <- 0

        if (min_events_per_level > 0 && any(event_n < min_events_per_level)) {
          small_event_levels <- names(event_n)[event_n < min_events_per_level]

          return(skip_model(
            model_id,
            paste0(
              "variable '", v, "' has level(s) with fewer than ",
              min_events_per_level, " events: ",
              paste(small_event_levels, collapse = ", "),
              "."
            )
          ))
        }
      }
    }

    list(
      data = tmp,
      n_used = nrow(tmp),
      n_events = n_events
    )
  }

  fit_one_model <- function(model_vars,
                            model_id,
                            model_type,
                            target_vars = model_vars,
                            adjusted_for = character(0)) {
    prepared <- prepare_model_data(model_vars, model_id)

    if (is.null(prepared)) {
      return(NULL)
    }

    tmp <- prepared$data
    formula_str <- make_formula(model_vars)

    separation_warning <- FALSE

    separation_pattern <- paste(
      "Loglik converged before variable",
      "coefficient may be infinite",
      "ran out of iterations",
      sep = "|"
    )

    fit <- tryCatch(
      withCallingHandlers(
        survival::coxph(stats::as.formula(formula_str), data = tmp),
        warning = function(w) {
          warning_msg <- conditionMessage(w)

          if (grepl(separation_pattern, warning_msg, ignore.case = TRUE)) {
            separation_warning <<- TRUE
            invokeRestart("muffleWarning")
          }
        }
      ),
      error = function(e) {
        if (isTRUE(verbose)) {
          message(
            "Model '", model_id, "' skipped: Cox model could not be fitted. ",
            "Reason: ", conditionMessage(e)
          )
        }
        return(NULL)
      }
    )

    if (is.null(fit)) {
      return(NULL)
    }

    if (isTRUE(separation_warning)) {
      return(skip_model(
        model_id,
        "possible perfect separation or infinite coefficient."
      ))
    }

    res <- tryCatch(
      broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE),
      error = function(e) {
        if (isTRUE(verbose)) {
          message(
            "Model '", model_id, "' skipped: model results could not be tidied. ",
            "Reason: ", conditionMessage(e)
          )
        }
        return(NULL)
      }
    )

    if (is.null(res) || nrow(res) == 0) {
      return(skip_model(model_id, "no valid model results."))
    }

    required_cols <- c("term", "estimate", "std.error", "conf.low", "conf.high", "p.value")

    if (!all(required_cols %in% colnames(res))) {
      return(skip_model(
        model_id,
        "model output does not contain the required columns."
      ))
    }

    if (any(!is.finite(res$estimate)) ||
        any(!is.finite(res$std.error)) ||
        any(!is.finite(res$conf.low)) ||
        any(!is.finite(res$conf.high)) ||
        any(is.na(res$p.value))) {
      return(skip_model(
        model_id,
        "non-finite HR, SE, confidence interval, or p-value."
      ))
    }

    res$variable <- vapply(
      res$term,
      guess_variable,
      character(1),
      candidate_vars = model_vars
    )

    res <- res[!is.na(res$variable), , drop = FALSE]

    if (nrow(res) == 0) {
      return(skip_model(model_id, "terms could not be mapped to variables."))
    }

    if (model_type == "adjusted" && !isTRUE(keep_adjust_terms)) {
      res <- res[res$variable %in% target_vars, , drop = FALSE]
    }

    if (nrow(res) == 0) {
      return(skip_model(model_id, "no target-variable terms remained."))
    }

    res$reference <- vapply(
      res$variable,
      function(v) {
        if (is.factor(tmp[[v]])) {
          levels(tmp[[v]])[1]
        } else {
          "continuous"
        }
      },
      character(1)
    )

    res$term_clean <- mapply(
      clean_term_label,
      res$term,
      res$variable,
      res$reference,
      USE.NAMES = FALSE
    )

    names(res)[names(res) == "estimate"] <- "HR"
    names(res)[names(res) == "conf.low"] <- "CI_low"
    names(res)[names(res) == "conf.high"] <- "CI_high"

    res$model        <- model_type
    res$model_id     <- model_id
    res$n_used       <- prepared$n_used
    res$n_events     <- prepared$n_events
    res$adjusted_for <- if (length(adjusted_for) > 0) {
      paste(adjusted_for, collapse = " + ")
    } else {
      NA_character_
    }

    preferred_cols <- c(
      "model", "model_id", "variable", "term", "term_clean", "reference",
      "HR", "CI_low", "CI_high", "p.value", "std.error", "statistic",
      "n_used", "n_events", "adjusted_for"
    )

    other_cols <- setdiff(colnames(res), preferred_cols)

    res[, c(intersect(preferred_cols, colnames(res)), other_cols), drop = FALSE]
  }

  # --- 4. Fit requested model type ---
  if (model == "univariable") {
    results_list <- lapply(vars, function(v) {
      fit_one_model(
        model_vars = v,
        model_id = v,
        model_type = "univariable",
        target_vars = v
      )
    })

  } else if (model == "multivariable") {
    results_list <- list(
      fit_one_model(
        model_vars = vars,
        model_id = "multivariable",
        model_type = "multivariable",
        target_vars = vars
      )
    )

  } else {
    results_list <- lapply(vars, function(v) {
      adj <- setdiff(adjust_vars, v)

      fit_one_model(
        model_vars = unique(c(v, adj)),
        model_id = paste0(v, " adjusted"),
        model_type = "adjusted",
        target_vars = v,
        adjusted_for = adj
      )
    })
  }

  results <- dplyr::bind_rows(results_list)

  if (is.null(results) || nrow(results) == 0) {
    stop("No Cox model could be fitted. Please check the data and variables.")
  }

  rownames(results) <- NULL
  results
}
