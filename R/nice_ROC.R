#####################
# Function nice_ROC #
#####################

#' Plot and Compare ROC Curves for Binary Classification Models
#'
#' @description
#' Generates publication-ready ROC curves for one or more binary classification
#' models evaluated on a held-out dataset. For each model the function computes
#' AUC with a 95% CI using the DeLong method. When exactly two models are
#' supplied, a DeLong paired test is performed to compare their AUCs and the
#' result is annotated on the plot.
#'
#' Models can be supplied either as fitted \code{glm} objects or as named
#' numeric vectors of predicted probabilities, making the function flexible
#' for any classifier that outputs scores.
#'
#' @param models
#'   A \emph{named} list. Each element is one of:
#'   \itemize{
#'     \item A fitted \code{glm} object (predictions generated internally via
#'           \code{predict(model, newdata = data, type = "response")}).
#'     \item A numeric vector of predicted probabilities of length
#'           \code{nrow(data)} (values in \[0, 1\]).
#'   }
#'   Names are used as legend labels (e.g.,
#'   \code{list("Clinical" = fit_A, "Clinical + Omics" = fit_B)}).
#' @param data
#'   A \code{data.frame} or \code{tibble}. Test / evaluation dataset. Required
#'   when any element of \code{models} is a \code{glm} object.
#' @param outcome
#'   Character (length 1). Name of the binary outcome column in \code{data}
#'   (must be \code{0}/\code{1} integer or numeric; \code{1} = event / positive
#'   class).
#' @param colors
#'   Character vector of colour codes for the ROC curves. Recycled if shorter
#'   than the number of models.
#'   Default: a four-colour colorblind-friendly palette
#'   (\code{c("#E63946", "#457B9D", "#2A9D8F", "#E9C46A")}).
#' @param show_ci
#'   Logical. Whether to display the DeLong 95% CI for each AUC in the legend.
#'   Default: \code{TRUE}.
#' @param show_delong
#'   Logical. When exactly two models are provided, annotate the plot with the
#'   DeLong test p-value and delta AUC. Default: \code{TRUE}.
#' @param smooth
#'   Logical. Apply kernel smoothing to the ROC curves
#'   (\code{pROC::smooth}). Useful for small samples; may hide real variability
#'   on large samples. Default: \code{FALSE}.
#' @param direction
#'   Character. Direction for \code{pROC::roc}. Usually \code{"auto"}
#'   (default). Set to \code{"<"} or \code{">"} if you need to enforce a
#'   direction.
#' @param plot_title
#'   Character string. Main title of the figure.
#'   Default: \code{"ROC Curve Comparison"}.
#' @param plot_subtitle
#'   Character string or \code{NULL}. Subtitle. Default: \code{NULL}.
#' @param theme_fn
#'   A \code{ggplot2} theme function (without parentheses). Applied to the
#'   plot. Default: \code{ggplot2::theme_classic}.
#' @param return_data
#'   Logical. If \code{TRUE}, returns a named list containing the ggplot object
#'   \emph{and} the underlying data/statistics. If \code{FALSE} (default), only
#'   the \code{ggplot} object is returned.
#'
#' @return
#' When \code{return_data = FALSE} (default): a \code{ggplot2} object.
#'
#' When \code{return_data = TRUE}: a named list with elements:
#' \describe{
#'   \item{\code{plot}}{The \code{ggplot2} ROC figure.}
#'   \item{\code{auc_table}}{A \code{tibble} with one row per model containing:
#'     \code{model}, \code{auc}, \code{ci_lower}, \code{ci_upper},
#'     \code{n_cases}, \code{n_controls}.}
#'   \item{\code{delong_test}}{Result of \code{pROC::roc.test} (DeLong) if
#'     exactly two models were supplied; otherwise \code{NULL}.}
#'   \item{\code{roc_objects}}{Named list of \code{pROC::roc} objects for
#'     further downstream analysis.}
#' }
#'
#' @details
#' ## AUC confidence intervals
#' CIs are computed with the DeLong method via \code{pROC::ci.auc} which
#' accounts for the correlation between paired measurements (same subjects).
#'
#' ## DeLong test
#' When two models are compared on the \emph{same} test set the observations
#' are paired, so the DeLong paired test is used
#' (\code{pROC::roc.test(method = "delong", paired = TRUE)}).
#'
#' ## Diagonal reference line
#' The grey dashed diagonal represents random performance (AUC = 0.50).
#'
#' @seealso \code{\link{get_glm}} for fitting the models fed into
#'   \code{nice_ROC}.
#'
#'
#'
#' @examples
#' \dontrun{
#' # -- Passing glm objects (predictions generated internally) -----------------
#' roc_result <- nice_ROC(
#'   models      = list("Clinical only"    = fit_A,
#'                      "Clinical + Omics" = fit_B),
#'   data        = test_data,
#'   outcome     = "stage_advanced",
#'   return_data = TRUE
#' )
#'
#' roc_result$plot
#' roc_result$auc_table
#' roc_result$delong_test
#'
#' # -- Passing probability vectors directly -----------------------------------
#' prob_A <- predict(fit_A, newdata = test_data, type = "response")
#' prob_B <- predict(fit_B, newdata = test_data, type = "response")
#'
#' nice_ROC(
#'   models  = list("Clinical" = prob_A, "Clinical + Omics" = prob_B),
#'   data    = test_data,
#'   outcome = "stage_advanced"
#' )
#'
#' # -- Single model (no comparison) -------------------------------------------
#' nice_ROC(
#'   models      = list("Clinical only" = fit_A),
#'   data        = test_data,
#'   outcome     = "stage_advanced",
#'   show_delong = FALSE
#' )
#' }
#'
#' @importFrom pROC roc auc ci.auc roc.test smooth
#' @importFrom ggplot2 ggplot aes geom_line geom_abline scale_color_manual labs theme theme_classic element_text element_rect annotate
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @importFrom stats predict
#'
#' @export
nice_ROC <- function(
    models,
    data         = NULL,
    outcome,
    colors       = c("#E63946", "#457B9D", "#2A9D8F", "#E9C46A"),
    show_ci      = TRUE,
    show_delong  = TRUE,
    smooth       = FALSE,
    direction    = "auto",
    plot_title   = "ROC Curve Comparison",
    plot_subtitle = NULL,
    theme_fn     = ggplot2::theme_classic,
    return_data  = FALSE
) {

  # -- 1. Validate inputs -----------------------------------------------------
  if (!is.list(models) || is.null(names(models))) {
    stop("`models` must be a *named* list.", call. = FALSE)
  }
  if (!is.character(outcome) || length(outcome) != 1) {
    stop("`outcome` must be a single character string.", call. = FALSE)
  }
  if (!is.null(data) && !outcome %in% names(data)) {
    stop(sprintf(
      "`outcome` column '%s' not found in `data`.", outcome
    ), call. = FALSE)
  }

  # -- 2. Extract true labels -------------------------------------------------
  # Determine true labels - either from data or from the first model element
  if (!is.null(data)) {
    true_labels <- data[[outcome]]
  } else {
    # All elements must be probability vectors of the same length; we need
    # true_labels supplied separately - error out.
    stop(
      "When `models` contains probability vectors, `data` must be provided ",
      "so that `outcome` labels can be extracted.",
      call. = FALSE
    )
  }

  if (!all(true_labels %in% c(0, 1, NA))) {
    stop(
      "`outcome` must be a binary 0/1 variable.",
      call. = FALSE
    )
  }

  # -- 3. Compute predicted probabilities for each model ---------------------
  pred_list <- lapply(seq_along(models), function(i) {
    mod  <- models[[i]]
    name <- names(models)[i]
    if (inherits(mod, "glm")) {
      if (is.null(data)) {
        stop(sprintf(
          "Model '%s' is a glm object but `data` is NULL.", name
        ), call. = FALSE)
      }
      probs <- stats::predict(mod, newdata = data, type = "response")
    } else if (is.numeric(mod)) {
      if (length(mod) != nrow(data)) {
        stop(sprintf(
          "Probability vector for model '%s' has length %d, expected %d.",
          name, length(mod), nrow(data)
        ), call. = FALSE)
      }
      probs <- mod
    } else {
      stop(sprintf(
        "Model '%s' must be a glm object or numeric probability vector.", name
      ), call. = FALSE)
    }
    return(probs)
  })
  names(pred_list) <- names(models)

  # -- 4. Build pROC::roc objects ---------------------------------------------
  roc_list <- lapply(names(pred_list), function(nm) {
    roc_obj <- pROC::roc(
      response  = true_labels,
      predictor = pred_list[[nm]],
      direction = direction,
      quiet     = TRUE
    )
    if (smooth) {
      roc_obj <- pROC::smooth(roc_obj, method = "density")
    }
    return(roc_obj)
  })
  names(roc_list) <- names(models)

  # -- 5. Compute AUC + DeLong 95% CI ----------------------------------------
  auc_rows <- lapply(names(roc_list), function(nm) {
    roc_obj <- roc_list[[nm]]
    auc_val <- as.numeric(pROC::auc(roc_obj))
    ci_val  <- pROC::ci.auc(roc_obj, method = "delong", conf.level = 0.95)
    tibble::tibble(
      model      = nm,
      auc        = round(auc_val, 4),
      ci_lower   = round(as.numeric(ci_val[1]), 4),
      ci_upper   = round(as.numeric(ci_val[3]), 4),
      n_cases    = sum(true_labels == 1, na.rm = TRUE),
      n_controls = sum(true_labels == 0, na.rm = TRUE)
    )
  })
  auc_table <- dplyr::bind_rows(auc_rows)

  # -- 6. DeLong test (only when exactly 2 models) ---------------------------
  delong_result <- NULL
  delong_label  <- NULL
  if (length(roc_list) == 2 && show_delong) {
    delong_result <- pROC::roc.test(
      roc1   = roc_list[[1]],
      roc2   = roc_list[[2]],
      method = "delong",
      paired = TRUE
    )
    delta_auc <- round(
      as.numeric(pROC::auc(roc_list[[1]])) -
        as.numeric(pROC::auc(roc_list[[2]])),
      4
    )
    p_delong <- delong_result$p.value
    delong_label <- sprintf(
      "DeLong test\n\u0394AUC = %.4f\np = %s",
      delta_auc,
      format.pval(p_delong, digits = 3, eps = 0.001)
    )
  }

  # -- 7. Build long data.frame for ggplot -----------------------------------
  roc_df <- dplyr::bind_rows(lapply(names(roc_list), function(nm) {
    roc_obj <- roc_list[[nm]]
    auc_row <- auc_table[auc_table$model == nm, ]

    if (show_ci) {
      legend_label <- sprintf(
        "%s\nAUC = %.3f [%.3f\u2013%.3f]",
        nm, auc_row$auc, auc_row$ci_lower, auc_row$ci_upper
      )
    } else {
      legend_label <- sprintf("%s  (AUC = %.3f)", nm, auc_row$auc)
    }

    data.frame(
      fpr   = 1 - roc_obj$specificities,
      tpr   = roc_obj$sensitivities,
      model = legend_label,
      stringsAsFactors = FALSE
    )
  }))

  # Preserve legend order
  roc_df$model <- factor(roc_df$model, levels = unique(roc_df$model))

  # -- 8. Recycle colours if needed ------------------------------------------
  n_models   <- length(models)
  pal_colors <- rep_len(colors, n_models)
  names(pal_colors) <- levels(roc_df$model)

  # -- 9. Build ggplot -------------------------------------------------------
  p <- ggplot2::ggplot(
    roc_df,
    ggplot2::aes(x = .data$fpr, y = .data$tpr, color = .data$model)
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype  = "dashed",
      color     = "grey55",
      linewidth = 0.5
    ) +
    ggplot2::geom_line(linewidth = 1.1, alpha = 0.92) +
    ggplot2::scale_color_manual(values = pal_colors) +
    ggplot2::labs(
      title    = plot_title,
      subtitle = plot_subtitle,
      x        = "False Positive Rate (1 \u2212 Specificity)",
      y        = "True Positive Rate (Sensitivity)",
      color    = NULL
    ) +
    theme_fn() +
    ggplot2::theme(
      legend.position  = c(0.75, 0.20),
      legend.text      = ggplot2::element_text(size = 9),
      legend.background = ggplot2::element_rect(
        fill    = "white",
        colour  = "grey80",
        linewidth = 0.4
      ),
      plot.title    = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey40"),
      axis.title    = ggplot2::element_text(size = 11)
    )

  # -- 10. Annotate DeLong test result ---------------------------------------
  if (!is.null(delong_label)) {
    p <- p +
      ggplot2::annotate(
        "label",
        x     = 0.68,
        y     = 0.35,
        label = delong_label,
        size  = 3.2,
        color = "grey20",
        fill  = "white",
        linewidth = 0.3
      )
  }

  # -- 11. Return ------------------------------------------------------------
  if (!return_data) {
    return(p)
  }

  return(list(
    plot         = p,
    auc_table    = auc_table,
    delong_test  = delong_result,
    roc_objects  = roc_list
  ))
}
