#' Draw a Forest Plot from Cox Model Results
#'
#' Creates a publication-ready forest plot from a tidy table of model results,
#' such as the output produced by \code{get_cox()}.
#'
#' @param data A \code{data.frame} containing model results. By default, it is
#'   expected to contain \code{HR}, \code{CI_low}, \code{CI_high}, and
#'   \code{p.value} columns.
#' @param estimate_col Character. Name of the column containing the point
#'   estimate. Default: \code{"HR"}.
#' @param ci_low_col Character. Name of the lower confidence interval column.
#'   Default: \code{"CI_low"}.
#' @param ci_high_col Character. Name of the upper confidence interval column.
#'   Default: \code{"CI_high"}.
#' @param p_col Character. Name of the p-value column. Default:
#'   \code{"p.value"}.
#' @param label_col Character or \code{NULL}. Name of a precomputed label
#'   column. If \code{NULL}, labels are built automatically from
#'   \code{variable_col}, \code{term_col}, and \code{reference_col}.
#' @param variable_col Character. Name of the variable column. Default:
#'   \code{"variable"}.
#' @param term_col Character. Name of the cleaned term column. Default:
#'   \code{"term_clean"}.
#' @param reference_col Character. Name of the reference-level column.
#'   Default: \code{"reference"}.
#' @param p_display Numeric. Only rows with p-value <= \code{p_display} are
#'   shown. Default: \code{1} (show all rows).
#' @param title Character. Plot title. If \code{NULL}, a default title is used.
#' @param xlab Character. X-axis label. Default:
#'   \code{"Hazard Ratio (log scale)"}.
#' @param log_scale Logical. If \code{TRUE}, the x-axis is shown on a log10
#'   scale. Default: \code{TRUE}.
#' @param vline Numeric. Reference line position. Default: \code{1}.
#' @param color_sig Color for significant points (p < 0.05). Default:
#'   \code{"#c0392b"}.
#' @param color_ns Color for non-significant points. Default:
#'   \code{"#7f8c8d"}.
#' @param ref_line_color Color for the vertical reference line. Default:
#'   \code{"#2c3e50"}.
#' @param sort_by Character. How to order the plot rows. One of
#'   \code{"estimate"}, \code{"p.value"}, or \code{"input"}.
#'   Default: \code{"estimate"}.
#' @param point_size Numeric. Point size. Default: \code{3.5}.
#' @param base_size Numeric. Base font size for the ggplot theme.
#'   Default: \code{12}.
#' @param return_table Logical. If \code{TRUE}, returns a list with
#'   \code{$plot} and \code{$table}. If \code{FALSE}, returns only the plot.
#'   Default: \code{FALSE}.
#'
#' @return A \code{ggplot} object if \code{return_table = FALSE}, or a named
#'   list with \code{$plot} and \code{$table} if \code{return_table = TRUE}.
#'
#' @details
#' This function is intentionally model-agnostic. It only requires a table with
#' point estimates, confidence intervals, and p-values. For Cox models generated
#' with \code{get_cox()}, the default columns work without modification.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_vline scale_x_log10 scale_color_manual labs theme_minimal theme element_text   element_blank
#'
#' @examples
#' \dontrun{
#' cox_tab <- get_cox(
#'   data = df,
#'   time_col = "PFI.time",
#'   event_col = "PFI",
#'   vars = c("ER_Status_nature2012", "PAM50Call_RNAseq"),
#'   model = "univariable"
#' )
#'
#' nice_forest(cox_tab)
#'
#' nice_forest(
#'   cox_tab,
#'   p_display = 0.05,
#'   title = "PFI — Significant Cox Terms"
#' )
#' }
#'
#' @seealso
#' \code{\link{get_cox}} for fitting Cox proportional hazards models and
#' returning tidy results that can be passed directly to \code{nice_forest()}.
#'
#' \code{\link{nice_KM}} for Kaplan-Meier survival curve visualization.
#'
#' \code{\link[ggplot2]{ggplot}} for the underlying plotting system.
#'
#' @export
nice_forest <- function(data,
                        estimate_col   = "HR",
                        ci_low_col     = "CI_low",
                        ci_high_col    = "CI_high",
                        p_col          = "p.value",
                        label_col      = NULL,
                        variable_col   = "variable",
                        term_col       = "term_clean",
                        reference_col  = "reference",
                        p_display      = 1,
                        title          = NULL,
                        xlab           = "Hazard Ratio (log scale)",
                        log_scale      = TRUE,
                        vline          = 1,
                        color_sig      = "#c0392b",
                        color_ns       = "#7f8c8d",
                        ref_line_color = "#2c3e50",
                        sort_by        = c("estimate", "p.value", "input"),
                        point_size     = 3.5,
                        base_size      = 12,
                        return_table   = FALSE) {

  sort_by <- match.arg(sort_by)

  # --- 1. Validations ---
  if (!is.data.frame(data)) {
    stop("'data' should be a data.frame.")
  }

  required_cols <- c(estimate_col, ci_low_col, ci_high_col, p_col)
  missing_cols <- setdiff(required_cols, colnames(data))

  if (length(missing_cols) > 0) {
    stop(
      "The following required columns were not found in 'data': ",
      paste(missing_cols, collapse = ", ")
    )
  }

  is_single_number <- function(x) {
    is.numeric(x) && length(x) == 1 && !is.na(x)
  }

  if (!is_single_number(p_display) || p_display < 0 || p_display > 1) {
    stop("'p_display' should be a single numeric value between 0 and 1.")
  }

  if (!is_single_number(vline)) {
    stop("'vline' should be a single numeric value.")
  }

  # --- 2. Prepare plot table ---
  plot_data <- data

  plot_data$.estimate_plot <- suppressWarnings(as.numeric(plot_data[[estimate_col]]))
  plot_data$.ci_low_plot   <- suppressWarnings(as.numeric(plot_data[[ci_low_col]]))
  plot_data$.ci_high_plot  <- suppressWarnings(as.numeric(plot_data[[ci_high_col]]))
  plot_data$.p_plot        <- suppressWarnings(as.numeric(plot_data[[p_col]]))

  keep <- !is.na(plot_data$.estimate_plot) &
    !is.na(plot_data$.ci_low_plot) &
    !is.na(plot_data$.ci_high_plot) &
    !is.na(plot_data$.p_plot) &
    is.finite(plot_data$.estimate_plot) &
    is.finite(plot_data$.ci_low_plot) &
    is.finite(plot_data$.ci_high_plot) &
    plot_data$.p_plot <= p_display

  plot_data <- plot_data[keep, , drop = FALSE]

  if (nrow(plot_data) == 0) {
    stop("No rows remained after filtering by finite values and 'p_display'.")
  }

  if (isTRUE(log_scale)) {
    positive <- plot_data$.estimate_plot > 0 &
      plot_data$.ci_low_plot > 0 &
      plot_data$.ci_high_plot > 0 &
      vline > 0

    if (!all(positive)) {
      stop(
        "All estimates, confidence intervals, and 'vline' must be > 0 ",
        "when 'log_scale = TRUE'."
      )
    }
  }

  # --- 3. Build labels ---
  if (!is.null(label_col)) {
    if (!label_col %in% colnames(plot_data)) {
      stop("Column '", label_col, "' not found in 'data'.")
    }

    plot_data$.label_plot <- as.character(plot_data[[label_col]])

  } else if (all(c(variable_col, term_col, reference_col) %in% colnames(plot_data))) {
    plot_data$.label_plot <- paste0(
      plot_data[[variable_col]],
      " [ref: ",
      plot_data[[reference_col]],
      "]\n",
      plot_data[[term_col]],
      "\n",
      "HR=",
      round(plot_data$.estimate_plot, 2),
      " [",
      round(plot_data$.ci_low_plot, 2),
      "-",
      round(plot_data$.ci_high_plot, 2),
      "]",
      "  p=",
      formatC(plot_data$.p_plot, digits = 2, format = "e")
    )

  } else if (variable_col %in% colnames(plot_data)) {
    plot_data$.label_plot <- as.character(plot_data[[variable_col]])

  } else {
    plot_data$.label_plot <- rownames(plot_data)
  }

  # --- 4. Sort rows ---
  if (sort_by == "estimate") {
    plot_data <- plot_data[order(plot_data$.estimate_plot), , drop = FALSE]
  } else if (sort_by == "p.value") {
    plot_data <- plot_data[order(plot_data$.p_plot), , drop = FALSE]
  }

  plot_data$.label_plot <- factor(
    plot_data$.label_plot,
    levels = plot_data$.label_plot
  )

  plot_data$.significance_plot <- ifelse(
    plot_data$.p_plot < 0.05,
    "p < 0.05",
    "p >= 0.05"
  )

  # --- 5. Plot ---
  if (is.null(title)) {
    title <- "Forest Plot"
  }

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$.estimate_plot,
      y = .data$.label_plot,
      color = .data$.significance_plot
    )
  ) +
    ggplot2::geom_vline(
      xintercept = vline,
      linetype = "dashed",
      color = ref_line_color,
      linewidth = 0.8
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$.ci_low_plot, xmax = .data$.ci_high_plot),
      orientation = "y",
      width = 0.3,
      linewidth = 0.6
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(
      values = c("p < 0.05" = color_sig, "p >= 0.05" = color_ns),
      name = "Significance"
    ) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      axis.text.y = ggplot2::element_text(size = base_size - 4),
      legend.position = "top",
      panel.grid.major.y = ggplot2::element_blank()
    )

  if (isTRUE(log_scale)) {
    p <- p +
      ggplot2::scale_x_log10(
        breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16),
        labels = c("0.125", "0.25", "0.5", "1", "2", "4", "8", "16")
      )
  }

  if (return_table) {
    return(list(plot = p, table = plot_data))
  }

  return(p)
}
