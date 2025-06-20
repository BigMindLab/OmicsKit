########################
# Function nice_KM    #
########################

#' Function to generate Kaplan Meier survival plots for a given binary gene variable
#'
#' This function fits Kaplan Meier survival curves stratified by status (e.g., 'No' vs 'Yes') for a specified gene column in a data frame. It allows:
#' * Automatic placement of the p-value at the midpoint of the time range (if coord = NULL).
#' * Automatic handling of cases where only one category ('No' or 'Yes') is present, returning an 'empty' placeholder plot with a warning message.
#'
#' @param data Data.frame. Must contain columns specified by `gene`, `time_var`, and `event_var`.
#' @param gene String. Name of the column in `data` indicating mutation status (values 'No'/'Yes').
#' @param time_var String. Name of the column for survival time.
#' @param event_var String. Name of the column for event indicator (0/1).
#' @param coord Numeric or NULL. X-coordinate at which to display the log-rank p-value. If NULL, placed at midpoint of time range.
#' @param title_prefix String. Prefix for the legend title. Default `"Mut "`.
#' @param colors Character vector of length 2. Two colors for the plot lines. Default `c("#1F8FFF", "#ED4D4D")`.
#' @param conf_int Logical. Whether to show confidence intervals. Default `TRUE`.
#' @param risk_table Logical. Whether to include risk-table below the plot. Default `FALSE`.
#' @param legend_pos Numeric vector of length 2 or NULL. X/Y inside-plot coordinates for legend placement. Default `NULL`.
#' @param title_size Numeric. Font size for the plot title.
#' @param axis_text Named numeric vector. Sizes for axis tick labels.
#' @param axis_title_size Named numeric vector. Sizes for axis titles.
#' @param legend_text Named numeric vector. Sizes for legend title and labels.
#' @param pval_size Numeric. Font size for the p-value annotation.
#' @param returnData Logical. If `TRUE`, returns a list with `km_fit` and `plot`; if `FALSE`, returns only the `ggplot` object. Default `FALSE`.
#' @return A `ggplot` object (or a list with `km_fit` and `plot` if `returnData = TRUE`).
#' @export

nice_KM <- function(data, gene, time_var, event_var, coord = NULL, title_prefix = "Mut ", colors = c("#1F8FFF", "#ED4D4D"), 
                    conf_int = TRUE, risk_table = FALSE,legend_pos = NULL, title_size = 24, axis_text = c(x = 18, y = 18),
                    axis_title_size = c(x = 20, y = 20), legend_text = c(title = 18, labels = 16), pval_size = 6,
                    returnData = FALSE) 
{
  
  # Compute midpoint for p-value
  times <- data[[time_var]]
  if (is.null(coord)) {
    coord <- mean(range(times, na.rm = TRUE))
  }
  
  # Count categories and handle singleton
  counts <- table(data[[gene]])
  if (length(counts) < 2 || any(counts == 0)) {
    warning("Solo una categoria en ", gene)
    empty_plot <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate(
        "text", x = .5, y = .5, hjust = .5, size = 6,
        label = paste0("Solo una\ncategoria en\n", sub("_muts$", "", gene))
      )
    
    if (returnData) {
      return(list(km_fit = NULL, plot = empty_plot))
    } else {
      return(empty_plot)
    }
  }
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop(
      "Package \"survival\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  # Fit Kaplan Meier model
  km_fit <- eval(substitute(
    survival::survfit(
      survival::Surv(t, e) ~ g,
      data = data
    ),
    list(
      t = as.name(time_var),
      e = as.name(event_var),
      g = as.name(gene)
    )
  )
  )
  
  # Legend labels with counts
  legend_labs <- paste0(names(counts), " (n=", as.integer(counts), ")")
  
  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop(
      "Package \"survminer\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  # Generate survminer plot
  km_out <- survminer::ggsurvplot(
    km_fit,
    data = data,
    xlim = c(0, max(times, na.rm = TRUE) + 50),
    pval = TRUE,
    pval.size = pval_size,
    pval.coord = c(coord, 0.9),
    conf.int = conf_int,
    risk.table = risk_table,
    legend.title = paste0(title_prefix, sub("_muts$", "", gene)),
    legend.labs = legend_labs,
    palette = colors
  )
  
  # Extract ggplot object
  km_plot <- km_out$plot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold"),
      axis.title.x = ggplot2::element_text(size = axis_title_size["x"]),
      axis.title.y = ggplot2::element_text(size = axis_title_size["y"]),
      axis.title = ggplot2::element_text(size = title_size["axis"]),
      axis.text.x = ggplot2::element_text(size = axis_text["x"]),
      axis.text.y = ggplot2::element_text(size = axis_text["y"]),
      legend.title = ggplot2::element_text(size = legend_text["title"]),
      legend.text = ggplot2::element_text(size = legend_text["labels"]),
      legend.position.inside = legend_pos,
      legend.direction = "horizontal"
    )
  
  if (returnData) {return(list(km_fit = km_fit, plot = km_plot))}
  km_plot
}