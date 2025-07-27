#########################
# Function nice_Volcano #
#########################

#' Function to draw Volcano plots.
#'
#' Volcano plot with configurable point shapes and threshold annotations:
#' * Automatic triangle shapes for points above a user-defined y-axis limit (and a matching legend entry).
#' * Horizontal dashed line at one or more significance thresholds, annotated with its value.
#' * Vertical dashed lines at log-fold-change cutoffs, shown as custom x-axis ticks.
#'
#' @param results A data frame containing at least one column of effect sizes (e.g. log₂FC) and one column of significance (e.g. FDR).
#' @param x_var Name of the column in `results` to plot on the x-axis (e.g. log₂FC).
#' @param y_var Name of the column in `results` to plot on the y-axis (e.g. FDR).
#' @param label_var to be defined.
#' @param legend Logical. Control legend display. Default: TRUE.
#' @param title title.
#' @param colors colors.
#' @param x_range X-axis range of values.
#' @param y_max Maximum values of y-axis.
#' @param cutoff_x to be defined.
#' @param cutoff_y to be defined.
#' @param nice_x to be defined.
#' @param nice_y to be defined.
#' @param genes Vector of genes to label in the plot. Default: NULL.
#' @import ggplot2
#' @importFrom rlang .data
#' @export

nice_Volcano <- function(results, x_range = 9, y_max = 8, cutoff_y = 0.05, cutoff_x = 1,
                         nice_y = NULL, nice_x = NULL, y_var, x_var, label_var, legend = TRUE,
                         title, colors = c("red", "grey70", "blue"), genes = NULL)
{
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "Package \"ggrepel\" must be installed to use this function.",
      call. = FALSE
    )
  }

  # Fixed variables
  if (is.null(nice_y)) sig_nice <- cutoff_y*0.2 else sig_nice <- nice_y

  if (is.null(nice_x)) chg_nice <- cutoff_x*4 else chg_nice <- nice_x

  # Manage data
  d.volcano <- data.frame(results)
  d.volcano <- d.volcano[!is.na(d.volcano[, y_var]), ]
  d.volcano[d.volcano[, y_var] < 10**-y_max, y_var] <- 10**-y_max

  d.volcano$colors <- rep("other", nrow(d.volcano))
  d.volcano$colors[(d.volcano[, x_var] >= cutoff_x & d.volcano[, y_var] < cutoff_y)] <- "over"
  d.volcano$colors[(d.volcano[, x_var] <= -cutoff_x & d.volcano[, y_var] < cutoff_y)] <- "under"

  d.volcano$shapes <- rep("nohits", nrow(d.volcano))
  d.volcano$shapes[d.volcano[, y_var] <= 10**-y_max] <- "hits"

  # Condition for labeling
  if (is.null(genes)) {
    cond <- d.volcano[, y_var] < sig_nice & abs(d.volcano[, x_var]) > chg_nice
  } else {
    cond <- d.volcano[, label_var] %in% genes
  }

  y_min <- min(-log10(d.volcano[[y_var]]))

  # Plot
  p.volcano <- ggplot() + theme_bw() +
    scale_y_continuous(expand = c(0, 0.2),
                       limits = c(y_min, y_max),
                       breaks = seq(round(y_min, digits = 0), y_max, by = 1),
                       minor_breaks = seq(round(y_min, digits = 0), y_max, by = 0.5)) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(-x_range, x_range),
                       breaks = seq(-x_range, x_range, by = x_range*0.5),
                       minor_breaks = seq(-x_range, x_range, by = x_range*0.1)) +

    geom_point(data = d.volcano,
               aes(x = .data[[x_var]], y = -log10(.data[[y_var]]), fill = colors, shape = .data[["shapes"]]),
               size = 5.5, color = "gray10", alpha = 0.4, show.legend = legend) +

    # Vertical lines and labels:
    geom_vline(xintercept = c(-cutoff_x, cutoff_x), color = "grey20", alpha = 0.9, linetype = 2, linewidth = 1.2) +

    # Horizontal line and label:
    geom_hline(yintercept = -log10(cutoff_y), color = "grey20", alpha = 0.9, linetype = 2, linewidth = 1.2) +
    geom_text(aes(x = -x_range*0.7, y = -log10(cutoff_y) + 0.5, label = paste("q =", cutoff_y)), color = "grey20", size = 6.2) +

    # Color settings
    scale_fill_manual(name = "Expression",
                      values = c("over" = colors[1], "other" = colors[2], "under" = colors[3]),
                      breaks = c("over", "other", "under"),
                      labels = c("Up-regulated", "No difference", "Down-regulated"),
                      guide = guide_legend(override.aes = aes(shape = 21, size = 5, alpha = 1))) +

    # Change outliers shapes
    scale_shape_manual(name = "Shape",
                       values = c("hits" = 24, "nohits" = 21),
                       guide = guide_legend(override.aes = list(fill = "white", size = 5, alpha = 1)),
                       labels = c(bquote(FDR <= 10^-.(y_max)), bquote(FDR > 10^-.(y_max)))) +

    # Labels for detectable genes
    ggrepel::geom_label_repel(data = d.volcano[cond, ],
                     aes(x = .data[[x_var]], y = -log10(.data[[y_var]]), label = .data[[label_var]]),
                     inherit.aes = FALSE, parse = FALSE, max.iter = 5000, color = "black", size = 5,
                     segment.alpha = 1.5, box.padding = grid::unit(0.8, "lines"), min.segment.length = grid::unit(0.01, "lines")) +

    # Axis labels:
    labs(x = expression("log"[2]*"(Fold Change)"),
         y = expression("-log"[10]*"(FDR)"),
         title = title) +

    # Labels size
    theme(plot.title = element_text(size = 22, hjust = 0.5),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          legend.text = element_text(size = 14),
          legend.title = element_text(size= 16),
          legend.box.background = element_rect(color = "black"),
          legend.background = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.15, 0.75))

  return(p.volcano)
}
