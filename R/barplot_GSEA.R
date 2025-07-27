#########################
# Function barplot_GSEA #
#########################

#' Create and save a customized barplot for GSEA results
#'
#' This function generates a customized barplot with:
#' * Grouped bars.
#' * Adjusted aesthetics.
#' * Personalized axis labels.
#' * Optionally save the result in SVG format.
#'
#' @param data A data frame containing GSEA results with columns such as `datatype`, `NES`, `-Log10FDR`, and `New_name`.
#' @param output_path The file path where the barplot will be saved (SVG format).
#' @param custom_labels A named vector of custom expressions for x-axis labels.
#' @param axis_y Name of the column to use for the y-axis aesthetic, as a string. Default: "NES".
#' @import ggplot2
#' @importFrom rlang .data
#' @export

barplot_GSEA <- function(data, output_path, custom_labels, axis_y = "NES")

{
  # Generate the barplot
  barplot <- ggplot(data, aes(x = .data$datatype, y = .data[[axis_y]], fill = .data[["-Log10FDR"]])) +
    geom_bar(stat = "identity", color = "black", size = 0.5, width = 0.6) +
    scale_fill_gradient(
      low = "white", high = "red", na.value = "white",
      limits = c(0, 5.5),                # Fixed legend limits
      oob = scales::squish,              # Squish out-of-bounds values
      guide = guide_colorbar(barwidth = 3, barheight = 18)
    ) +
    labs(x = "Comparisons", y = axis_y) +
    theme_bw() +
    theme(
      axis.line.x = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      axis.title.x = element_text(size = 45, face = "bold", margin = margin(t = 25)),
      axis.title.y = element_text(size = 45, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 30),
      axis.text.y.left = element_text(size = 50),
      legend.title = element_text(size = 50),
      legend.text = element_text(size = 50),
      panel.spacing = grid::unit(0.6, "lines"),
      panel.border = element_blank(),
      strip.background = element_rect(fill = "white", color = "white"),
      strip.text.y.left = element_text(size = 50, angle = 0, hjust = 0.5),
      strip.placement = "outside"
    ) +
    geom_hline(yintercept = 0, color = "black", size = 2) +
    expand_limits(y = 0) +
    ylim(-3.4, 3.2) +
    facet_wrap(~ .data$New_name, ncol = 2, strip.position = "left", scales = "free_y") +
    scale_x_discrete(labels = custom_labels)

  return(barplot)
}
