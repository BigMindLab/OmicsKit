#' Create and save a customized barplot for GSEA results
#'
#' @param data A data frame containing GSEA results with columns such as `datatype`, `NES`, `-Log10FDR`, and `New_name`.
#' @param output_path The file path where the barplot will be saved (SVG format).
#' @param custom_labels A named vector of custom expressions for x-axis labels.
#' @param height Height of the saved plot in inches. Default: 49.
#' @param width Width of the saved plot in inches. Default: 30.
#' @param axis_y Name of the column to use for the y-axis aesthetic, as a string. Default: "NES".
#'
#' @details This function generates a customized barplot with grouped bars,
#' adjusted aesthetics, and personalized axis labels, saving the result as an SVG file.
#'
#' @import ggplot2
#' @importFrom scales squish
#' @importFrom rlang .data
#' @importFrom grid unit
#' @export
create_barplot <- function(data, output_path, custom_labels, height = 49, width = 30, axis_y = "NES") {
  # Generate the barplot
  barplot <- ggplot2::ggplot(data, aes(x = .data$datatype, y = .data[[axis_y]], fill = .data[["-Log10FDR"]])) +
    ggplot2::geom_bar(stat = "identity", color = "black", size = 0.5, width = 0.6) +
    ggplot2::scale_fill_gradient(
      low = "white", high = "red", na.value = "white",
      limits = c(0, 5.5),                # Fixed legend limits
      oob = scales::squish,              # Squish out-of-bounds values
      guide = ggplot2::guide_colorbar(barwidth = 3, barheight = 18)
    ) +
    ggplot2::labs(x = "Comparisons", y = axis_y) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.line.x     = ggplot2::element_blank(),
      axis.line       = ggplot2::element_line(color = "black", size = 0.5),
      axis.title.x    = ggplot2::element_text(size = 45, face = "bold", margin = ggplot2::margin(t = 25)),
      axis.title.y    = ggplot2::element_text(size = 45, face = "bold"),
      axis.text.x     = ggplot2::element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 30),
      axis.text.y.left= ggplot2::element_text(size = 50),
      legend.title    = ggplot2::element_text(size = 50),
      legend.text     = ggplot2::element_text(size = 50),
      panel.spacing   = unit(0.6, "lines"),
      panel.border    = ggplot2::element_blank(),
      strip.background= ggplot2::element_rect(fill = "white", color = "white"),
      strip.text.y.left = ggplot2::element_text(size = 50, angle = 0, hjust = 0.5),
      strip.placement = "outside"
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "black", size = 2) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::ylim(-3.4, 3.2) +
    ggplot2::facet_wrap(~ .data$New_name, ncol = 2, strip.position = "left", scales = "free_y") +
    ggplot2::scale_x_discrete(labels = custom_labels)
  
  # Save the barplot
  ggplot2::ggsave(filename = output_path, plot = barplot, height = height, width = width, limitsize = FALSE)
}
