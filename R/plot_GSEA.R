utils::globalVariables(c(
  "NAME","GENES","SIZE","tags","L.EDGE_size",
  "...12","numeric_cols","LEADING EDGE","signal",
  "FDR q-val","Log10FDR","FWER p-val","Comparison"
))

#' Unified GSEA plotting function with theme configuration
#'
#' Creates either a global GSEA plot or a faceted barplot depending on the number of unique
#' comparisons in the `Comparison` column. Allows customizing all previously hard-coded theme
#' parameters via a single `theme_params` list.
#'
#' @param data A data frame containing GSEA results.
#' @param Comparison Name of the column defining different comparisons. Default: "Comparison".
#' @param custom_labels Named vector of labels for the x-axis discrete scale (barplot mode only). Default: NULL.
#' @param axis_y Name of the column to use for the y-axis aesthetic. Default: "NES".
#' @param fdr_col Name of the column containing FDR values. Default: "FDR".
#' @param logFDR Logical; if TRUE, compute -log10(FDR) from `fdr_col`, otherwise use `fdr_col` directly. Default: TRUE.
#' @param geneset_col Name of the column containing the geneset labels (single comparison mode).
#' @param collection_col Name of the column containing the MSigDB collections (single comparison mode).
#' @param nes_col Name of the column containing NES values (single comparison mode).
#' @param logfdr_col Name of the column containing -log10(FDR) or similar (single comparison mode).
#' @param order One of "desc" or "asc"; order of `axis_y` values. Default: c("desc","asc").
#' @param ncol_wrap Number of columns for `facet_wrap` in barplot mode. Default: 2.
#' @param free_y Logical; if TRUE, allow free y scales in facets. Default: TRUE.
#' @param fill_limits Numeric vector of length 2 to set fill gradient limits (barplot mode). Default: NULL.
#' @param fill_palette Character vector of two colors for fill gradient. Default: c("white","red").
#' @param theme_params Named list to override default theme parameters (see details).
#' @details theme_params may include:
#'   \describe{
#'     \item{side_label_size}{Size for side panel labels (default 35)}
#'     \item{geneset_text_size}{Text size for geneset labels (default 5)}
#'     \item{collection_text_size}{Text size for collection labels (default 5)}
#'     \item{panel_widths}{Patchwork widths vector (default c(4,25,15,3,10,3))}
#'     \item{bar_col}{Bar/col border color (default "black")}  
#'     \item{bar_size}{Border size for bars (default 0.5)}
#'     \item{bar_width}{Width for bars (default 0.6)}
#'     \item{col_size}{Border size for geom_col (default 1)}
#'     \item{hline_size}{Size for horizontal line at y=0 (default 2)}
#'     \item{axis_title_size}{Font size for axis titles (default 45)}
#'     \item{axis_text_size_x}{Font size for x-axis text (default 30)}
#'     \item{axis_text_size_y}{Font size for y-axis text (default 50)}
#'     \item{tick_size}{Size for axis ticks (default 1.5)}
#'     \item{tick_length}{Length for axis ticks in cm (default 0.3)}
#'     \item{strip_text_size}{Font size for strip text (default 50)}
#'     \item{panel_spacing_single}{Facet spacing single mode (default 4)}
#'     \item{panel_spacing_multi}{Facet spacing multi mode (default 0.6)}
#'   }
#' @return A ggplot or patchwork object for the GSEA plot.
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom patchwork plot_layout
#' @importFrom cowplot get_legend
#' @importFrom utils modifyList
#' @export
plot_GSEA <- function(
    data,
    Comparison = "Comparison",
    custom_labels = NULL,
    axis_y = "NES",
    fdr_col = "FDR",
    logFDR = TRUE,
    geneset_col,
    collection_col,
    nes_col,
    logfdr_col,
    order = c("desc", "asc"),
    ncol_wrap = 2,
    free_y = TRUE,
    fill_limits = NULL,
    fill_palette = c("white", "red"),
    theme_params = list()
) {
  defaults <- list(
    side_label_size      = 35,
    geneset_text_size    = 5,
    collection_text_size = 5,
    panel_widths         = c(4,25,15,3,10,3),
    bar_col              = "black",
    bar_size             = 0.5,
    bar_width            = 0.6,
    col_size             = 1,
    hline_size           = 2,
    axis_title_size      = 45,
    axis_text_size_x     = 30,
    axis_text_size_y     = 50,
    tick_size            = 1.5,
    tick_length          = 0.3,
    strip_text_size      = 50,
    panel_spacing_single = 4,
    panel_spacing_multi  = 0.6
  )
  params <- utils::modifyList(defaults, theme_params)
  if (logFDR) data$logFDR <- -log10(data[[fdr_col]]) else data$logFDR <- data[[fdr_col]]
  order <- match.arg(order)
  data <- data[order(data[[axis_y]], decreasing = (order == "desc")), ]
  if (length(unique(data[[Comparison]])) == 1) {
    if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork required", call. = FALSE)
    if (!requireNamespace("cowplot", quietly = TRUE)) stop("cowplot required", call. = FALSE)
    df <- data[, c(geneset_col, collection_col, nes_col, logfdr_col)]
    colnames(df) <- c("Geneset", "Collection", "NES", "logFDR")
    df$Geneset    <- factor(df$Geneset, levels = rev(unique(df$Geneset)))
    df$Collection <- factor(df$Collection, levels = unique(df$Collection))
    plot_text_pathways <- ggplot() +
      annotate("text", label = "Pathways", fontface = "bold.italic", angle = 90,
               size = params$side_label_size, x = 0, y = 0.5) + theme_void()
    plot_left <- ggplot(df, aes(y = .data$Geneset, x = 0, label = .data$Geneset)) +
      geom_text(hjust = 1, size = params$geneset_text_size) + theme_void() +
      theme(axis.text.y = element_blank(), plot.margin = margin(0, 0, 0, -50))
    plot_center <- ggplot(df, aes(x = .data$NES, y = .data$Geneset, fill = .data$logFDR)) +
      geom_col(color = params$bar_col, size = params$col_size) +
      scale_fill_gradient(low = fill_palette[1], high = fill_palette[2],
                          limits = fill_limits, breaks = scales::pretty_breaks()) +
      scale_y_discrete(position = "right") + facet_grid(Collection ~ ., scales = "free_y", space = "free_y") +
      theme_bw() + labs(x = "NES", y = "") +
      theme(axis.text.y = element_blank(), strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
            axis.ticks.y = element_line(size = params$tick_size),
            axis.ticks.length = grid::unit(params$tick_length, "cm"),
            strip.text.y = element_text(size = 1), legend.position = "none",
            axis.title.x = element_text(size = params$axis_title_size),
            axis.text.x = element_text(size = params$axis_text_size_x),
            panel.spacing = grid::unit(params$panel_spacing_single, "lines"))
    plot_text_msigdb <- ggplot() +
      annotate("text", label = "MSigDB", fontface = "bold.italic", angle = 90,
               size = params$side_label_size, x = 0, y = 0.5) + theme_void()
    plot_right <- ggplot(df, aes(y = .data$Geneset, x = 1.5, label = .data$Collection)) +
      geom_text(aes(label = ifelse(duplicated(.data$Collection), "", .data$Collection)),
                hjust = 0.5, size = params$collection_text_size, fontface = "bold") +
      facet_grid(Collection ~ ., scales = "free_y", space = "free", switch = "y") + theme_void() +
      theme(strip.text.y = element_text(size = params$collection_text_size),
            panel.spacing = grid::unit(1, "lines"))
    plot_legend <- ggplot(df, aes(x = .data$NES, y = .data$Geneset, fill = .data$logFDR)) +
      geom_tile() + scale_fill_gradient(low = fill_palette[1], high = fill_palette[2],
                                        name = expression(-log[10] ~ FDR), limits = fill_limits,
                                        guide = guide_colorbar(ticks.colour = "black", ticks.linewidth = 1.5,
                                                               draw.ulim = TRUE, draw.llim = TRUE)) + theme_bw() +
      theme(legend.title = element_text(size = 44, face = "bold"),
            legend.text = element_text(size = 30),
            legend.key.size = grid::unit(1.5, "cm"),
            legend.key.height = grid::unit(2, "cm"),
            legend.spacing = grid::unit(3.5, "cm"),
            legend.box.margin = margin(10, 20, 10, 10))
    plot_right_legend <- cowplot::get_legend(plot_legend)
    final_plot <- plot_text_pathways + plot_left + plot_center + plot_right +
      plot_text_msigdb + plot_right_legend +
      patchwork::plot_layout(ncol = 6, widths = params$panel_widths)
  } else {
    final_plot <- ggplot(data, aes(x = .data[[Comparison]], y = .data[[axis_y]], fill = .data$logFDR)) +
      geom_bar(stat = "identity", color = params$bar_col, size = params$bar_size,
               width = params$bar_width) +
      scale_fill_gradient(low = fill_palette[1], high = fill_palette[2],
                          limits = fill_limits, oob = scales::squish,
                          guide = guide_colorbar(barwidth = 3, barheight = 18)) +
      labs(x = "Comparisons", y = axis_y) + theme_bw() +
      theme(axis.line.x = element_blank(), axis.line = element_line(size = 0.5),
            axis.title.x = element_text(size = params$axis_title_size),
            axis.title.y = element_text(size = params$axis_title_size),
            axis.text.x = element_text(size = params$axis_text_size_x),
            axis.text.y = element_text(size = params$axis_text_size_y),
            axis.ticks = element_line(size = params$tick_size),
            axis.ticks.length = grid::unit(params$tick_length, "cm"),
            strip.text = element_text(size = params$strip_text_size),
            panel.spacing = grid::unit(params$panel_spacing_multi, "lines")) +
      geom_hline(yintercept = 0, size = params$hline_size) + expand_limits(y = 0) +
      facet_wrap(~ .data$New_name, ncol = ncol_wrap,
                 scales = if (free_y) "free_y" else "fixed")
    if (!is.null(custom_labels)) final_plot <- final_plot +
        scale_x_discrete(labels = custom_labels)
  }
  return(final_plot)
}
