#' Plot global GSEA results
#'
#' Generates a composite plot displaying NES values, pathway labels,
#' and a \emph{logFDR} legend, organized by MSigDB collections.
#'
#' @param data Data frame containing the GSEA results.
#' @param geneset_col Name of the column containing the genesets.
#' @param collection_col Name of the column containing the collections.
#' @param nes_col Name of the column containing the NES values.
#' @param logfdr_col Name of the column containing \eqn{-\log_{10}(FDR)} values.
#' @param output_path_base Path and base name for the output files (without extension).
#' @param width_output Plot width in inches.
#' @param height_output Plot height in inches.
#' @param text_size_genesets Text size for the geneset labels.
#' @param text_size_collection Text size for the collection labels.
#' @return Saves two files: \code{output_path_base}.jpg and \code{output_path_base}.pdf.
#' @import ggplot2 patchwork cowplot scales
#' @export
plot_GSEA <- function(data, geneset_col, collection_col, nes_col, logfdr_col, output_path_base, width_output, height_output, text_size_genesets = 5, text_size_collection = 5) {
  
  # Rename columns dynamically
  data <- data[, c(geneset_col, collection_col, nes_col, logfdr_col)]
  colnames(data) <- c("Geneset", "Collection", "NES", "logFDR")
  
  # Order data by NES value (descending)
  data <- data[order(data$NES, decreasing = TRUE), ]
  
  # Ensure Geneset and Collection are factors with ordered levels
  data$Geneset <- factor(data$Geneset, levels = rev(unique(data$Geneset)))
  data$Collection <- factor(data$Collection, levels = unique(data$Collection))
  
  # Right-side label: "MSigDB" vertically centered, in bold and italic
  plot_text_msigdb <- ggplot() + 
    annotate("text", label = "MSigDB", fontface = "bold.italic", angle = 90, size = 35, x = 0, y = 0.5)+
    theme_void()
  
  # Lef-side label: "Pathways" vertically centered, in bold and italic
  plot_text_pathways <- ggplot() + 
    annotate("text", label = "Pathways", fontface = "bold.italic", angle = 90, size = 35, x = 0, y = 0.5)+
    theme_void()
  
  # Right panel: Collection labels (without repetition)
  plot_right <- ggplot(data, aes(y = Geneset, x = 1.5, label = Collection)) +
    geom_text(aes(label = ifelse(duplicated(Collection), "", Collection)), 
              hjust = 0.5, size = 0, fontface = "bold") +  
    facet_grid(Collection ~ ., scales = "free_y", space = "free", switch = "y") + 
    theme_void() + 
    theme(strip.text.y = element_text(angle = 0, hjust = 1, size = text_size_collection),
          panel.spacing = unit(1, "lines"))
  
  # Center panel: NES bar plot
  plot_center <- ggplot(data, aes(x = NES, y = Geneset, fill = logFDR)) +
    geom_col(color = "black", size = 1) +
    scale_fill_gradient(low = "white", high = "red", 
                        limits = c(0,3), breaks = seq(0,3,1)) +
    scale_y_discrete(position = "right") +  
    facet_grid(Collection ~ ., scales = "free_y", space = "free_y") +
    theme_bw() +
    labs(x = "NES", y = "") +
    theme(axis.text.y = element_blank(),
          strip.background = element_rect(fill = "white", color = "black",linewidth = 1 ),
          axis.ticks.y = element_line(color = "black", size = 1.5),
          axis.ticks.length = unit(0.3, "cm"),
          strip.text.y = element_text(size = 1, margin = margin(0, 0, 0, 0)),# element_blank(),
          legend.position = "none",
          axis.title.x = element_text(size = 49),
          axis.text.x = element_text(size = 45),
          panel.spacing = unit(4, "lines")
    )
  
  # Left panel: Pathays labels
  plot_left <- ggplot(data, aes(y = Geneset, x = 0, label = Geneset)) +
    geom_text(hjust = 1, size = text_size_genesets) +  
    theme_void() + 
    theme(axis.text.y = element_blank(),
          plot.margin = margin(0, 0, 0, -50))
  
  # Legend panel
  plot_legend <- ggplot(data, aes(x = NES, y = Geneset, fill = logFDR)) +
    geom_tile() + 
    scale_fill_gradient(low = "white", high = "red",
                        name = expression(-log[10] ~ FDR),  # log10FDR with subscrip,
                        limits = c(0,3), breaks = seq(0,3,1),
                        guide = guide_colorbar(ticks.colour = "black",   # Make ticks black
                                               ticks.linewidth = 1.5,    # Make ticks thicker
                                               draw.ulim = TRUE,       # Draw upper limit tick
                                               draw.llim = TRUE)) +    # Draw lower limit tick  
    theme_bw() +
    theme(legend.position = "right",
          legend.box = "vertical",
          legend.title = element_text(size = 44, hjust = 0.5, face = "bold"),  # Bigger title
          legend.text = element_text(size = 30),  # Bigger legend text
          legend.key.size = unit(1.5, "cm"),  # Bigger color key size
          legend.key.height = unit(2, "cm"),  # Increase the height of the legend box
          legend.spacing = unit(3.5, "cm"),  # More space between title and legend
          legend.box.margin = margin(10, 20, 10, 10))#5, 5, 10, 5))  # Adjust internal spacing
  
  plot_legend <- plot_legend + theme(legend.box = "vertical")
  plot_right_legend <- get_legend(plot_legend)
  
  # Extract legend
  #plot_right_legend <- get_legend(plot_legend)
  
  # Combine all plots
  final_plot <- plot_text_pathways + plot_left +  plot_center + plot_right + plot_text_msigdb + plot_right_legend + 
    plot_layout(ncol = 6, widths = c(4, 25, 15, 3, 10, 3))  
  
  # Display the plot
  #print(final_plot)
  
  # Save the plot as JPG and PDF
  ggsave(paste0(output_path_base, ".jpg"), final_plot, width = width_output, height = height_output, units = "in", dpi = 300)
  ggsave(paste0(output_path_base, ".pdf"), final_plot, width = width_output, height = height_output, units = "in")
}
