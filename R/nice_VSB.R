#####################
# Function nice_VSB #
#####################

#' Function to make Violin-Scatter-Box plots.
#'
#' This function will make a Boxplot, using a DEseq object.
#' It will show the data points on top with a small deviation (jitter) for a better visualization.
#'
#' @param object A DEseq object already transformed with the variance stabilizing or rlog transformations.
#' @param annotations Data frame with annotations.
#' @param variables To indicate the variables to be used as Shape and Fill of the markers.
#' @param genename The gene name to be used for the plot.
#' @param symbol The gene symbol to be used for the plot.
#' @param labels A vector containing the x-labels of the box-plot. Default: c("N", "P", "R", "M").
#' @param categories A vector containing the labels for the legend. Default: c("normal", "primary", "recurrence", "metastasis").
#' @param colors Vector of colors to be used for the categories of the variable assigned as Marker Fill.
#' @param shapes Vector of shapes to be used for the categories of the variable assigned as Marker Shape.
#' @param markersize Size of the marker.
#' @param alpha Transparency of the marker, which goes from 0 (transparent) to 1 (no transparent). Default: 0.8.
#' @param jitter Random deviation added to the dots. Default: 0.2.
#' @param title_size Font of the title and axis names. Default: c(axis = 20, fig = 24).
#' @param label_size Font of the labels (x-axis) and numbers (y-axis). Default: c(x = 20, y = 16).
#' @param legend_size Font of the title and elements of the legend. Default: c(title = 14, elements = 12).
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

nice_VSB <- function (object = NULL, annotations, variables = c(fill = "VarFill", shape = "VarShape"),
                      genename = NULL, symbol = NULL, labels = c("N", "P", "R", "M"),
                      categories = c("normal", "primary", "recurrence", "metastasis"),
                      colors = NULL, shapes = NULL, markersize = NULL, alpha = 0.8, jitter = 0.2,
                      title_size = c(axis = 20, fig = 24), label_size = c(x = 20, y = 16),
                      legend_size = c(title = 14, elements = 12))
{
  # Extracting the vector of counts for that gene
  gene_counts <- object[genename, ]
  log2_gc <- log2(gene_counts)

  # Making a dataframe for the plot
  df.box <- cbind(annotations, log2_gc)

  # Re-ordering sample_type for the plot
  df.box[, "sample_type"] <- factor(df.box[, "sample_type"],
                                    levels = categories,
                                    labels = labels)

  # Plot
  p.bs <- ggplot(df.box, aes(x = .data$sample_type, y = log2_gc)) + theme_bw() +
    geom_violin(alpha = 0.1, scale = "width", fill = "yellow", color = "peru",
		show.legend = FALSE, trim = TRUE) +
    geom_boxplot(width = 0.6, fill = "gray90") +
    labs(title = paste("Gene:", genename, symbol),
         x = expression("Sample Type"),
         y = expression("log"[2]*"(Normalized Gene Counts)")) +
    theme(plot.title = element_text(size = title_size["fig"]),
          axis.title = element_text(size = title_size["axis"]),
          axis.text.x = element_text(size = label_size["x"]),
          axis.text.y = element_text(size = label_size["y"]),
          legend.title = element_text(size = legend_size["title"]),
          legend.text=element_text(size = legend_size["elements"]))

  if (length(variables) == 1) {
    p.bs <- p.bs + geom_point(data = df.box, aes_string(fill = variables["fill"]), shape = 21,
                              size = markersize, alpha = alpha, color = "black",
                              position = position_jitter(width = jitter))

  } else if (length(variables) == 2) {
    p.bs <- p.bs + geom_point(data = df.box, aes_string(fill = variables["fill"], shape = variables["shape"]),
                              size = markersize, alpha = alpha, color = "black",
                              position = position_jitter(width = jitter)) +
      scale_shape_manual(name = variables["shape"], values = shapes)

  } else { return("Up to two variables allowed") }

  p.bs <- p.bs + scale_fill_manual(name = variables[1], values = colors,
                                   guide = guide_legend(override.aes = aes(shape = 21, size = 7)))

  return(p.bs)
}
