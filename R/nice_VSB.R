#####################
# Function nice_VSB #
#####################

#' Function to make Violin-Scatter-Box plots from data frames.
#'
#' This function will make a Boxplot, using a DEseq object.
#' It will show the data points on top with a small deviation (jitter) for a better visualization.
#'
#' @param object A data frame object with normalized counts genes(in rows) across samples(in columns).
#' @param annotations Data frame with annotations.
#' @param variables To indicate the variables to be used as Shape and Fill of the markers.
#' @param genename The gene name to be used for the plot.
#' @param symbol The gene symbol to display in the plot title. To obtain
#'  gene symbols from Ensembl IDs, use [get_annotations()].
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
#' @param box_width Width of the boxplot relative to the category spacing. Default: 0.25.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A ggplot2 object.
#'
#' @examples
#' data(brca_rna_expr_tumor_normal_filtered)
#' data(brca_rna_metadata_tumor_normal)
#'
#' nice_VSB(
#'   object      = brca_rna_expr_tumor_normal_filtered,
#'   annotations = brca_rna_metadata_tumor_normal,
#'   variables   = c(fill = "sample_type"),
#'   genename    = rownames(brca_rna_expr_tumor_normal_filtered)[1],
#'   categories  = c("normal", "tumor"),
#'   labels      = c("Normal", "Tumor"),
#'   colors      = c("steelblue", "firebrick"),
#'   shapes      = 21,
#'   markersize  = 3
#' )
#'
#' @seealso [nice_Volcano()] for genome-wide visualization; [detectability_filter()]
#'   to identify reliably expressed genes; [get_stars()] to add significance
#'   annotations; [brca_rna_expr_tumor_normal_filtered] for an example input matrix.
#'
#' @export

nice_VSB <- function (object = NULL, annotations, variables = c(fill = "VarFill", shape = "VarShape"),
                      genename = NULL, symbol = NULL, labels = c("N", "P", "R", "M"),
                      categories = c("normal", "primary", "recurrence", "metastasis"),
                      colors = NULL, shapes = NULL, markersize = NULL, alpha = 0.8, jitter = 0.2,
                      title_size = c(axis = 20, fig = 24), label_size = c(x = 20, y = 16),
                      legend_size = c(title = 14, elements = 12), box_width = 0.25)
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

  # Plot: layers are drawn back-to-front in the order they are added, so the
  # violin (background density) goes first, the jittered points (raw data) on
  # top of it, and the boxplot (summary statistics) last so it reads cleanly
  # on top of the point cloud. The boxplot uses a light neutral fill (#E5E5E5)
  # so the summary shape is legible without competing with the fill-colored
  # data points, and has no outlier points of its own, since individual points
  # are already shown by geom_point().
  p.bs <- ggplot(df.box, aes(x = .data$sample_type, y = log2_gc)) + theme_bw() +
    geom_violin(alpha = 0.1, scale = "width", fill = "yellow", color = "peru",
		show.legend = FALSE, trim = TRUE)

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

  p.bs <- p.bs +
    geom_boxplot(width = box_width, fill = "#E5E5E5", alpha = 0.8, outlier.shape = NA, color = "black") +
    labs(title = paste("Gene:", genename, symbol),
         x = expression("Sample Type"),
         y = expression("log"[2]*"(Normalized Gene Counts)")) +
    theme(plot.title = element_text(size = title_size["fig"]),
          axis.title = element_text(size = title_size["axis"]),
          axis.text.x = element_text(size = label_size["x"]),
          axis.text.y = element_text(size = label_size["y"]),
          legend.title = element_text(size = legend_size["title"]),
          legend.text=element_text(size = legend_size["elements"]))

  # Clean legend title: strip underscores from the raw column name (e.g.
  # "ER_group" -> "ER group") instead of showing the variable name verbatim.
  legend_label <- gsub("_", " ", variables[1])
  p.bs <- p.bs + scale_fill_manual(name = legend_label, values = colors,
                                   guide = guide_legend(override.aes = aes(shape = 21, size = 7)))

  return(p.bs)
}
