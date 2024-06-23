#####################
# Function nice_PCA #
#####################

#' Function to make nice PCA plots
#'
#' This was inspired on the plotPCA function from DESeq2, made by Wolfgang Huber
#' But including some improvements made by David Requena. Now it allows:
#' - To choose which PCs to plot.
#' - To use one or two features to represent as the fill or shape of the markers.
#' - To provide the colors, shapes and fonts.
#'
#' @param object A DEseq object already transformed with the variance stabilizing or rlog transformations.
#' @param PCs A vector indicating the two Principal Components to plot. Default: c(1,2).
#' @param ntop Number of top genes to use for principal components, selected by highest row variance.
#' @param returnData Indicates if the function should return the data (TRUE) or the plot (FALSE). Default: FALSE.
#' @param variables To indicate the variables to be used as Shape and Fill of the markers.
#' @param legend_names The names to be used for the legend of the Shape and Fill.
#' @param size Size of the marker. Default: 7.
#' @param alpha Transparency of the marker, which goes from 0 (transparent) to 1 (no transparent). Default: 1.
#' @param colors Vector of colors to be used for the categories of the variable assigned as Marker Fill.
#' @param shapes Vector of shapes to be used for the categories of the variable assigned as Marker Shape.
#' @param legend_title Font of the legend title. Default: 20.
#' @param legend_elements Font of the elements of the legend Default: 16.
#' @param legend_pos Position of the legend in the plot. Default: c(0.80, 0.80, "right").
#' @param labels A vector containing the variable to be used as labels (name inside the marker), and the label size. Example: c(var = "patient", size = 2). Default: NULL (no labels).
#' @param name_tags A vector containing the variable to be used as name tags (name outside the marker), tag size, minimum distance in order to add an arrow connecting the tag and the marker, and minimum distance from the tag and the center of the marker. Example: c(var = "label", size = 2, minlen = 2, box = 0.6). Default: NULL (no name tags).
#' @export

nice_PCA <- function(object, PCs = c(1,2), ntop = 200, returnData = FALSE,
                     variables = c(fill = "VarFill", shape = "VarShape"),
                     legend_names = c(fill = "Sample Type", shape = "Library"),
                     size = 7, alpha = 1, colors = NULL, shapes = NULL,
                     legend_title = 20, legend_elements = 16, legend_pos = c(0.80, 0.80, "right"),
                     labels = NULL, # c(var = "patient", size = 2)
                     name_tags = NULL) # c(var = "label", size = 2, minlen = 2, box = 0.6)

{
  require("matrixStats")
  require("ggplot2")

  # Estimate the variance in each row (gene or transcript)
  variances <- rowVars(assay(object))

  # Select the top variances (value passed as ntop)
  top.variances <- order(variances, decreasing = TRUE)[1:min(ntop, length(variances))]

  # Principal Component Analysis
  pca <- prcomp(t(assay(object)[top.variances, ]))

  # Calculate the percent of variance per component
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)

  if (!all(variables %in% names(colData(object)))) {
    stop("The variable names should be columns of the DESeq2 dataset, examine: colData(object)")
  }

  # Prepare a data frame with the PCs to plot and the metadata variable selected
  d <- data.frame(PCx = pca$x[, PCs[1]], PCy = pca$x[, PCs[2]],
                  colData(object)[, variables, drop = FALSE])

  if (length(variables) == 2) {

    p.nicePCA <- ggplot(data = d, aes_string(x = "PCx", y = "PCy", fill = variables[1], shape = variables[2])) +
      geom_point(size = size, alpha = alpha) + scale_shape_manual(values = shapes) +
      scale_fill_manual(values = colors, guide = guide_legend(override.aes = aes(shape = 21))) +
      labs(fill = legend_names[1], shape = legend_names[2])

  } else if (length(variables) == 1) {

    p.nicePCA <- ggplot(data = d, aes_string(x = "PCx", y = "PCy", fill = variables)) +
      geom_point(size = size, alpha = alpha, shape = 21) +
      scale_fill_manual(values = colors, guide = guide_legend(override.aes = aes(shape = 21))) +
      labs(fill = legend_names)
  }

  p.nicePCA <- p.nicePCA + coord_fixed() + theme_bw() +
    theme(legend.title = element_text(size = legend_title),
          legend.text = element_text(size = legend_elements),
          legend.background = element_rect(color = "black"),
          legend.position = c(legend_pos[1], legend_pos[2]), legend.box.just = legend_pos[3]) +
    xlab(paste0("PC", PCs[1], ": ", round(percentVar[PCs[1]] * 100), "% variance")) +
    ylab(paste0("PC", PCs[2], ": ", round(percentVar[PCs[2]] * 100), "% variance"))

  if (is.null(labels) == FALSE) {

    # Add the column of labels to the data frame
    d <- data.frame(d, colData(object)[,c(labels[1]), drop = FALSE])

    # Add the labels to the plot
    p.nicePCA <- p.nicePCA +
      geom_text(aes(label = d[,labels[1]]),
                color = "black", size = as.numeric(labels[2]))
  }

  if (is.null(name_tags) == FALSE) {

    require("ggrepel")

    # Add the column of name tags to the data frame
    d <- data.frame(d, colData(object)[,c(name_tags[1]), drop = FALSE])

    # Add the name tags to the plot
    p.nicePCA <- p.nicePCA +
      geom_text_repel(aes(label = d[,name_tags[1]]),
                      color = "black", cex = as.numeric(name_tags[2]),
                      min.segment.length = unit(as.numeric(name_tags[3]), "lines"),
                      box.padding = unit(as.numeric(name_tags[4]), "lines"))
  }

  if (returnData) {
    attr(d, "percentVar") <- percentVar[PCs[1]:PCs[2]]
    return(d)
    # return(pca$x)
    # return(pca$x[, c(PCs[1], PCs[2])])
  }

  return(p.nicePCA)
}
