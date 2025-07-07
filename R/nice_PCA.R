#####################
# Function nice_PCA #
#####################

#' Function to make nice PCA plots
#'
#' This was inspired on the plotPCA function from DESeq2, made by Wolfgang Huber
#' But including some improvements made by David Requena. Now it allows:
#' * To choose which PCs to plot.
#' * To use one or two features to represent as the fill or shape of the markers.
#' * To provide the colors, shapes and fonts.
#'
#' @param object A matrix of counts with genes as rows and sample ids as columns.
#' @param annotations A data frame of annotation, including sample ids and variables to plot. Default: NULL.
#' @param PCs A vector indicating the two Principal Components to plot. Default: c(1,2).
#' @param ntop Number of top genes to use for principal components, selected by highest row variance. Default: NULL.
#' @param variables To indicate the variables to be used as Shape and Fill of the markers.
#' @param legend_names The names to be used for the legend of the Shape and Fill.
#' @param size Size of the marker. Default: 5.
#' @param alpha Transparency of the marker, which goes from 0 (transparent) to 1 (no transparent). Default: 1.
#' @param colors Vector of colors to be used for the categories of the variable assigned as Marker Fill.
#' @param shapes Vector of shapes to be used for the categories of the variable assigned as Marker Shape.
#' @param title Plot title. Default: NULL.
#' @param legend_title Font of the legend title. Default: 16.
#' @param legend_elements Font of the elements of the legend Default: 14.
#' @param legend_pos Position of the legend inside the plot. Example: c(0.80, 0.80). Default: NULL.
#' @param returnData Indicates if the function should return the data (TRUE) or the plot (FALSE). Default: FALSE.
#' @param labels A vector containing the variable to be used as labels (name inside the marker), and the label size. Example: c(var = "patient", size = 2). Default: NULL (no labels).
#' @param name_tags A vector containing the variable to be used as name tags (name outside the marker), tag size, minimum distance in order to add an arrow connecting the tag and the marker, and minimum distance from the tag and the center of the marker. Example: c(var = "label", size = 3, minlen = 2, box = 0.5). Default: NULL (no name tags).
#' @param cluster_data Indicates if the function generates the clusters (TRUE) or not (FALSE). This new cluster variable can be used as fill or shape. Default: FALSE.
#' @param n_clusters Number of cluster categories. Default: 3.
#' @param transform Logical. Indicates whether to log2 transform the input `object` or not. Default: FALSE.
#' @param outPCs Number of Principal Components to keep if `returnData` is TRUE. Default: 50.
#' @param returnData Indicates if the function should return the data (TRUE) or the plot (FALSE). Default: FALSE.
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

nice_PCA <- function(object, annotations = NULL, PCs = c(1,2), ntop = NULL,
                     variables = c(fill = "VarFill", shape = "VarShape"),
                     legend_names = c(fill = "Sample Type", shape = "Library"),
                     size = 5, alpha = 1, colors = NULL, shapes = NULL, title = NULL,
                     legend_title = 16, legend_elements = 14, legend_pos = NULL,
                     labels = NULL, name_tags = NULL, cluster_data = FALSE,
                     n_clusters = 3, transform = FALSE, outPCs = 50, returnData = FALSE)

{
  expr <- if (transform) log2(object + 0.001) else object

  if (!is.null(ntop)) {
    # Estimate and select the top variances
    top.variances <- order(matrixStats::rowVars(expr), decreasing = TRUE)[1:min(ntop, nrow(expr))]

    # Principal Component Analysis
    pca <- stats::prcomp(t(expr[top.variances, , drop = FALSE]), scale = TRUE)

  } else {
    pca <- stats::prcomp(t(expr), scale = TRUE)
  }

  # Calculate the percent of variance per component
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)

  if (returnData) {
    n_take  <- min(outPCs, ncol(pca$x))
    pca_data <- pca$x[, seq_len(n_take), drop = FALSE]

    colnames(pca_data) <- paste0("PC", seq_len(n_take))
    attr(pca_data, "percentVar") <- percentVar[seq_len(n_take)]

    return(pca_data)

  }

  if (is.null(annotations)) {
    stop("`annotations` must be supplied when returnData = FALSE", call. = FALSE)
  }

  # Extract data for clustering if `cluster_data` is TRUE
  if (cluster_data) {
    if (!requireNamespace("stats", quietly = TRUE)) {
      stop(
        "Package \"stats\" must be installed to perform clustering.",
        call. = FALSE
      )
    }

    expr_clust <- if (!is.null(ntop)) expr[top.variances, , drop = FALSE] else expr
    hc <- stats::hclust(stats::dist(t(expr_clust)), method = "ward.D2")

    # Add cluster labels to annotations
    annotations$cluster <- factor(stats::cutree(hc, k = n_clusters))
  }

  df.pca <- data.frame(PCx = pca$x[, PCs[1]], PCy = pca$x[, PCs[2]]) %>%
    tibble::rownames_to_column(var = "id") %>%
    dplyr::inner_join(annotations, by = "id")

  if (length(variables) == 3) {

    if (!requireNamespace("ggnewscale", quietly = TRUE)) {
      stop(
        "Package \"ggnewscale\" must be installed for third variable aesthetics.",
        call. = FALSE
      )
    }

    p.nicePCA <- ggplot(data = df.pca, aes(x = .data[["PCx"]], y = .data[["PCy"]], fill = .data[[variables[3]]])) +
      stat_ellipse(geom = "polygon", alpha = 0.2) + scale_fill_manual(name = legend_names[3], values = c("orange", "thistle4", "yellow")) +

      ggnewscale::new_scale_fill() +

      geom_point(aes(fill = .data[[variables[1]]], shape = .data[[variables[2]]]), size = size, alpha = alpha) +
      scale_shape_manual(name = legend_names[2], values = shapes,
                         guide  = guide_legend(override.aes = list(size = 7), keyheight = 1.7))

  } else if (length(variables) == 2) {
    p.nicePCA <- ggplot(data = df.pca, aes(x = .data[["PCx"]], y = .data[["PCy"]], fill = .data[[variables[1]]], shape = .data[[variables[2]]])) +
      geom_point(size = size, alpha = alpha) + labs(fill = legend_names[1], shape = legend_names[2]) +
      scale_shape_manual(values = shapes, guide = guide_legend(override.aes = list(size = 7), keyheight = 1.7))

    if (is.numeric(df.pca[[variables[1]]])) {
      p.nicePCA <- p.nicePCA + scale_fill_gradient2(low = "blue", mid = "white", high = "red")
    } else {
      p.nicePCA <- p.nicePCA + scale_fill_manual(values = colors, guide = guide_legend(override.aes = aes(shape = 21, size = 7)))
    }

  } else if (length(variables) == 1) {
    p.nicePCA <- ggplot(data = df.pca, aes(x = .data[["PCx"]], y = .data[["PCy"]], fill = .data[[variables[1]]])) +
      geom_point(size = size, alpha = alpha, shape = 21) + labs(fill = legend_names)

    if (is.numeric(df.pca[[variables]])) {
      p.nicePCA <- p.nicePCA + scale_fill_gradient2(low = "blue", mid = "white", high = "red")
    } else {
      p.nicePCA <- p.nicePCA + scale_fill_manual(values = colors, guide = guide_legend(override.aes = aes(size = 7)))
    }
  }

  p.nicePCA <- p.nicePCA + coord_fixed() + theme_bw()

  if (!is.null(title)) {
    p.nicePCA <- p.nicePCA + labs(title = title) + theme(plot.title = element_text(size = 18, hjust = 0.5))
  }

  p.nicePCA <- p.nicePCA +
    xlab(paste0("PC", PCs[1], ": ", round(percentVar[PCs[1]] * 100), "% variance")) +
    ylab(paste0("PC", PCs[2], ": ", round(percentVar[PCs[2]] * 100), "% variance")) +
    theme(axis.title = element_text(size = legend_title),
          axis.text = element_text(size = legend_elements),
          legend.title = element_text(size = legend_title),
          legend.text = element_text(size = legend_elements))

  if (is.null(legend_pos) == FALSE) {
    p.nicePCA <- p.nicePCA + theme(legend.background = element_rect(color = "black"), legend.position = "inside",
                                   legend.position.inside = c(legend_pos[1], legend_pos[2]))
  }

  if (is.null(labels) == FALSE) {
    if (length(labels) == 1) {
      lab_size <- as.numeric(labels)
      df.pca$label <- seq_len(nrow(df.pca))
    } else {
      lab_var <- labels[1]
      if (!lab_var %in% names(df.pca)) {
        stop("`labels[1]` must be a column in `annotations`")
      }
      lab_size <- labels[2]
      df.pca$label <- df.pca[[lab_var]]
    }

    # Add the labels to the plot
    p.nicePCA <- p.nicePCA +
      geom_text(data = df.pca, aes(label = .data[["label"]]), color = "black", size = lab_size)
  }

  if (is.null(name_tags) == FALSE) {

    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      stop(
        "Package \"ggrepel\" must be installed to use this function.",
        call. = FALSE
        )
    }

    if (length(name_tags) == 1) {
      tag_var <- "id"
      tag_size <- as.numeric(name_tags)
      minlen <- 2
      boxpad <- 0.5
    } else {
      tag_var <- name_tags[1]
      if (!tag_var %in% names(df.pca)) {
        stop("`name_tags[1]` must be a column in `annotations`")
      }
      tag_size <- as.numeric(name_tags[2])
      minlen <- as.numeric(name_tags[3])
      boxpad <- as.numeric(name_tags[4])
    }

    # Add the column of name tags to the data frame
    df.pca$tag <- df.pca[[tag_var]]

    # Add the name tags to the plot
    p.nicePCA <- p.nicePCA +
      ggrepel::geom_text_repel(data = df.pca, aes(label = .data[["tag"]]), color = "black", size = tag_size,
                               min.segment.length = unit(minlen, "lines"), box.padding = unit(boxpad, "lines"))
  }

  return(p.nicePCA)
}
