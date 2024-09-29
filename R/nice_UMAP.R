######################
# Function nice_UMAP #
######################

#' Function to make UMAP plots.
#'
#' @param object A DEseq object already transformed with the variance stabilizing or rlog transformations.
#' @param neighbors Number of nearest neighbors to consider. Default: 4.
#' @param components Number of components to be consider for the dimensional reduction. Default: 3.
#' @param epochs Number of iterations. Default: 20000.
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
#' @param name_tags A vector containing the variable to be used as name tags (name outside the marker), tag size, minimum distance in order to add an arrow connecting the tag and the marker, and minimum distance from the tag and the center of the marker. Example: c(var = "label", size = 3, minlen = 2, box = 0.5). Default: NULL (no name tags).
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom SummarizedExperiment assay colData
#' @export

nice_UMAP <- function(object, neighbors = 4, components = 3, epochs = 20000, returnData = FALSE,
                      variables = c(fill = "VarFill", shape = "VarShape"),
                      colors = NULL, shapes = NULL, size = 7, alpha = 1,
                      legend_names = c(fill = "Label Fill", shape = "Label Shape"),
                      legend_title = 20, legend_elements = 16, legend_pos = c(0.80, 0.80, "right"),
                      labels = NULL, # c(var = "patient", size = 2)
                      name_tags = NULL) # c(var = "label", size = 3, minlen = 2, box = 0.5)

{

  if (!requireNamespace("umap", quietly = TRUE)) {
    stop(
      "Package \"umap\" must be installed to use this function.",
      call. = FALSE
      )
  }

  # Initialize umap
  umap.params = umap::umap.defaults
  umap.params$n_neighbors=neighbors
  umap.params$n_components=components
  umap.params$n_epochs=epochs
  umap.params$random_state=1
  umap.params$transform_state=1
  umap.params$verbose=TRUE

  # Create data frame
  umap_data <- umap::umap(t(assay(object)))
  df.umap <- data.frame(umap_data$layout) %>%
    tibble::rownames_to_column("id") %>%
    dplyr::inner_join(colData(object), by = "id")

  # Create plot
  if (length(variables) == 2) {

    p.umap <- ggplot(data = df.umap, aes_string(x = "X1", y = "X2", fill = variables[1], shape = variables[2])) +
      geom_point(size = size, alpha = alpha) + labs(fill = legend_names[1], shape = legend_names[2]) +
      scale_fill_manual(values = colors, guide = guide_legend(override.aes = aes(shape = 21, size = 9))) +
      scale_shape_manual(values = shapes, guide = guide_legend(override.aes = list(size = 7), keyheight = 1.7))

  } else if (length(variables) == 1) {

    p.umap <- ggplot(data = df.umap, aes_string(x = "X1", y = "X2", fill = variables[1])) +
      geom_point(size = size, alpha = alpha, shape = 21) +
      scale_fill_manual(values = colors, guide = guide_legend(override.aes = aes(shape = 21, size = 9)))

  }

  p.umap <- p.umap + coord_fixed() + theme_bw() +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
          legend.title = element_text(size=legend_title),
          legend.text=element_text(size=legend_elements),
          legend.background = element_rect(color = "black"),
          legend.box.just = legend_pos[3], legend.position = c(legend_pos[1], legend_pos[2]))

  # Manage labels and tags
  if (is.null(labels) == FALSE) {

    # Add the column of labels to the data frame
    df.umap <- data.frame(df.umap, colData(object)[, labels[1], drop = FALSE])

    # Add the labels to the plot
    p.umap <- p.umap +
      geom_text(aes(label = df.umap[,labels[1]]), color = "black", size = as.numeric(labels[2]))
  }

  if(is.null(name_tags) == FALSE) {

    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      stop(
        "Package \"ggrepel\" must be installed to use this function.",
        call. = FALSE
        )
    }

    # Add the column of name tags to the data frame
    df.umap <- data.frame(df.umap, colData(object)[, name_tags[1], drop = FALSE])

    # Add the name tags to the plot
    p.umap <- p.umap +
      ggrepel::geom_text_repel(aes(label = df.umap[,name_tags[1]]),
                               color = "black", cex = as.numeric(name_tags[2]),
                               min.segment.length = unit(as.numeric(name_tags[3]), "lines"),
                               box.padding = unit(as.numeric(name_tags[4]), "lines"))
  }

  if (returnData) {
    return(df.umap)
  } else {
    return(p.umap)
  }
}
