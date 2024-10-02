% Edit documentation in R/nice_UMAP.R
\name{nice_UMAP}
\alias{nice_UMAP}
\title{Function to Make UMAP Plots}
\usage{
nice_UMAP(object,
	  neighbors = 4,
	  components = 3,
	  epochs = 10000,
	  returnData = FALSE,
          variables = c(fill = "VarFill", shape = "VarShape"),
          colors = NULL,
	  shapes = NULL,
	  size = 7,
	  alpha = 1,
          legend_names = c(fill = "Label Fill", shape = "Label Shape"),
          legend_title = 20,
	  legend_elements = 16,
	  legend_pos = c(0.80, 0.80, "right"),
          labels = NULL, # c(var = "patient", size = 2)
          name_tags = NULL) # c(var = "label", size = 3, minlen = 2, box = 0.5)
}
\arguments{
  \item{object}{
    A DESeq object already transformed with the variance stabilizing or rlog transformations.}

  \item{neighbors}{
    Number of nearest neighbors to consider. Default is 4.}

  \item{components}{
    Number of components to be considered for the dimensional reduction. Default is 3.}

  \item{epochs}{
    Number of iterations for the UMAP algorithm. Default is 10000.}

  \item{returnData}{
    Indicates if the function should return the data (TRUE) or the plot (FALSE). Default is FALSE.}

  \item{variables}{
    To indicate the variables to be used as Shape and Fill of the markers.}

  \item{legend_names}{
    The names to be used for the legend of the Shape and Fill.}

  \item{size}{
    Size of the marker. Default is 7.}

  \item{alpha}{
    Transparency of the marker, ranging from 0 (transparent) to 1 (no transparency). Default is 1.}

  \item{colors}{
    Vector of colors to be used for the categories of the variable assigned as Marker Fill.}

  \item{shapes}{
    Vector of shapes to be used for the categories of the variable assigned as Marker Shape.}

  \item{legend_title}{
    Font size of the legend title. Default is 20.}

  \item{legend_elements}{
    Font size of the elements of the legend. Default is 16.}

  \item{legend_pos}{
    Position of the legend in the plot. Default is c(0.80, 0.80, "right").}

  \item{labels}{
    A vector containing the variable to be used as labels (name inside the marker), and the label size. Example: c(var = "patient", size = 2). Default is NULL (no labels).}

  \item{name_tags}{
    A vector containing the variable to be used as name tags (name outside the marker), tag size, minimum distance in order to add an arrow connecting the tag and the marker, and minimum distance from the tag and the center of the marker. Example: c(var = "label", size = 3, minlen = 2, box = 0.5). Default is NULL (no name tags).}
}
\description{
Function to make UMAP plots.
}
