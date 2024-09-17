% Edit documentation in R/nice_BSV.R
\name{nice_BSV}
\alias{nice_BSV}
\title{Function to make Box-Scatter-Violin plots.}
\usage{
nice_BSV(
  object = NULL,
  variables = c(fill = "VarFill", shape = "VarShape"),
  genename = NULL,
  symbol = NULL,
  labels = c("N", "P", "R", "M"),
  categories = c("normal", "primary", "recurrence", "metastasis"),
  colors = NULL,
  shapes = NULL,
  markersize = NULL,
  alpha = 0.8,
  width = NULL,
  height = NULL,
  jitter = 0.2,
  dpi = 150,
  save = FALSE,
  title_size = c(axis = 20, fig = 24),
  label_size = c(x = 20, y = 16),
  legend_size = c(title = 14, elements = 12)
)
}
\arguments{
\item{object}{A DEseq object already transformed with the variance stabilizing or rlog transformations.}

\item{variables}{To indicate the variables to be used as Shape and Fill of the markers.}

\item{genename}{The gene name to be used for the plot.}

\item{symbol}{The gene symbol to be used for the plot.}

\item{labels}{A vector containing the x-labels of the box-plot. Default: c("N", "P", "R", "M").}

\item{categories}{A vector containing the labels for the legend. Default: c("normal", "primary", "recurrence", "metastasis").}

\item{colors}{Vector of colors to be used for the categories of the variable assigned as Marker Fill.}

\item{shapes}{Vector of shapes to be used for the categories of the variable assigned as Marker Shape.}

\item{markersize}{Size of the marker.}

\item{alpha}{Transparency of the marker, which goes from 0 (transparent) to 1 (no transparent). Default: 0.8.}

\item{width}{Width of the plot.}

\item{height}{Height of the plot.}

\item{jitter}{Random deviation added to the dots. Default: 0.2.}

\item{dpi}{DPI of the plot. Default: 150.}

\item{save}{To save the plot. Default: FALSE.}

\item{title_size}{Font of the title and axis names. Default: c(axis = 20, fig = 24).}

\item{label_size}{Font of the labels (x-axis) and numbers (y-axis). Default: c(x = 20, y = 16).}

\item{legend_size}{Font of the title and elements of the legend. Default: c(title = 14, elements = 12).}
}
\description{
This function will make a Boxplot, using a DEseq object.
It will show the data points on top with a small deviation (jitter) for a better visualization.
}
