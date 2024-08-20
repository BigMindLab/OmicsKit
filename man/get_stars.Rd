% Edit documentation in R/get_stars.R
\name{get_stars}
\alias{get_stars}
\title{Obtain annotations to display significance}
\usage{
get_stars(
  geneID,
  object,
  thresholds = c(0.001, 0.01, 0.1, 0.25)
)
}
\arguments{
  \item{geneID}{Ensembl ID of the gene of interest.}

  \item{object}{DESeq2 results object of a comparison.}

  \item{thresholds}{Vector with 4 values of significance. Default c(0.001, 0.01, 0.1, 0.25).}
}
\description{
  This function will create asteriscs (*) from DESeq2 results objects to represent the significance of comparisons.
}
