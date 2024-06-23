% Edit documentation in R/save_results.R
\name{save_results}
\alias{save_results}
\title{This function will save the results of a Differential Expression analysis}
\usage{
save_results(df, name, l2fc = 0)
}
\arguments{
\item{df}{A dataframe with the results of a Differential Expression analysis.}

\item{name}{The name to be used to save the tables, without file extension.}

\item{l2fc}{The cut-off of Log2(Fold Change) for the over- and under-expressed tables. Default = 0.}
}
\description{
This function takes as input the output of the function "results()" of DEseq2.
And will save 3 tables:
\itemize{
\item A table with all genes
\item A table including only the over-expressed genes
\item A table including only the under-expressed genes
}
}
