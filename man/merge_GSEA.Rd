% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_GSEA.R
\name{merge_GSEA}
\alias{merge_GSEA}
\title{Merge GSEA results data frames.}
\usage{
merge_GSEA(input_directory, output_file)
}
\arguments{
\item{input_directory}{The directory containing the GSEA collection results in TSV format.}

\item{output_file}{The output file to save the merged data. If not provided, the file will be saved in the input directory.}
}
\description{
After run GSEA_all.sh from GSEA.sh, move the .tsv files to a single directory
}
