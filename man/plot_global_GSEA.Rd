% Generated by roxygen2: do not edit by hand
\name{plot_global_GSEA}
\alias{plot_global_GSEA}
\title{Plot global GSEA results}
\usage{
plot_global_GSEA(data, geneset_col, collection_col, nes_col,
                 logfdr_col, output_path_base, width_output,
                 height_output, text_size_genesets = 5,
                 text_size_collection = 5)
}
\arguments{
  \item{data}{Data frame containing the GSEA results.}
  \item{geneset_col}{Name of the column containing the genesets.}
  \item{collection_col}{Name of the column containing the MSigDB collections.}
  \item{nes_col}{Name of the column containing NES values.}
  \item{logfdr_col}{Name of the column containing \eqn{-\log_{10}(FDR)} values.}
  \item{output_path_base}{Path and base name for the output files (no extension).}
  \item{width_output}{Width of the output plot in inches.}
  \item{height_output}{Height of the output plot in inches.}
  \item{text_size_genesets}{Text size for the geneset labels.}
  \item{text_size_collection}{Text size for the collection labels.}
}
\description{
Generates a multi-panel plot summarizing GSEA results, showing pathways,
normalized enrichment scores (NES), and log-transformed FDR values
grouped by MSigDB collection.
}
\details{
This function creates a visual summary of GSEA results by combining:
\itemize{
  \item \strong{Pathway labels} on the left,
  \item \strong{NES bars} filled by \emph{log10 FDR} in the center,
  \item \strong{Collection labels} on the right,
  \item Vertical titles on both ends ("Pathways" and "MSigDB").
}
The final figure is saved as both PDF and JPG.
}
\examples{
\dontrun{
df <- readr::read_tsv("collections_merged_gsea_data.tsv")
plot_global_GSEA(df, "Geneset", "COLLECTION", "NES", "Log10FDR",
                 "output/plotGSEA", 12, 8)
}
}
