################
# Function tmp #
################

#' Function to calculate the TPMs from a table of raw gene counts.
#'
#' TPM: Transcript per million. See https://www.biostars.org/p/273537/
#' The input table is numeric:
#' - The row names are the gene identifiers (ensembl ID).
#' - The column names represent the samples.
#' The gene lengths are in a column of a dataframe with the same row order.
#'
#' @param raw_counts A table with the gene counts.
#' @param gene_lengths A column with the gene lengths.
#' @export

tpm <- function(raw_counts, gene_lengths) {

  x <- raw_counts*1e3 / gene_lengths
  return(t(t(x)*1e6 / colSums(x)))

}
