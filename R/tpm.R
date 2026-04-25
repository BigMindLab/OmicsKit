################
# Function tpm #
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
#'
#' @return A numeric matrix of the same dimensions as `raw_counts` with TPM
#'   values. Column sums equal 1,000,000 by definition.
#'
#' @note TPM normalizes for both sequencing depth and gene length, making
#'   values comparable between genes within a sample. It is not appropriate
#'   for differential expression analysis, use DESeq2 normalized counts
#'   ([norm_counts]) for that purpose. Gene lengths from [get_annotations()]
#'   are genomic lengths (including introns); for higher accuracy use
#'   transcript-level lengths.
#'
#' @examples
#' \dontrun{
#' data(raw_counts)
#' data(deseq2_results)
#'
#' # Gene lengths are needed, retrieve from get_annotations() or use
#' # pre-fetched lengths. Here we use the gene_length column if available.
#' annotations <- get_annotations(
#'   ensembl_ids = rownames(raw_counts),
#'   mode        = "genes"
#' )
#'
#' # Match gene lengths to raw_counts row order
#' gene_lengths <- annotations$gene_length[
#'   match(rownames(raw_counts), annotations$geneID)
#' ]
#'
#' # Calculate TPM
#' tpm_matrix <- tpm(raw_counts, gene_lengths)
#'
#' # Check: column sums should all be 1,000,000
#' round(colSums(tpm_matrix)[1:3])
#' }
#'
#' @seealso [get_annotations()] to obtain `gene_lengths`;
#'   [norm_counts] for DESeq2 size-factor normalized counts.
#'
#' @export

tpm <- function(raw_counts, gene_lengths) {

  x <- raw_counts*1e3 / gene_lengths
  return(t(t(x)*1e6 / colSums(x)))

}
