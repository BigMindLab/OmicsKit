########################
# Add gene annotations #
########################

#' A function to add annotations to a table of gene counts.
#'
#' @param object A table of gene counts (rows: genes, columns: samples, rownames are ENSEMBL gene IDs).
#' @param reference A reference table with the annotations including a column named "geneID".
#' @param variables Character vector of columns in `reference` to add. If NULL (default), all columns except geneID are used.
#' @param data_frame Logical; if TRUE, coerce `object` to a data.frame first. Default: FALSE.
#'
#' @return The input `object` as a data frame with additional columns from
#'   `reference` joined by Ensembl gene ID. A `geneID` column is added
#'   containing the row names of the original object.
#'
#' @examples
#' \dontrun{
#' data(norm_counts)
#'
#' # Requires a reference table with a "geneID" column.
#' # Use get_annotations() to generate it:
#' annotations <- get_annotations(
#'   ensembl_ids = rownames(norm_counts),
#'   mode        = "genes"
#' )
#'
#' # Add gene symbol and biotype columns to the counts matrix
#' norm_counts_annot <- add_annotations(
#'   object    = norm_counts,
#'   reference = annotations,
#'   variables = c("symbol", "biotype")
#' )
#'
#' # Inspect result
#' head(norm_counts_annot[, c("geneID", "symbol", "biotype")])
#'
#' # Add all annotation columns (variables = NULL uses everything)
#' norm_counts_full <- add_annotations(
#'   object    = norm_counts,
#'   reference = annotations
#' )
#' }
#'
#' @seealso [get_annotations()] to generate the `reference` table;
#'   [norm_counts] for an example input matrix.
#'
#' @export

add_annotations <- function(object, reference, variables = NULL, data_frame = FALSE){

  df <- if (data_frame) as.data.frame(object, stringsAsFactors = FALSE) else object
  df$geneID <- rownames(df)

  if (is.null(variables)) {
    variables <- setdiff(colnames(reference), "geneID")
  }

  index <- match(df$geneID, reference$geneID)
  df[, variables] <- reference[index, variables]

  return(df)
}
