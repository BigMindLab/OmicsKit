########################
# Add gene annotations #
########################

#' A function to add annotations to a table of gene counts.
#'
#' @param object A table of gene counts (rows: genes, columns: samples, rownames are ENSEMBL gene IDs).
#' @param reference A reference table with the annotations including a column named "geneID".
#' @param variables Character vector of columns in `reference` to add. If NULL (default), all columns except geneID are used.
#' @param data_frame Logical; if TRUE, coerce `object` to a data.frame first. Default: FALSE.
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
