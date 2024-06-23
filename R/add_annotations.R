########################
# Add gene annotations #
########################

#' A function to add annotations to a table of gene counts.
#'
#' @param object A table of gene counts (rows: genes, columns:samples, rownames: ENSEMBL gene IDs).
#' @param reference A reference table with the annotations. Default: geneID.details.
#' @param variables Annotations (columns) from the reference table to add, i.e. gene symbol, gene_length, description.
#' @export

add_annotations <- function(object, reference = geneID.details, variables = NULL){

  df <- as.data.frame(object)
  df$geneID <- rownames(df)
  index <- match(df$geneID, geneID.details$geneID)

  # Add a new column for each variable
  for(i in 1:length(variables)){
    df[,variables[i]] <- geneID.details[index, variables[i]]
  }

  return(df)
}
