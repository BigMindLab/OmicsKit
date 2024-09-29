#########################
# Function save_results #
#########################

#' Save the results of a Differential Expression analysis
#'
#' This function takes as input the output of the function "results()" of DEseq2.
#' And will save 3 tables:
#' - A table with all genes
#' - A table including only the over-expressed genes
#' - A table including only the under-expressed genes
#'
#' @param df A dataframe with the results of a Differential Expression analysis.
#' @param name The name to be used to save the tables, without file extension.
#' @param l2fc The cut-off of Log2(Fold Change) for the over- and under-expressed tables. Default = 0.
#' @export

save_results <- function(df, name, l2fc = 0){

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop(
      "Package \"openxlsx\" must be installed to use this function.",
      call. = FALSE
      )
  }

  names(df)[names(df) == "padj"] <- "FDR"

  # Saving all genes:
  openxlsx::write.xlsx(df, colNames = T, rowNames = F, append = F,
                       file = paste0(name, "_full.xlsx"), overwrite = T)

  #Saving over-expressed genes:
  df.sig.fold_over <- subset(df, ((FDR < cutoff_alpha) & !is.na(FDR)) &
                               log2FoldChange >= l2fc)
  openxlsx::write.xlsx(df.sig.fold_over, colNames = T, rowNames = F, append = F,
                       file = paste0(name, "_Overexp.xlsx"), overwrite = T)

  #Saving under-expressed genes:
  df.sig.fold_under <- subset(df, ((FDR < cutoff_alpha) & !is.na(FDR)) &
                                log2FoldChange <= -l2fc)
  openxlsx::write.xlsx(df.sig.fold_under, colNames = T, rowNames = F, append = F,
                       file = paste0(name, "_Underexp.xlsx"), overwrite = T)
}
