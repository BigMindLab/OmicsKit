######################
# Function get_stars #
######################

#' Obtain annotations to display significance.
#'
#' This function will create asteriscs (*) from DESeq2 results objects to represent the significance of comparisons.
#'
#' @param geneID Ensembl ID of the gene of interest.
#' @param object DESeq2 results object of a comparison.
#' @param thresholds Vector with 4 values of significance. Default c(0.001, 0.01, 0.1, 0.25).
#'
#' @return A single character string: `"****"`, `"***"`, `"**"`, `"*"`,
#'   `"ns"` (not significant), or `"Gene ID not found"` if the gene is absent
#'   from `object`.
#'
#' @examples
#' data(deseq2_results)
#'
#' # get_stars expects a column named "ensembl"
#' res <- deseq2_results
#' colnames(res)[colnames(res) == "gene_id"] <- "ensembl"
#'
#' # Get significance stars for the most significant gene
#' get_stars(
#'   geneID = res$ensembl[1],
#'   object = res
#' )
#'
#' # Custom thresholds
#' get_stars(
#'   geneID     = res$ensembl[1],
#'   object     = res,
#'   thresholds = c(0.001, 0.01, 0.05, 0.10)
#' )
#'
#' # Non-significant gene
#' get_stars(
#'   geneID = res$ensembl[nrow(res)],
#'   object = res
#' )
#'
#' @seealso [detect_filter()] to identify detectable genes before annotating;
#'   [nice_VSB()] where significance stars can be added to plots;
#'   [deseq2_results] for an example input.
#'
#' @export

get_stars <- function(geneID, object, thresholds = c(0.001, 0.01, 0.1, 0.25))
{
	# Ensure the geneID column is of the correct type of matching
	object <- as.data.frame(object)
	object$ensembl <- as.character(object$ensembl)

	# Find the row with the marching geneID and get de padj value
	padj <- object$padj[object$ensembl == geneID]

	# Check if gene ID was found
	if (length(padj) == 0) {
		return("Gene ID not found")
	}

	# Order thresholds vector increasingly
	thresholds <- sort(thresholds, decreasing = FALSE)

	# Retrieve stars
	if (is.na(padj) || padj > thresholds[4]) {
		stars <- "ns" # Not significant
	} else if (padj <= thresholds[1]) {
		stars <- "****"
	} else if (padj <= thresholds[2]) {
		stars <- "***"
	} else if (padj <= thresholds[3]) {
		stars <- "**"
	} else {
		stars <- "*"
	}

	return(stars)
}
