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
