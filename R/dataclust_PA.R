####################
# Example datasets #
####################

#' Example gene set list for pathway analysis clustering
#'
#' A named list of 40 curated gene sets spanning four biological themes:
#' apoptosis & cell death, cell cycle & DNA damage, immune response &
#' inflammation, and metabolism. Gene set names follow standard database
#' conventions (`KEGG_`, `HALLMARK_`, `GO_`) and gene symbols are real human
#' genes. Designed to be used as input to [calc_jaccard()].
#'
#' @format A named list of 40 elements. Each element is a character vector of
#'   human gene symbols (HGNC) belonging to that gene set. Gene set sizes range
#'   from 11 to 20 genes.
#'
#' @source Curated manually for OmicsKit examples, based on KEGG, MSigDB
#'   Hallmark, and Gene Ontology gene set collections.
#'
#' @examples
#' data(geneset_list)
#'
#' # How many gene sets?
#' length(geneset_list)
#'
#' # Inspect one gene set
#' geneset_list[["KEGG_APOPTOSIS"]]
#'
#' # Use with calc_jaccard()
#' data(camera_results)
#' jac <- calc_jaccard(geneset_list, camera_results, fdr_th = 0.05)
#'
#' @seealso [calc_jaccard()], [camera_results]
"geneset_list"


#' Example CAMERA enrichment results for pathway analysis clustering
#'
#' A data frame simulating the output of a CAMERA differential expression
#' analysis, containing significance values for the 40 gene sets in
#' [geneset_list]. Approximately 60% of gene sets have FDR < 0.05, providing
#' enough significant sets for meaningful clustering. Designed to be used
#' alongside [geneset_list] as input to [calc_jaccard()].
#'
#' @format A data frame with 40 rows and 4 columns:
#'   \describe{
#'     \item{GeneSet}{Character. Gene set name, matching the names in
#'       [geneset_list].}
#'     \item{Direction}{Character. Enrichment direction: `"Up"` or `"Down"`.}
#'     \item{PValue}{Numeric. Raw p-value from the simulated CAMERA test.}
#'     \item{FDR}{Numeric. Benjamini-Hochberg adjusted p-value.}
#'   }
#'
#' @source Simulated with `set.seed(1905)` in `data-raw/example_PA.R` for
#'   OmicsKit examples.
#'
#' @examples
#' data(camera_results)
#'
#' # Overview
#' head(camera_results)
#'
#' # How many gene sets are significant?
#' sum(camera_results$FDR < 0.05)
#'
#' # Use with calc_jaccard()
#' data(geneset_list)
#' jac <- calc_jaccard(geneset_list, camera_results, fdr_th = 0.05)
#'
#' @seealso [calc_jaccard()], [geneset_list]
"camera_results"
