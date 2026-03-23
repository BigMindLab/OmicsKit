######################
# gsea_results data  #
######################

#' Simulated GSEA pathway analysis results for TCGA-LUAD
#'
#' A simulated data frame representing the output of [merge_PA()] for three
#' pairwise comparisons of TCGA-LUAD samples across 40 gene sets from three
#' MSigDB collections (HALLMARK, KEGG, GO). Gene sets and gene memberships are
#' derived from [geneset_list]. NES values and FDR are simulated with
#' `set.seed(174)` to produce realistic enrichment patterns, where ~60% of
#' gene sets per comparison are significant (FDR < 0.05).
#'
#' This dataset is designed to demonstrate [splot_PA()], [multiplot_PA()],
#' [getgenesPA()], [addgenesPA()], and [heatmap_PA()] without requiring
#' external GSEA output files.
#'
#' @format A data frame with 120 rows (40 gene sets x 3 comparisons) and 15
#'   columns:
#'   \describe{
#'     \item{NAME}{Character. Gene set name, matching the names in
#'       [geneset_list].}
#'     \item{SIZE}{Integer. Number of genes in the gene set.}
#'     \item{ES}{Numeric. Enrichment score.}
#'     \item{NES}{Numeric. Normalized enrichment score.}
#'     \item{NOM p-val}{Numeric. Nominal p-value.}
#'     \item{FDR}{Numeric. False discovery rate. Approximately 60% of gene
#'       sets per comparison have FDR < 0.05.}
#'     \item{FWER p-val}{Numeric. Family-wise error rate.}
#'     \item{RANK AT MAX}{Integer. Gene rank at maximum enrichment score.}
#'     \item{Log10FDR}{Numeric. `-log10(FDR)`.}
#'     \item{tags}{Numeric. Fraction of gene set in the leading edge (0-1).}
#'     \item{list}{Numeric. Fraction of the ranked list used (0-1).}
#'     \item{signal}{Numeric. Enrichment signal strength (0-1).}
#'     \item{LEADING EDGE}{Character. Leading edge string in GSEA format
#'       (e.g., `"tags=20%, list=35%, signal=15%"`).}
#'     \item{COLLECTION}{Character. MSigDB collection name: `"HALLMARK"`,
#'       `"KEGG"`, or `"GO"`.}
#'     \item{COMPARISON}{Character. Comparison name: `"TumorVsNormal"`,
#'       `"MetastasisVsNormal"`, or `"MetastasisVsTumor"`.}
#'   }
#'
#' @source Simulated with `set.seed(174)` in `data-raw/gsea_results.R`.
#'   Gene set names and memberships derived from [geneset_list]. NES values
#'   and significance are simulated to reflect realistic GSEA output patterns.
#'
#' @examples
#' data(gsea_results)
#'
#' # Overview
#' dim(gsea_results)
#' table(gsea_results$COMPARISON)
#' table(gsea_results$COLLECTION)
#'
#' # How many gene sets are significant per comparison?
#' tapply(gsea_results$FDR < 0.05, gsea_results$COMPARISON, sum)
#'
#' # Single comparison plot
#' single <- gsea_results[gsea_results$COMPARISON == "TumorVsNormal", ]
#' \dontrun{
#' splot_PA(
#'   data           = single,
#'   geneset_col    = "NAME",
#'   collection_col = "COLLECTION",
#'   nes_col        = "NES",
#'   fdr_col        = "FDR"
#' )
#' }
#'
#' # Multi-comparison plot
#' \dontrun{
#' multiplot_PA(
#'   data             = gsea_results,
#'   comparison_col   = "COMPARISON",
#'   facet_col        = "NAME",
#'   fdr_col          = "FDR",
#'   comparison_order = c("TumorVsNormal", "MetastasisVsNormal",
#'                        "MetastasisVsTumor")
#' )
#' }
#'
#' @seealso [merge_PA()] which produces this format from real GSEA output;
#'   [splot_PA()], [multiplot_PA()] for visualization;
#'   [getgenesPA()], [addgenesPA()] for gene-level annotation;
#'   [geneset_list] for the gene set memberships used here.
"gsea_results"
