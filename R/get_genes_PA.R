######################
# Function getgenesPA #
######################

#' Extract gene members from pathway analysis results
#'
#' For each gene set in a pathway analysis results table, retrieves either the
#' leading edge genes (the subset that drives enrichment, in GSEA output only),
#' all genes in the gene set, or both. Leading edge genes are ordered by their
#' rank in the provided ranked gene list.
#'
#' This function is the recommended first step before [heatmap_PA()],
#' [addgenesPA()], and gene-level visualization with [nice_VSB()].
#'
#' @param pa_data A data frame of pathway analysis results. Must contain:
#'   * `NAME`: gene set name.
#'   * `SIZE`: number of genes in the gene set (GSEA output only).
#'   * `LEADING EDGE`: GSEA leading edge string
#'     (e.g., `"tags=20%, list=30%, signal=15%"`) (GSEA output only).
#'   Typically the output of [merge_PA()].
#' @param geneset_list A named list of gene sets, where each element is a
#'   character vector of gene symbols. Typically the output of [list_gmts()].
#' @param ranked_genes A character vector of gene symbols ordered by their
#'   ranking metric (e.g., stat, log2FC or signal-to-noise ratio), from most
#'   positive to most negative, not-significant in the middle of the list.
#'   Used to order leading edge genes by rank.
#' @param genes Character vector specifying which gene sets to return. One or
#'   both of `"le"` (leading edge genes only) and `"all"` (all genes in the
#'   gene set). Default: `c("all", "le")`.
#'
#' @return Depends on `genes`:
#'   * `genes = "le"` : A named list where each element is a character vector
#'     of leading edge gene symbols ordered by rank. The list has an attribute
#'     `type = "le"` used by [addgenesPA()].
#'   * `genes = "all"` : A named list where each element is a character vector
#'     of all gene symbols in that gene set. The list has an attribute
#'     `type = "all"` used by [addgenesPA()].
#'   * `genes = c("all", "le")` : A named list with two elements:
#'     * `$all`: named list of all gene set members.
#'     * `$le`: named list of leading edge genes ordered by rank.
#'
#' @examples
#' \dontrun{
#' gsl     <- list_gmts("path/to/gmt_folder/")
#' pa_data <- read.delim("gsea_results_merged.tsv")
#' ranked  <- deseq2_results$gene_id[order(deseq2_results$stat,
#'                                         decreasing = TRUE)]
#'
#' # Get both leading edge and all genes
#' gene_lists <- getgenesPA(pa_data, gsl, ranked, genes = c("all", "le"))
#' pa_annot   <- addgenesPA(pa_data, gene_lists)
#' head(pa_annot[, c("NAME", "all_genes", "le_genes")])
#'
#' # Get only leading edge genes
#' le_only  <- getgenesPA(pa_data, gsl, ranked, genes = "le")
#' pa_annot <- addgenesPA(pa_data, le_only)
#' head(pa_annot[, c("NAME", "le_genes")])
#'
#' # Access genes for a specific gene set
#' gene_lists$le[["KEGG_APOPTOSIS"]]
#' gene_lists$all[["KEGG_APOPTOSIS"]]
#' }
#'
#' @seealso [addgenesPA()] to append gene columns to pa_data;
#'   [heatmap_PA()] for heatmap visualization;
#'   [list_gmts()] to generate `geneset_list`;
#'   [merge_PA()] to generate `pa_data`.
#'
#' @export

getgenesPA <- function(pa_data, geneset_list, ranked_genes,
                       genes = c("all", "le")) {

  genes <- match.arg(genes, choices = c("all", "le"), several.ok = TRUE)

  if (!is.data.frame(pa_data)) {
    stop("`pa_data` must be a data frame.", call. = FALSE)
  }
  if (!"NAME" %in% colnames(pa_data)) {
    stop("`pa_data` must contain a column named 'NAME'.", call. = FALSE)
  }
  if (!is.list(geneset_list) || is.null(names(geneset_list))) {
    stop("`geneset_list` must be a named list. Use list_gmts() to generate it.",
         call. = FALSE)
  }
  if (!is.character(ranked_genes) || length(ranked_genes) == 0) {
    stop("`ranked_genes` must be a non-empty character vector.", call. = FALSE)
  }

  common_sets <- intersect(pa_data$NAME, names(geneset_list))
  if (length(common_sets) == 0) {
    stop(
      "No gene set names in `pa_data$NAME` match `geneset_list`. ",
      "Check that both use the same naming convention.",
      call. = FALSE
    )
  }

  # Internal helper: compute leading edge size from LEADING EDGE string
  .get_le_size <- function(size, tags_str) {
    tag_pct <- as.numeric(sub("%", "", sub("tags=", "", tags_str))) / 100
    round(size * tag_pct)
  }

  result <- list()

  if ("all" %in% genes) {
    result$all <- stats::setNames(
      lapply(common_sets, function(gs) geneset_list[[gs]]),
      common_sets
    )
  }

  if ("le" %in% genes) {
    le_list <- lapply(common_sets, function(gs) {
      row      <- pa_data[pa_data$NAME == gs, ][1, ]
      all_g    <- geneset_list[[gs]]
      tags_str <- strsplit(as.character(row[["LEADING EDGE"]]), ",")[[1]][1]
      le_size  <- .get_le_size(row$SIZE, tags_str)
      if (is.na(le_size) || le_size <= 0) return(character(0))
      gene_ranks <- match(all_g, ranked_genes)
      ordered_g  <- all_g[order(gene_ranks, na.last = TRUE)]
      utils::head(ordered_g, le_size)
    })
    names(le_list) <- common_sets
    result$le <- le_list
  }

  # If only one type requested: simplify but tag with attribute
  # so addgenesPA() knows which column name to use
  if (length(genes) == 1) {
    out <- result[[genes]]
    attr(out, "genes_type") <- genes
    return(out)
  }

  return(result)
}


######################
# Function addgenesPA #
######################

#' Add gene columns to pathway analysis results
#'
#' Appends `all_genes` and/or `le_genes` (GSEA output only) columns to a pathway
#' analysis results data frame based on the output of [getgenesPA()].
#' Gene symbols within each cell are comma-separated. Automatically detects
#' which column(s) to add based on the structure of the input.
#'
#' @param pa_data A data frame of pathway analysis results containing a `NAME`
#'   column. Typically the output of [merge_PA()].
#' @param gene_lists Output of [getgenesPA()]. Can be:
#'   * A list with `$all` and/or `$le` slots : when
#'     `getgenesPA(..., genes = c("all", "le"))`. Adds both `all_genes` and
#'     `le_genes` columns.
#'   * A flat named list with attribute `genes_type = "le"` : when
#'     `getgenesPA(..., genes = "le")`. Adds `le_genes` column only.
#'   * A flat named list with attribute `genes_type = "all"` : when
#'     `getgenesPA(..., genes = "all")`. Adds `all_genes` column only.
#'
#' @return The input `pa_data` data frame with one or two additional columns:
#'   * `all_genes`: comma-separated string of all gene set members.
#'   * `le_genes`: comma-separated string of leading edge genes ordered by
#'     rank (GSEA output only).
#'
#'   Gene sets not found in `gene_lists` receive `NA`.
#'
#' @examples
#' \dontrun{
#' gsl     <- list_gmts("path/to/gmt_folder/")
#' pa_data <- read.delim("gsea_results_merged.tsv")
#' ranked  <- deseq2_results$gene_id[order(deseq2_results$stat,
#'                                         decreasing = TRUE)]
#'
#' # Add both columns at once
#' gene_lists <- getgenesPA(pa_data, gsl, ranked, genes = c("all", "le"))
#' pa_annot   <- addgenesPA(pa_data, gene_lists)
#' head(pa_annot[, c("NAME", "all_genes", "le_genes")])
#'
#' # Add only leading edge genes
#' le_only  <- getgenesPA(pa_data, gsl, ranked, genes = "le")
#' pa_annot <- addgenesPA(pa_data, le_only)
#' head(pa_annot[, c("NAME", "le_genes")])
#' }
#'
#' @seealso [getgenesPA()] to generate `gene_lists`;
#'   [heatmap_PA()] for heatmap visualization;
#'   [save_results()] to export the annotated results.
#'
#' @export

addgenesPA <- function(pa_data, gene_lists) {

  if (!is.data.frame(pa_data)) {
    stop("`pa_data` must be a data frame.", call. = FALSE)
  }
  if (!"NAME" %in% colnames(pa_data)) {
    stop("`pa_data` must contain a column named 'NAME'.", call. = FALSE)
  }
  if (!is.list(gene_lists)) {
    stop("`gene_lists` must be a list : output of getgenesPA().", call. = FALSE)
  }

  # Helper: collapse gene vector to comma-separated string per gene set
  .collapse <- function(named_list, set_names) {
    vapply(set_names, function(gs) {
      if (!gs %in% names(named_list)) return(NA_character_)
      g <- named_list[[gs]]
      if (length(g) == 0) return(NA_character_)
      paste(g, collapse = ",")
    }, character(1))
  }

  # Detect input structure
  has_all_slot <- "all" %in% names(gene_lists) && is.list(gene_lists$all)
  has_le_slot  <- "le"  %in% names(gene_lists) && is.list(gene_lists$le)
  genes_type   <- attr(gene_lists, "genes_type")  # set by getgenesPA for single type

  if (has_all_slot || has_le_slot) {
    # Full output: list with $all and/or $le
    if (has_all_slot) pa_data$all_genes <- .collapse(gene_lists$all, pa_data$NAME)
    if (has_le_slot)  pa_data$le_genes  <- .collapse(gene_lists$le,  pa_data$NAME)

  } else if (!is.null(genes_type)) {
    # Single-type output tagged with attribute
    col_name <- if (genes_type == "le") "le_genes" else "all_genes"
    pa_data[[col_name]] <- .collapse(gene_lists, pa_data$NAME)

  } else {
    stop(
      "`gene_lists` structure not recognized. ",
      "Pass the direct output of getgenesPA().",
      call. = FALSE
    )
  }

  return(pa_data)
}
