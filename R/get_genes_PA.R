######################
# Function getgenesPA #
######################

#' Extract gene members from pathway analysis results
#'
#' For each gene set in a pathway analysis results table, retrieves leading
#' edge genes, a user-defined top fraction of genes, all genes in the gene
#' set, or any combination. All gene lists are ordered by their rank in the
#' provided ranked gene list.
#'
#' **Three extraction modes:**
#'
#' * `"le"`: **GSEA output only.** Leading edge genes: the subset of genes
#'   that drives the enrichment signal. Size is computed as
#'   `round(SIZE * tags)`, where `tags` is the fraction of gene hits before
#'   (positive ES) or after (negative ES) the peak in the running enrichment
#'   score. Definition from the GSEA User Guide: *"The percentage of gene hits
#'   before (for positive ES) or after (for negative ES) the peak in the
#'   running enrichment score. This gives an indication of the percentage of
#'   genes contributing to the enrichment score."*
#'   (\url{https://docs.gsea-msigdb.org/#GSEA/GSEA_User_Guide/}).
#'   Requires columns `SIZE` and `tags` in `pa_data`, produced automatically
#'   by [merge_PA()].
#'
#' * `"top"`: **Any enrichment result (GSEA, CAMERA, PADOG, etc.).**
#'   A user-defined top fraction of genes ordered by rank. Size is computed as
#'   `round(SIZE * top)`, where `top` is a numeric value between 0 and 1
#'   provided in a `top` column of `pa_data`. This does **not** represent true
#'   leading edge genes: it is a flexible, rank-based selection suitable for
#'   exploratory visualization with any pathway analysis method.
#'   Requires columns `SIZE` and `top` in `pa_data`.
#'
#' * `"all"`: All genes in the gene set, ordered by rank.
#'
#' @param pa_data A data frame of pathway analysis results. Must always
#'   contain:
#'   * `NAME`: gene set name.
#'
#'   Additionally required depending on `genes`:
#'   * `SIZE`: number of genes in the gene set. Required for `"le"` and
#'     `"top"`.
#'   * `tags`: numeric fraction (0-1) of genes contributing to the
#'     enrichment score (GSEA leading edge). Produced automatically by
#'     [merge_PA()]. Required for `genes = "le"`.
#'   * `top`: numeric fraction (0-1) defining the proportion of top-ranked
#'     genes to extract. Set manually by the user (e.g.,
#'     `pa_data$top <- 0.25` for the top 25%). Required for `genes = "top"`.
#'
#'   Typically the output of [merge_PA()].
#' @param geneset_list A named list of gene sets, where each element is a
#'   character vector of gene symbols. Typically the output of [list_gmts()],
#'   or use the built-in [geneset_list] for quick testing.
#' @param ranked_genes A character vector of gene symbols ordered by their
#'   ranking metric (e.g., DESeq2 `stat`, log2FC, or signal-to-noise ratio),
#'   from most positive to most negative. Non-significant genes fall in the
#'   middle of the list. Used to order genes within each extracted set.
#' @param genes Character vector specifying which extraction mode(s) to use.
#'   Any combination of `"all"`, `"le"`, and `"top"`. Default:
#'   `c("all", "le")`.
#'
#' @return Depends on `genes`:
#'   * Single mode (e.g., `genes = "le"`): a named list where each element
#'     is a character vector of gene symbols. The list has an attribute
#'     `genes_type` used by [addgenesPA()] to name the output column.
#'   * Multiple modes (e.g., `genes = c("all", "le", "top")`): a named list
#'     with one element per requested mode:
#'     * `$all`: named list of all gene set members.
#'     * `$le`: named list of leading edge genes (GSEA only).
#'     * `$top`: named list of top-ranked genes.
#'
#' @examples
#' \dontrun{
#' data(gsea_results)
#' data(geneset_list)
#' data(deseq2_results)
#'
#' #or
#' gsl <- list_gmts("path/to/gmt_folder/")
#'
#' ranked    <- deseq2_results$gene_id[order(deseq2_results$stat,
#'                                           decreasing = TRUE)]
#' pa_single <- gsea_results[gsea_results$COMPARISON == "TumorVsNormal", ]
#'
#' # ── GSEA results: all three modes available
#' gene_lists <- getgenesPA(pa_single, geneset_list, ranked,
#'                          genes = c("all", "le", "top"))
#'
#' # But first add the top column (e.g. top 30% of genes by rank)
#' pa_single$top <- 0.30
#' gene_lists <- getgenesPA(pa_single, geneset_list, ranked,
#'                          genes = c("all", "le", "top"))
#'
#' gene_lists$le[["KEGG_APOPTOSIS"]]    # leading edge genes
#' gene_lists$top[["KEGG_APOPTOSIS"]]   # top 30% by rank
#' gene_lists$all[["KEGG_APOPTOSIS"]]   # all genes
#'
#' pa_annot <- addgenesPA(pa_single, gene_lists)
#' head(pa_annot[, c("NAME", "all_genes", "le_genes", "top_genes")])
#'
#' # ── CAMERA results: use "top" (no leading edge available) ───
#' data(camera_results)
#' camera_pa      <- camera_results
#' colnames(camera_pa)[colnames(camera_pa) == "GeneSet"] <- "NAME"
#' camera_pa$SIZE <- sapply(camera_pa$NAME,
#'                          function(x) length(geneset_list[[x]]))
#' camera_pa$top  <- 0.25   # top 25% by rank
#'
#' gene_lists_cam <- getgenesPA(camera_pa, geneset_list, ranked,
#'                              genes = c("all", "top"))
#' pa_annot_cam   <- addgenesPA(camera_pa, gene_lists_cam)
#' head(pa_annot_cam[, c("NAME", "all_genes", "top_genes")])
#' }
#'
#' @seealso [addgenesPA()] to append gene columns to pa_data;
#'   [heatmap_PA()] for heatmap visualization;
#'   [list_gmts()] to generate `geneset_list`;
#'   [merge_PA()] to generate `pa_data` with the required `tags` column.
#'
#' @export

getgenesPA <- function(pa_data, geneset_list, ranked_genes,
                       genes = c("all", "le")) {

  genes <- match.arg(genes, choices = c("all", "le", "top"), several.ok = TRUE)

  # --- Input validation ----
  if (!is.data.frame(pa_data)) {
    stop("`pa_data` must be a data frame.", call. = FALSE)
  }
  if (!"NAME" %in% colnames(pa_data)) {
    stop("`pa_data` must contain a column named 'NAME'.", call. = FALSE)
  }
  if ("le" %in% genes && !all(c("SIZE", "tags") %in% colnames(pa_data))) {
    stop(
      "`pa_data` must contain columns 'SIZE' and 'tags' when genes = 'le'. ",
      "These are produced by merge_PA() from GSEA output. ",
      "For non-GSEA results, use genes = 'top' and add: ",
      "pa_data$SIZE <- ...; pa_data$top <- 0.25",
      call. = FALSE
    )
  }
  if ("top" %in% genes && !all(c("SIZE", "top") %in% colnames(pa_data))) {
    stop(
      "`pa_data` must contain columns 'SIZE' and 'top' when genes = 'top'. ",
      "Add them manually: pa_data$SIZE <- ...; pa_data$top <- 0.25",
      call. = FALSE
    )
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

  # --- Internal helper: order genes by rank and take top n -
  .extract_top_n <- function(all_g, n, ranked_genes) {
    if (is.na(n) || n <= 0) return(character(0))
    gene_ranks <- match(all_g, ranked_genes)
    ordered_g  <- all_g[order(gene_ranks, na.last = TRUE)]
    utils::head(ordered_g, n)
  }

  result <- list()

  # --- all: every gene in the gene set ordered by rank
  if ("all" %in% genes) {
    result$all <- stats::setNames(
      lapply(common_sets, function(gs) {
        all_g      <- geneset_list[[gs]]
        gene_ranks <- match(all_g, ranked_genes)
        all_g[order(gene_ranks, na.last = TRUE)]
      }),
      common_sets
    )
  }

  # --- le: leading edge genes (GSEA only, uses tags column)
  if ("le" %in% genes) {
    result$le <- stats::setNames(
      lapply(common_sets, function(gs) {
        row    <- pa_data[pa_data$NAME == gs, ][1, ]
        all_g  <- geneset_list[[gs]]
        # Leading edge size from tags (fraction contributing to ES peak)
        le_size <- round(as.numeric(row$SIZE) * as.numeric(row$tags))
        .extract_top_n(all_g, le_size, ranked_genes)
      }),
      common_sets
    )
  }

  # --- top: user-defined top fraction by rank (any enrichment method) -
  if ("top" %in% genes) {
    result$top <- stats::setNames(
      lapply(common_sets, function(gs) {
        row      <- pa_data[pa_data$NAME == gs, ][1, ]
        all_g    <- geneset_list[[gs]]
        top_size <- round(as.numeric(row$SIZE) * as.numeric(row$top))
        .extract_top_n(all_g, top_size, ranked_genes)
      }),
      common_sets
    )
  }

  # --- Simplify output if only one mode requested -----
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
#' Appends `all_genes`, `le_genes`, and/or `top_genes` columns to a pathway
#' analysis results data frame based on the output of [getgenesPA()].
#' Gene symbols within each cell are comma-separated. Automatically detects
#' which column(s) to add based on the structure of the input.
#'
#' @param pa_data A data frame of pathway analysis results containing a `NAME`
#'   column. Typically the output of [merge_PA()].
#' @param gene_lists Output of [getgenesPA()]. Can be:
#'   * A list with `$all`, `$le`, and/or `$top` slots: when multiple modes
#'     are requested (e.g., `getgenesPA(..., genes = c("all", "le", "top"))`).
#'     Adds the corresponding columns.
#'   * A flat named list with attribute `genes_type`: when a single mode is
#'     requested. Adds the corresponding column (`all_genes`, `le_genes`, or
#'     `top_genes`).
#'
#' @return The input `pa_data` data frame with one or more additional columns:
#'   * `all_genes`: comma-separated string of all gene set members ordered by
#'     rank.
#'   * `le_genes`: comma-separated string of leading edge genes (GSEA only),
#'     ordered by rank.
#'   * `top_genes`: comma-separated string of top-ranked genes based on the
#'     user-defined `top` fraction.
#'
#'   Gene sets not found in `gene_lists` receive `NA`.
#'
#' @examples
#' \dontrun{
#' data(gsea_results)
#' data(geneset_list)
#' data(deseq2_results)
#'
#' ranked    <- deseq2_results$gene_id[order(deseq2_results$stat,
#'                                           decreasing = TRUE)]
#' pa_single <- gsea_results[gsea_results$COMPARISON == "TumorVsNormal", ]
#' pa_single$top <- 0.30
#'
#' # Add all three columns
#' gene_lists <- getgenesPA(pa_single, geneset_list, ranked,
#'                          genes = c("all", "le", "top"))
#' pa_annot   <- addgenesPA(pa_single, gene_lists)
#' head(pa_annot[, c("NAME", "all_genes", "le_genes", "top_genes")])
#'
#' # Add only leading edge genes
#' le_only  <- getgenesPA(pa_single, geneset_list, ranked, genes = "le")
#' pa_annot <- addgenesPA(pa_single, le_only)
#' head(pa_annot[, c("NAME", "le_genes")])
#'
#' # CAMERA: add only top and all (no leading edge)
#' data(camera_results)
#' camera_pa      <- camera_results
#' colnames(camera_pa)[colnames(camera_pa) == "GeneSet"] <- "NAME"
#' camera_pa$SIZE <- sapply(camera_pa$NAME,
#'                          function(x) length(geneset_list[[x]]))
#' camera_pa$top  <- 0.25
#' gene_lists_cam <- getgenesPA(camera_pa, geneset_list, ranked,
#'                              genes = c("all", "top"))
#' pa_annot_cam   <- addgenesPA(camera_pa, gene_lists_cam)
#' head(pa_annot_cam[, c("NAME", "all_genes", "top_genes")])
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
    stop("`gene_lists` must be a list, output of getgenesPA().", call. = FALSE)
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

  # Column name mapping per genes_type
  col_map <- c(all = "all_genes", le = "le_genes", top = "top_genes")

  has_all_slot <- "all" %in% names(gene_lists) && is.list(gene_lists$all)
  has_le_slot  <- "le"  %in% names(gene_lists) && is.list(gene_lists$le)
  has_top_slot <- "top" %in% names(gene_lists) && is.list(gene_lists$top)
  genes_type   <- attr(gene_lists, "genes_type")

  if (has_all_slot || has_le_slot || has_top_slot) {
    if (has_all_slot) pa_data$all_genes <- .collapse(gene_lists$all, pa_data$NAME)
    if (has_le_slot)  pa_data$le_genes  <- .collapse(gene_lists$le,  pa_data$NAME)
    if (has_top_slot) pa_data$top_genes <- .collapse(gene_lists$top, pa_data$NAME)
  } else if (!is.null(genes_type) && genes_type %in% names(col_map)) {
    pa_data[[col_map[[genes_type]]]] <- .collapse(gene_lists, pa_data$NAME)
  } else {
    stop(
      "`gene_lists` structure not recognized. ",
      "Pass the direct output of getgenesPA().",
      call. = FALSE
    )
  }

  return(pa_data)
}
