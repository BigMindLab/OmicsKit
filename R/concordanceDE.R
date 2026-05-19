##############################
# Function concordanceDE      #
##############################

#' Classify differential results by cross-layer concordance.
#'
#' Compares two differential-analysis tables from paired omics layers, such as
#' RNA-seq and total-protein RPPA, by merging shared genes and classifying each
#' gene according to statistical significance and log-fold-change direction.
#'
#' Genes are assigned to one of five categories: `concordant`, `discordant`,
#' `only_x`, `only_y`, or `not_significant`. A concordance score is also
#' computed as the fraction of genes significant in both layers that have the
#' same effect direction.
#'
#' @param de_x Data frame for layer X. Must contain columns defined by
#'   `gene_col`, `logfc_col`, and `padj_col`.
#' @param de_y Data frame for layer Y. Must contain columns defined by
#'   `gene_col`, `logfc_col`, and `padj_col`.
#' @param gene_col Column name containing gene identifiers. Default: `"gene"`.
#' @param logfc_col Column name containing log-fold changes. Default: `"logFC"`.
#' @param padj_col Column name containing adjusted p-values/FDR values.
#'   Default: `"padj"`.
#' @param padj_threshold Numeric threshold for adjusted p-values. Use a single
#'   value for both layers or a length-two vector for layer X and layer Y,
#'   respectively. Default: `0.05`.
#' @param logfc_threshold Numeric absolute log-fold-change threshold. Use a
#'   single value for both layers or a length-two vector for layer X and layer Y,
#'   respectively. Default: `1`.
#'
#' @return A list with:
#'   \itemize{
#'   \item `table`: merged data frame with one `category` column.
#'   \item `concordance_score`: fraction of genes significant in both layers
#'     with matching log-fold-change direction.
#'   \item `summary`: table with counts per concordance category.
#'   \item `n_shared_genes`: number of genes shared by both layers.
#'   \item `n_both_significant`: number of genes significant in both layers.
#'   }
#'
#' @examples
#' de_rna <- data.frame(
#'   gene = c("ESR1", "PGR", "EGFR", "MKI67", "GATA3"),
#'   logFC = c(3.2, 2.1, -1.6, -1.4, 1.8),
#'   padj = c(1e-6, 0.002, 0.01, 0.03, 0.04)
#' )
#'
#' de_protein <- data.frame(
#'   gene = c("ESR1", "PGR", "EGFR", "MKI67", "GATA3"),
#'   logFC = c(0.7, 0.4, -0.5, 0.3, 0.05),
#'   padj = c(1e-4, 0.01, 0.02, 0.04, 0.40)
#' )
#'
#' res <- concordanceDE(
#'   de_x = de_rna,
#'   de_y = de_protein,
#'   logfc_threshold = c(1, 0.2)
#' )
#' res$summary
#' res$concordance_score
#'
#' @importFrom dplyr case_when
#' 
#' @seealso
#' [nice_ConcordanceScatter()] to visualize the output of `concordanceDE()`;
#' [crossLayerCorr()] to estimate sample-level concordance between two omics
#' layers.
#'
#' @export
concordanceDE <- function(
  de_x,
  de_y,
  gene_col = "gene",
  logfc_col = "logFC",
  padj_col = "padj",
  padj_threshold = 0.05,
  logfc_threshold = 1
) {
  de_x <- as.data.frame(de_x)
  de_y <- as.data.frame(de_y)

  required_cols <- c(gene_col, logfc_col, padj_col)
  missing_x <- setdiff(required_cols, names(de_x))
  missing_y <- setdiff(required_cols, names(de_y))

  if (length(missing_x) > 0) {
    stop(
      "`de_x` is missing required columns: ",
      paste(missing_x, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_y) > 0) {
    stop(
      "`de_y` is missing required columns: ",
      paste(missing_y, collapse = ", "),
      call. = FALSE
    )
  }

  de_x <- de_x[!is.na(de_x[[gene_col]]) & nzchar(trimws(de_x[[gene_col]])), , drop = FALSE]
  de_y <- de_y[!is.na(de_y[[gene_col]]) & nzchar(trimws(de_y[[gene_col]])), , drop = FALSE]

  if (anyDuplicated(de_x[[gene_col]]) > 0) {
    stop("`de_x` must contain one row per gene in `gene_col`.", call. = FALSE)
  }

  if (anyDuplicated(de_y[[gene_col]]) > 0) {
    stop("`de_y` must contain one row per gene in `gene_col`.", call. = FALSE)
  }

  if (!is.numeric(padj_threshold) || length(padj_threshold) > 2) {
    stop("`padj_threshold` must be a numeric vector of length 1 or 2.", call. = FALSE)
  }

  if (!is.numeric(logfc_threshold) || length(logfc_threshold) > 2) {
    stop("`logfc_threshold` must be a numeric vector of length 1 or 2.", call. = FALSE)
  }

  padj_threshold <- rep(padj_threshold, length.out = 2)
  logfc_threshold <- rep(logfc_threshold, length.out = 2)

  shared <- merge(de_x, de_y, by = gene_col, suffixes = c("_x", "_y"))

  if (nrow(shared) == 0) {
    stop("No shared genes were found between `de_x` and `de_y`.", call. = FALSE)
  }

  lfc_x_col <- paste0(logfc_col, "_x")
  lfc_y_col <- paste0(logfc_col, "_y")
  padj_x_col <- paste0(padj_col, "_x")
  padj_y_col <- paste0(padj_col, "_y")

  lfc_x <- suppressWarnings(as.numeric(shared[[lfc_x_col]]))
  lfc_y <- suppressWarnings(as.numeric(shared[[lfc_y_col]]))
  padj_x <- suppressWarnings(as.numeric(shared[[padj_x_col]]))
  padj_y <- suppressWarnings(as.numeric(shared[[padj_y_col]]))

  sig_x <- !is.na(lfc_x) & !is.na(padj_x) &
    abs(lfc_x) >= logfc_threshold[1] & padj_x < padj_threshold[1]
  sig_y <- !is.na(lfc_y) & !is.na(padj_y) &
    abs(lfc_y) >= logfc_threshold[2] & padj_y < padj_threshold[2]

  dir_x <- sign(lfc_x)
  dir_y <- sign(lfc_y)

  shared$category <- dplyr::case_when(
    sig_x & sig_y & dir_x == dir_y ~ "concordant",
    sig_x & sig_y & dir_x != dir_y ~ "discordant",
    sig_x & !sig_y ~ "only_x",
    !sig_x & sig_y ~ "only_y",
    TRUE ~ "not_significant"
  )

  category_levels <- c(
    "concordant",
    "discordant",
    "only_x",
    "only_y",
    "not_significant"
  )

  shared$category <- factor(shared$category, levels = category_levels)

  both_sig <- sig_x & sig_y
  concordance_score <- if (any(both_sig)) {
    mean(dir_x[both_sig] == dir_y[both_sig], na.rm = TRUE)
  } else {
    NA_real_
  }

  out <- list(
    table = shared,
    concordance_score = round(concordance_score, 4),
    summary = table(shared$category),
    n_shared_genes = nrow(shared),
    n_both_significant = sum(both_sig, na.rm = TRUE)
  )

  class(out) <- c("omicskit_concordanceDE", "list")
  out
}
