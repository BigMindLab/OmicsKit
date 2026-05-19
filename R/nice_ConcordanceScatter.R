####################################
# Function nice_ConcordanceScatter  #
####################################

#' Plot cross-layer log-fold-change concordance.
#'
#' Builds a scatter plot comparing log-fold changes from two omics layers after
#' classification with [concordanceDE()]. Points are colored by concordance
#' category, reference lines divide the plot into directional quadrants, and an
#' optional regression line is drawn for concordant genes.
#'
#' @param concordance_result Output from [concordanceDE()]. A data frame with a
#'   `category` column can also be supplied through `concordance_result$table`.
#' @param x_label X-axis label. Default: `"log2FC (RNA)"`.
#' @param y_label Y-axis label. Default: `"log2FC (Protein)"`.
#' @param genes_label Optional character vector of genes to label with
#'   `ggrepel`.
#' @param method Correlation method used by `ggpubr::stat_cor()`. One of
#'   `"pearson"` or `"spearman"`. Default: `"pearson"`.
#' @param ci_level Confidence level for the linear-model confidence interval.
#'   Default: `0.95`.
#' @param point_size Point size. Default: `1.5`.
#' @param alpha Point transparency. Default: `0.65`.
#' @param gene_col Optional gene column name. If `NULL`, the function uses the
#'   first column containing `"gene"`, ignoring case.
#' @param logfc_col Base log-fold-change column name used in [concordanceDE()].
#'   Default: `"logFC"`, which expects columns `logFC_x` and `logFC_y`.
#' @param category_colors Optional named character vector with colors for
#'   `concordant`, `discordant`, `only_x`, `only_y`, and `not_significant`.
#'
#' @return A `ggplot2` object.
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
#'
#' if (
#'   requireNamespace("ggpubr", quietly = TRUE) &&
#'     requireNamespace("ggrepel", quietly = TRUE)
#' ) {
#'   nice_ConcordanceScatter(
#'     res,
#'     genes_label = c("ESR1", "EGFR"),
#'     method = "spearman"
#'   )
#' }
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @seealso
#' [concordanceDE()] to classify genes by cross-layer differential concordance;
#' [crossLayerCorr()] to estimate sample-level concordance between two omics
#' layers.
#' 
#' @export
nice_ConcordanceScatter <- function(
  concordance_result,
  x_label = "log2FC (RNA)",
  y_label = "log2FC (Protein)",
  genes_label = NULL,
  method = c("pearson", "spearman"),
  ci_level = 0.95,
  point_size = 1.5,
  alpha = 0.65,
  gene_col = NULL,
  logfc_col = "logFC",
  category_colors = NULL
) {
  method <- match.arg(method)

  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop(
      "Package \"ggpubr\" must be installed to use `nice_ConcordanceScatter()`.",
      call. = FALSE
    )
  }

  if (!is.null(genes_label) && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "Package \"ggrepel\" must be installed when `genes_label` is used.",
      call. = FALSE
    )
  }

  if (!is.list(concordance_result) || is.null(concordance_result$table)) {
    stop("`concordance_result` must be the output of `concordanceDE()`.", call. = FALSE)
  }

  df <- as.data.frame(concordance_result$table)

  if (!"category" %in% names(df)) {
    stop("`concordance_result$table` must contain a `category` column.", call. = FALSE)
  }

  lfc_x <- paste0(logfc_col, "_x")
  lfc_y <- paste0(logfc_col, "_y")

  if (!all(c(lfc_x, lfc_y) %in% names(df))) {
    lfc_candidates_x <- grep("_x$", names(df), value = TRUE)
    lfc_candidates_y <- grep("_y$", names(df), value = TRUE)

    lfc_x <- lfc_candidates_x[grepl("log", lfc_candidates_x, ignore.case = TRUE)][1]
    lfc_y <- lfc_candidates_y[grepl("log", lfc_candidates_y, ignore.case = TRUE)][1]
  }

  if (is.na(lfc_x) || is.na(lfc_y) || !all(c(lfc_x, lfc_y) %in% names(df))) {
    stop(
      "Could not identify log-fold-change columns. Expected `",
      paste0(logfc_col, "_x"), "` and `", paste0(logfc_col, "_y"), "`.",
      call. = FALSE
    )
  }

  if (is.null(gene_col)) {
    gene_col <- grep("gene", names(df), ignore.case = TRUE, value = TRUE)[1]
  }

  default_colors <- c(
    concordant = "#1D9E75",
    discordant = "#D85A30",
    only_x = "#7F77DD",
    only_y = "#EF9F27",
    not_significant = "#B4B2A9"
  )

  if (is.null(category_colors)) {
    category_colors <- default_colors
  } else {
    category_colors <- utils::modifyList(default_colors, as.list(category_colors))
    category_colors <- unlist(category_colors, use.names = TRUE)
  }

  df$category <- factor(
    as.character(df$category),
    levels = names(default_colors)
  )

  category_labels <- c(
    concordant = paste0("Concordant (score=", concordance_result$concordance_score, ")"),
    discordant = "Discordant",
    only_x = "Only layer X",
    only_y = "Only layer Y",
    not_significant = "NS"
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[lfc_x]],
      y = .data[[lfc_y]],
      color = .data[["category"]]
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.4
    ) +
    ggplot2::geom_point(alpha = alpha, size = point_size) +
    ggpubr::stat_cor(
      ggplot2::aes(x = .data[[lfc_x]], y = .data[[lfc_y]]),
      method = method,
      label.x.npc = 0.05,
      inherit.aes = FALSE,
      size = 3
    ) +
    ggplot2::scale_color_manual(
      values = category_colors,
      labels = category_labels,
      drop = FALSE
    ) +
    ggplot2::labs(x = x_label, y = y_label, color = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "right")

  df_concordant <- df[df$category == "concordant", , drop = FALSE]

  var_concordant_x <- stats::var(df_concordant[[lfc_x]], na.rm = TRUE)
  var_concordant_y <- stats::var(df_concordant[[lfc_y]], na.rm = TRUE)

  if (
    nrow(df_concordant) >= 2 &&
      is.finite(var_concordant_x) && var_concordant_x > 0 &&
      is.finite(var_concordant_y) && var_concordant_y > 0
  ) {
    p <- p +
      ggplot2::geom_smooth(
        data = df_concordant,
        ggplot2::aes(x = .data[[lfc_x]], y = .data[[lfc_y]]),
        method = "lm",
        se = TRUE,
        level = ci_level,
        color = "#085041",
        fill = "#9FE1CB",
        linewidth = 0.8,
        inherit.aes = FALSE
      )
  }

  if (!is.null(genes_label)) {
    if (is.na(gene_col) || !gene_col %in% names(df)) {
      stop("Could not identify a gene column for `genes_label`.", call. = FALSE)
    }

    df_lbl <- df[df[[gene_col]] %in% genes_label, , drop = FALSE]

    if (nrow(df_lbl) > 0) {
      p <- p +
        ggrepel::geom_label_repel(
          data = df_lbl,
          ggplot2::aes(
            x = .data[[lfc_x]],
            y = .data[[lfc_y]],
            label = .data[[gene_col]]
          ),
          inherit.aes = FALSE,
          size = 3,
          max.overlaps = 20,
          fill = scales::alpha("white", 0.8)
        )
    }
  }

  p
}
