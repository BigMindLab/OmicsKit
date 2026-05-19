##############################
# Function crossLayerCorr     #
##############################

#' Correlate paired sample profiles between two omics layers.
#'
#' Computes one correlation per shared sample between two normalized matrices.
#' Matrices are first matched by shared sample names and shared feature names,
#' making the function suitable for RNA-vs-protein comparisons after mapping
#' protein features to gene symbols.
#'
#' By default, the function uses the top variable shared features to reduce noise
#' and returns both a sorted correlation table and an optional bar plot.
#'
#' @param mat_x Numeric matrix for layer X with features as rows and samples as
#'   columns. Row names and column names are required.
#' @param mat_y Numeric matrix for layer Y with features as rows and samples as
#'   columns. Row names and column names are required.
#' @param method Correlation method. One of `"spearman"` or `"pearson"`.
#'   Default: `"spearman"`.
#' @param top_n Number of most variable shared features to use. If `NULL`, all
#'   shared features are used. Default: `500`.
#' @param plot Logical; if `TRUE`, returns a `ggplot2` bar plot in the `plot`
#'   element. Default: `TRUE`.
#'
#' @return A list with:
#'   \itemize{
#'   \item `correlations`: data frame with `sample` and `correlation` columns.
#'   \item `median_r`: median sample-level correlation.
#'   \item `n_shared_samples`: number of shared samples used.
#'   \item `n_shared_features`: number of shared features before top-variable
#'     filtering.
#'   \item `n_features_used`: number of features used for correlation.
#'   \item `plot`: `ggplot2` object or `NULL`.
#'   }
#'
#' @examples
#' set.seed(1)
#' mat_rna <- matrix(rnorm(60), nrow = 10)
#' mat_protein <- mat_rna + matrix(rnorm(60, sd = 0.4), nrow = 10)
#' rownames(mat_rna) <- rownames(mat_protein) <- paste0("GENE", seq_len(10))
#' colnames(mat_rna) <- colnames(mat_protein) <- paste0("S", seq_len(6))
#'
#' res <- crossLayerCorr(mat_rna, mat_protein, top_n = 8, plot = FALSE)
#' head(res$correlations)
#' res$median_r
#'
#' @import ggplot2
#' @importFrom rlang .data
#' 
#' @seealso
#' [concordanceDE()] to classify genes by cross-layer differential concordance;
#' [nice_ConcordanceScatter()] to visualize RNA-protein differential
#' concordance.
#' 
#' @export
crossLayerCorr <- function(
  mat_x,
  mat_y,
  method = c("spearman", "pearson"),
  top_n = 500,
  plot = TRUE
) {
  method <- match.arg(method)

  mat_x <- as.matrix(mat_x)
  mat_y <- as.matrix(mat_y)

  if (is.null(colnames(mat_x)) || is.null(colnames(mat_y))) {
    stop("Both matrices must have sample names as column names.", call. = FALSE)
  }

  if (is.null(rownames(mat_x)) || is.null(rownames(mat_y))) {
    stop("Both matrices must have feature names as row names.", call. = FALSE)
  }

  storage.mode(mat_x) <- "numeric"
  storage.mode(mat_y) <- "numeric"

  shared_samples <- intersect(colnames(mat_x), colnames(mat_y))

  if (length(shared_samples) == 0) {
    stop("No shared sample names between `mat_x` and `mat_y`.", call. = FALSE)
  }

  shared_features <- intersect(rownames(mat_x), rownames(mat_y))

  if (length(shared_features) < 2) {
    stop(
      "At least two shared feature names are required between `mat_x` and `mat_y`.",
      call. = FALSE
    )
  }

  mat_x <- mat_x[shared_features, shared_samples, drop = FALSE]
  mat_y <- mat_y[shared_features, shared_samples, drop = FALSE]

  var_x <- apply(mat_x, 1, stats::var, na.rm = TRUE)
  var_y <- apply(mat_y, 1, stats::var, na.rm = TRUE)
  combined_var <- rowMeans(cbind(var_x, var_y), na.rm = TRUE)

  keep_features <- is.finite(combined_var) & combined_var > 0

  if (sum(keep_features) < 2) {
    stop("Too few variable shared features for correlation.", call. = FALSE)
  }

  mat_x <- mat_x[keep_features, , drop = FALSE]
  mat_y <- mat_y[keep_features, , drop = FALSE]
  combined_var <- combined_var[keep_features]

  n_shared_features <- length(shared_features)

  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 2) {
      stop("`top_n` must be `NULL` or a numeric value greater than or equal to 2.", call. = FALSE)
    }

    top_idx <- order(combined_var, decreasing = TRUE)[seq_len(min(top_n, length(combined_var)))]
    mat_x <- mat_x[top_idx, , drop = FALSE]
    mat_y <- mat_y[top_idx, , drop = FALSE]
  }

  cors <- vapply(shared_samples, function(sample_id) {
    x <- mat_x[, sample_id]
    y <- mat_y[, sample_id]
    keep <- is.finite(x) & is.finite(y)

    if (sum(keep) < 2) {
      return(NA_real_)
    }

    stats::cor(x[keep], y[keep], method = method, use = "complete.obs")
  }, numeric(1))

  result <- data.frame(
    sample = shared_samples,
    correlation = cors,
    stringsAsFactors = FALSE
  )

  result <- result[order(result$correlation, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  result$sample_ordered <- factor(result$sample, levels = result$sample)

  median_r <- stats::median(cors, na.rm = TRUE)

  p <- NULL

  if (isTRUE(plot)) {
    p <- ggplot2::ggplot(
      result,
      ggplot2::aes(
        x = .data[["sample_ordered"]],
        y = .data[["correlation"]]
      )
    ) +
      ggplot2::geom_col(fill = "#7F77DD", alpha = 0.85, width = 0.75) +
      ggplot2::geom_hline(
        yintercept = median_r,
        linetype = "dashed",
        color = "#EF9F27"
      ) +
      ggplot2::labs(
        x = "Sample",
        y = paste0(tools::toTitleCase(method), " r"),
        caption = paste0(
          "Median r = ", round(median_r, 3),
          " | n = ", length(shared_samples), " samples",
          " | features = ", nrow(mat_x)
        )
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7)
      )
  }

  list(
    correlations = result,
    median_r = median_r,
    n_shared_samples = length(shared_samples),
    n_shared_features = n_shared_features,
    n_features_used = nrow(mat_x),
    plot = p
  )
}
