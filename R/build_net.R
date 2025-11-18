###############################
# Build enrichment clustering #
###############################

#' Calculate gene set clusters from GMT files and enrichment results
#'
#' This function parses one or more GMT files, filters an enrichment results
#' data frame by an FDR column, computes a Jaccard similarity matrix between
#' selected gene sets, performs hierarchical clustering, and determines the
#' optimal number of clusters using the silhouette index.
#'
#' @param results_df A \code{data.frame} or \code{tibble} containing enrichment
#'   results. Must contain at least the columns specified by \code{fdr_col} and
#'   \code{pathway_col}.
#' @param gmt_path A character scalar. Path to a directory containing one or
#'   more \code{.gmt} files.
#' @param fdr_threshold A numeric scalar specifying the FDR cutoff used to
#'   filter \code{results_df}. Only pathways with FDR values less than this
#'   threshold will be considered. Default is \code{0.25}.
#' @param fdr_col A character scalar specifying the column name in
#'   \code{results_df} that contains FDR values.
#' @param pathway_col A character scalar specifying the column name in
#'   \code{results_df} that contains pathway identifiers (names matching the
#'   GMT gene set names).
#' @param max_k_ratio A numeric scalar used to determine the maximum number of
#'   clusters to evaluate with the silhouette index. The maximum k is computed
#'   as \code{max(2, floor(n / max_k_ratio))}, where \code{n} is the number of
#'   selected gene sets. Default is \code{2.5}.
#'
#' @details
#' GMT files are expected to be tab-delimited, with the first column being the
#' gene set name, the second column a description (which is ignored), and the
#' remaining columns listing gene identifiers.
#'
#' Jaccard similarity between gene sets \eqn{A} and \eqn{B} is defined as:
#' \deqn{J(A,B) = |A \cap B| / |A \cup B|}.
#'
#' Hierarchical clustering is performed using Euclidean distances derived from
#' \code{1 - Jaccard_similarity} and the \code{"ward.D2"} linkage method.
#' The optimal number of clusters is chosen as the \code{k} that maximizes the
#' mean silhouette width for \code{k} in \code{2:max_k}.
#'
#' This function does not perform any file I/O beyond reading GMT files; it is
#' intended to be used inside packages and downstream plotting functions.
#'
#' @return A named list with the following elements:
#' \itemize{
#'   \item \code{clusters}: A \code{tibble} with columns \code{pathway} and
#'     \code{cluster} (integer cluster memberships).
#'   \item \code{jaccard_matrix}: A symmetric numeric matrix of Jaccard
#'     similarities between gene sets (rows and columns named by pathways).
#'   \item \code{hclust_obj}: An \code{\link[stats]{hclust}} object containing
#'     the hierarchical clustering result.
#'   \item \code{optimal_k}: An integer giving the selected number of clusters.
#'   \item \code{silhouette_data}: A \code{tibble} with columns \code{k} and
#'     \code{silhouette}, suitable for plotting the silhouette index as a
#'     function of \code{k}.
#' }
#'
#' @examples
#' \dontrun{
#' res <- read.csv("pathway_results.csv")
#' clustering <- build_net(
#'   results_df   = res,
#'   gmt_path     = "gmt_dir/",
#'   fdr_threshold = 0.25,
#'   fdr_col       = "FDR_gsea",
#'   pathway_col   = "Pathway",
#'   max_k_ratio   = 2.5
#' )
#' }
#'
#' @importFrom stats hclust dist cutree
#' @importFrom utils read.delim
#' @importFrom cluster silhouette
#' @importFrom dplyr filter pull distinct
#' @importFrom dplyr .data
#' @importFrom rlang .data
#' @importFrom tibble tibble
#' @export
build_net <- function(results_df,
                                       gmt_path,
                                       fdr_threshold = 0.25,
                                       fdr_col,
                                       pathway_col,
                                       max_k_ratio = 2.5) {
  # Basic input checks
  if (!dir.exists(gmt_path)) {
    stop("The provided 'gmt_path' does not exist: ", gmt_path)
  }
  if (!is.data.frame(results_df)) {
    stop("'results_df' must be a data.frame or tibble.")
  }
  if (!is.character(fdr_col) || length(fdr_col) != 1L) {
    stop("'fdr_col' must be a single character string.")
  }
  if (!is.character(pathway_col) || length(pathway_col) != 1L) {
    stop("'pathway_col' must be a single character string.")
  }
  if (!fdr_col %in% colnames(results_df)) {
    stop("Column '", fdr_col, "' not found in 'results_df'.")
  }
  if (!pathway_col %in% colnames(results_df)) {
    stop("Column '", pathway_col, "' not found in 'results_df'.")
  }
  if (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1L || is.na(fdr_threshold)) {
    stop("'fdr_threshold' must be a single numeric value.")
  }
  if (!is.numeric(max_k_ratio) || length(max_k_ratio) != 1L || is.na(max_k_ratio) || max_k_ratio <= 0) {
    stop("'max_k_ratio' must be a single positive numeric value.")
  }
  
  # Find GMT files
  gmt_files <- list.files(gmt_path, pattern = "\\.gmt$", full.names = TRUE)
  if (length(gmt_files) == 0L) {
    stop("No .gmt files found in directory: ", gmt_path)
  }
  
  # Read GMT files into a named list of gene vectors
  geneset_list <- list()
  for (f in gmt_files) {
    gmt <- utils::read.delim(f,
                             header = FALSE,
                             stringsAsFactors = FALSE,
                             sep = "\t",
                             quote = "",
                             comment.char = "")
    if (ncol(gmt) < 3L) {
      next
    }
    for (i in seq_len(nrow(gmt))) {
      name <- as.character(gmt[i, 1])
      genes <- as.character(gmt[i, 3:ncol(gmt)])
      genes <- genes[genes != "" & !is.na(genes)]
      geneset_list[[name]] <- unique(genes)
    }
  }
  
  if (length(geneset_list) == 0L) {
    stop("No gene sets could be parsed from the GMT files in '", gmt_path, "'.")
  }
  
  # Filter enrichment results by FDR and select pathways
  results_tbl <- if (inherits(results_df, "tbl_df")) results_df else tibble::as_tibble(results_df)
  
  results_filtered <- dplyr::filter(
    results_tbl,
    .data[[fdr_col]] < fdr_threshold
  )
  
  if (nrow(results_filtered) == 0L) {
    stop("No pathways passed the FDR threshold (", fdr_threshold, ").")
  }
  
  selected_sets <- dplyr::pull(results_filtered, .data[[pathway_col]])
  selected_sets <- unique(as.character(selected_sets))
  
  # Intersect with gene sets available in GMTs
  selected_sets <- intersect(selected_sets, names(geneset_list))
  if (length(selected_sets) == 0L) {
    stop(
      "No overlap between pathways in 'results_df' and gene set names in GMT files. ",
      "Check that '", pathway_col, "' matches the GMT gene set names."
    )
  }
  
  geneset_list <- geneset_list[selected_sets]
  
  # Build Jaccard similarity matrix
  n <- length(geneset_list)
  jaccard_sim <- matrix(
    0,
    nrow = n,
    ncol  = n,
    dimnames = list(names(geneset_list), names(geneset_list))
  )
  
  for (i in seq_along(geneset_list)) {
    for (j in seq_along(geneset_list)) {
      a <- geneset_list[[i]]
      b <- geneset_list[[j]]
      inter <- length(intersect(a, b))
      union <- length(unique(c(a, b)))
      if (union == 0) {
        jaccard_sim[i, j] <- 0
      } else {
        jaccard_sim[i, j] <- inter / union
      }
    }
  }
  
  # Distance matrix for hierarchical clustering
  dist_mat <- stats::as.dist(1 - jaccard_sim)
  
  # Hierarchical clustering using Ward.D2
  hc <- stats::hclust(dist_mat, method = "ward.D2")
  
  # Determine maximum k for silhouette evaluation
  max_k <- max(2L, floor(n / max_k_ratio))
  max_k <- min(max_k, n - 1L) # cannot have k >= n for silhouette
  
  if (max_k < 2L) {
    stop("Not enough gene sets to form at least two clusters.")
  }
  
  ks <- 2:max_k
  sil_scores <- numeric(length(ks))
  
  for (idx in seq_along(ks)) {
    k <- ks[idx]
    cl <- stats::cutree(hc, k = k)
    sil <- cluster::silhouette(cl, dist_mat)
    sil_scores[idx] <- mean(sil[, "sil_width"])
  }
  
  optimal_k <- ks[which.max(sil_scores)]
  
  silhouette_data <- tibble::tibble(
    k = ks,
    silhouette = sil_scores
  )
  
  # Cluster membership at optimal k
  cluster_assignments <- stats::cutree(hc, k = optimal_k)
  clusters_tbl <- tibble::tibble(
    pathway = names(cluster_assignments),
    cluster = as.integer(cluster_assignments)
  )
  
  # Return structured output
  list(
    clusters        = clusters_tbl,
    jaccard_matrix  = jaccard_sim,
    hclust_obj      = hc,
    optimal_k       = optimal_k,
    silhouette_data = silhouette_data
  )
}