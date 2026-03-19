# =============================================================================
# doclust_PA.R
# Pathway Analysis — Jaccard similarity, hierarchical clustering,
# community detection, and super-term generation.
#
# Functions:
#   calc_jaccard            — Compute Jaccard similarity & distance matrices
#   do_clust                — Hierarchical clustering with silhouette selection
#   get_superterm           — TF-IDF super-term labels for gene set communities
#   get_network_communities — Community detection + super-terms in one call
# =============================================================================

########################
# Function calc_jaccard #
########################

#' Compute Jaccard similarity and distance matrices for gene sets
#'
#' Filters a named list of gene sets by a significance threshold and computes
#' pairwise Jaccard similarity and distance matrices for the retained sets.
#' The output object can be passed directly to [do_clust()],
#' [get_network_communities()], [network_clust()], or [network_clust_gg()],
#' or its individual slots can be used independently (e.g., `$dist_mat` for
#' UMAP, `$jaccard_sim` for custom visualizations).
#'
#' @param geneset_list A named list where each element is a character vector of
#'   gene symbols belonging to that gene set. Typically the output of
#'   [list_gmts()].
#' @param results A data frame with at least two columns: `GeneSet` (gene set
#'   names) and `FDR` (adjusted p-values).
#' @param fdr_th Numeric. FDR cutoff to retain significant gene sets.
#'   Default: `0.05`.
#'
#' @return An object of class `JaccardResult`, a named list with three slots:
#'   * `$jaccard_sim`: Numeric matrix of pairwise Jaccard similarities.
#'   * `$dist_mat`: A `dist` object of 1 - Jaccard similarity, suitable for
#'     clustering or UMAP.
#'   * `$geneset_list`: Named list of gene sets retained after FDR filtering.
#'
#' @examples
#' geneset_list <- list(
#'   KEGG_APOPTOSIS      = c("TP53", "BCL2", "CASP3", "BAX"),
#'   KEGG_CELL_CYCLE     = c("CDK2", "CCND1", "TP53", "RB1"),
#'   HALLMARK_HYPOXIA    = c("HIF1A", "VEGFA", "LDHA", "BNIP3"),
#'   HALLMARK_GLYCOLYSIS = c("LDHA", "ENO1", "PKM", "HIF1A")
#' )
#'
#' results <- data.frame(
#'   GeneSet = names(geneset_list),
#'   FDR     = c(0.01, 0.03, 0.04, 0.20)
#' )
#'
#' # Only the first three gene sets pass the FDR threshold
#' jac <- calc_jaccard(geneset_list, results, fdr_th = 0.05)
#'
#' jac$jaccard_sim   # similarity matrix
#' jac$dist_mat      # distance object (usable in UMAP, clustering, etc.)
#' jac$geneset_list  # filtered gene sets
#'
#' @seealso [list_gmts()], [do_clust()], [get_network_communities()],
#'   [network_clust()], [network_clust_gg()]
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

calc_jaccard <- function(geneset_list, results, fdr_th = 0.05) {

  if (!is.list(geneset_list) || is.null(names(geneset_list))) {
    stop("`geneset_list` must be a named list of character vectors.", call. = FALSE)
  }
  if (!is.data.frame(results)) {
    stop("`results` must be a data frame.", call. = FALSE)
  }
  if (!all(c("GeneSet", "FDR") %in% colnames(results))) {
    stop("`results` must contain columns named `GeneSet` and `FDR`.", call. = FALSE)
  }
  if (!is.numeric(fdr_th) || fdr_th <= 0 || fdr_th > 1) {
    stop("`fdr_th` must be a numeric value between 0 and 1.", call. = FALSE)
  }

  # Filter gene sets by FDR threshold
  selected_sets <- results %>%
    dplyr::filter(.data[["FDR"]] < fdr_th) %>%
    dplyr::pull(.data[["GeneSet"]])

  selected_sets <- intersect(selected_sets, names(geneset_list))

  if (length(selected_sets) == 0) {
    stop(
      "No gene sets remain after FDR filtering. ",
      "Consider increasing `fdr_th`.",
      call. = FALSE
    )
  }

  geneset_list <- geneset_list[selected_sets]
  n            <- length(geneset_list)

  # Build pairwise Jaccard similarity matrix
  jaccard_sim <- matrix(
    0,
    nrow     = n,
    ncol     = n,
    dimnames = list(names(geneset_list), names(geneset_list))
  )

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      a <- geneset_list[[i]]
      b <- geneset_list[[j]]
      u <- length(union(a, b))
      jaccard_sim[i, j] <- if (u == 0) 0 else length(intersect(a, b)) / u
    }
  }

  dist_mat <- stats::as.dist(1 - jaccard_sim)

  result <- list(
    jaccard_sim  = jaccard_sim,
    dist_mat     = dist_mat,
    geneset_list = geneset_list
  )
  class(result) <- "JaccardResult"
  return(result)
}


#####################
# Function do_clust #
#####################

#' Hierarchical clustering of gene sets with silhouette-based optimization
#'
#' Performs hierarchical clustering on a Jaccard distance matrix, selects the
#' optimal number of clusters by maximizing average silhouette width, and
#' returns cluster assignments, a silhouette ggplot2 object, and a
#' ComplexHeatmap with dendrogram.
#'
#' @param x A `JaccardResult` object (output of [calc_jaccard()]) or an
#'   object of class `dist`.
#' @param method Agglomeration method passed to [stats::hclust()].
#'   Default: `"ward.D2"`.
#' @param max_k Maximum number of clusters to evaluate in silhouette analysis.
#'   Default: `NULL`, which sets it automatically to `max(1, floor(n / 2))`.
#'
#' @return A named list with five elements:
#'   * `$hclust`: The [stats::hclust()] object.
#'   * `$cluster_assignments`: A [tibble::tibble()] with columns `NAME` and
#'     `cluster`.
#'   * `$optimal_k`: Integer. The optimal number of clusters.
#'   * `$silhouette_plot`: A ggplot2 object of average silhouette width vs. k.
#'   * `$heatmap`: A `ComplexHeatmap` object. Display with
#'     `ComplexHeatmap::draw(result$heatmap)`.
#'
#' @examples
#' \dontrun{
#' # Requires ComplexHeatmap and cluster packages
#' geneset_list <- list(
#'   KEGG_APOPTOSIS      = c("TP53", "BCL2", "CASP3", "BAX"),
#'   KEGG_CELL_CYCLE     = c("CDK2", "CCND1", "TP53", "RB1"),
#'   HALLMARK_HYPOXIA    = c("HIF1A", "VEGFA", "LDHA", "BNIP3"),
#'   HALLMARK_GLYCOLYSIS = c("LDHA", "ENO1", "PKM", "HIF1A"),
#'   KEGG_P53_PATHWAY    = c("TP53", "MDM2", "CDKN1A", "BAX")
#' )
#'
#' results <- data.frame(
#'   GeneSet = names(geneset_list),
#'   FDR     = c(0.01, 0.02, 0.03, 0.04, 0.01)
#' )
#'
#' jac   <- calc_jaccard(geneset_list, results)
#' clust <- do_clust(jac)
#'
#' clust$silhouette_plot               # ggplot2 silhouette curve
#' ComplexHeatmap::draw(clust$heatmap) # Jaccard heatmap with dendrogram
#' clust$optimal_k                     # selected number of clusters
#' clust$cluster_assignments           # tibble: NAME | cluster
#' }
#'
#' @seealso [calc_jaccard()], [get_network_communities()],
#'   [network_clust()], [network_clust_gg()]
#' @import ggplot2
#' @importFrom rlang .data
#' @export

do_clust <- function(x, method = "ward.D2", max_k = NULL) {

  if (!requireNamespace("cluster",        quietly = TRUE)) {
    stop("Package \"cluster\" must be installed to use this function.",        call. = FALSE)
  }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package \"ComplexHeatmap\" must be installed to use this function.", call. = FALSE)
  }
  if (!requireNamespace("circlize",       quietly = TRUE)) {
    stop("Package \"circlize\" must be installed to use this function.",       call. = FALSE)
  }

  # Extract distance and similarity matrices
  if (inherits(x, "JaccardResult")) {
    dist_mat    <- x$dist_mat
    jaccard_sim <- x$jaccard_sim
  } else if (inherits(x, "dist")) {
    dist_mat    <- x
    jaccard_sim <- 1 - as.matrix(x)
  } else {
    stop(
      "`x` must be a `JaccardResult` object (output of `calc_jaccard()`) ",
      "or an object of class `dist`.",
      call. = FALSE
    )
  }

  n <- attr(dist_mat, "Size")
  if (n < 3) stop("At least 3 gene sets are needed for clustering.", call. = FALSE)

  if (is.null(max_k)) max_k <- max(1L, floor(n / 2L))
  if (max_k >= n)     max_k <- n - 1L
  if (max_k < 2L)     max_k <- 2L

  # Hierarchical clustering
  hc <- stats::hclust(dist_mat, method = method)

  # Silhouette optimisation across k = 2 ... max_k
  sil_scores <- vapply(2:max_k, function(k) {
    cl <- stats::cutree(hc, k = k)
    mean(cluster::silhouette(cl, dist_mat)[, 3])
  }, numeric(1))

  optimal_k <- which.max(sil_scores) + 1L

  # Silhouette ggplot2 object
  df_sil <- tibble::tibble(k = 2:max_k, silhouette = sil_scores)

  p_sil <- ggplot(df_sil, aes(x = .data[["k"]], y = .data[["silhouette"]])) +
    geom_line(color = "gray40") +
    geom_point(color = "gray40") +
    geom_point(
      data  = df_sil[df_sil$k == optimal_k, ],
      color = "red", size = 3
    ) +
    geom_text(
      data  = df_sil[df_sil$k == optimal_k, ],
      aes(label = paste0("k = ", .data[["k"]])),
      vjust = -1, color = "red", size = 4
    ) +
    labs(
      x     = "Number of clusters (k)",
      y     = "Average silhouette width",
      title = "Silhouette analysis"
    ) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust = 0.5))

  # Cluster assignments tibble
  cluster_assignments <- tibble::tibble(
    NAME    = labels(dist_mat),
    cluster = stats::cutree(hc, k = optimal_k)
  )

  # Heatmap with dendrogram
  ht <- ComplexHeatmap::Heatmap(
    as.matrix(jaccard_sim),
    name              = "Jaccard",
    cluster_rows      = hc,
    cluster_columns   = hc,
    show_row_names    = FALSE,
    show_column_names = FALSE,
    col               = circlize::colorRamp2(c(0, 1), c("white", "blue"))
  )

  return(list(
    hclust              = hc,
    cluster_assignments = cluster_assignments,
    optimal_k           = optimal_k,
    silhouette_plot     = p_sil,
    heatmap             = ht
  ))
}


########################
# Function get_superterm #
########################

#' Generate representative super-term labels for gene set communities
#'
#' For each community in a gene set network, applies **TF-IDF** (Term
#' Frequency-Inverse Document Frequency) weighting to the words present in gene
#' set names to produce a short, representative label called a *super-term*.
#'
#' **How TF-IDF works here:** each gene set name is treated as a document and
#' each word as a term. TF-IDF upweights words that are frequent within a
#' community but rare across all communities, making the resulting label
#' specific to that cluster rather than generic. A frequency-based fallback is
#' used when TF-IDF returns no terms (e.g., very small communities).
#'
#' Common pathway words (`"pathway"`, `"signaling"`, `"regulation"`, etc.) and
#' standard English stopwords are removed before scoring.
#'
#' **Note:** this function is most easily used through
#' [get_network_communities()], which handles community detection and calls
#' `get_superterm()` internally. If you prefer to call it directly, you need a
#' community membership vector:
#'
#' ```r
#' adj        <- (jac$jaccard_sim > 0.3) * 1
#' g          <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
#' comm       <- igraph::cluster_louvain(g)
#' membership <- igraph::membership(comm)
#'
#' st <- get_superterm(
#'   geneset_names        = names(membership),
#'   community_membership = membership
#' )
#' ```
#'
#' @param geneset_names Character vector of gene set names (nodes in the
#'   network).
#' @param community_membership A named numeric or integer vector mapping each
#'   gene set to its community ID. Typically the output of
#'   [igraph::membership()] applied to a community detection result (e.g.,
#'   [igraph::cluster_louvain()]). Must have the same length as
#'   `geneset_names`. See [get_network_communities()] for a simpler workflow.
#' @param n_terms Integer. Number of top TF-IDF terms to include in each label.
#'   Default: `3`.
#' @param remove_prefix Logical. If `TRUE`, removes the text before the first
#'   underscore in gene set names (e.g., strips the `"KEGG_"` prefix from
#'   `"KEGG_GLYCOLYSIS"`). Default: `TRUE`.
#'
#' @return A named list with two elements:
#'   * `$mapping`: A [tibble::tibble()] with columns `geneset`, `community`,
#'     and `superterm` — one row per gene set, sorted by community.
#'   * `$summary`: A [tibble::tibble()] with columns `community`, `superterm`,
#'     and `n_genesets` — one row per community, sorted by decreasing size.
#'
#' @examples
#' \dontrun{
#' # Recommended: use get_network_communities() which calls this internally
#' net <- get_network_communities(jac, threshold = 0.3)
#' net$superterms$mapping
#' net$superterms$summary
#'
#' # Direct usage with a pre-computed membership vector
#' adj        <- (jac$jaccard_sim > 0.3) * 1
#' g          <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
#' comm       <- igraph::cluster_louvain(g)
#' membership <- igraph::membership(comm)
#'
#' st <- get_superterm(
#'   geneset_names        = names(membership),
#'   community_membership = membership,
#'   n_terms              = 3,
#'   remove_prefix        = TRUE
#' )
#'
#' st$mapping   # per-gene-set labels
#' st$summary   # per-community summary
#' }
#'
#' @seealso [get_network_communities()], [network_clust()], [network_clust_gg()]
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

get_superterm <- function(geneset_names, community_membership,
                          n_terms = 3, remove_prefix = TRUE) {

  if (!requireNamespace("tm", quietly = TRUE)) {
    stop("Package \"tm\" must be installed to use this function.", call. = FALSE)
  }
  if (!is.character(geneset_names) || length(geneset_names) == 0) {
    stop("`geneset_names` must be a non-empty character vector.", call. = FALSE)
  }
  if (!is.numeric(community_membership) && !is.integer(community_membership)) {
    stop("`community_membership` must be a numeric or integer vector.", call. = FALSE)
  }
  if (length(geneset_names) != length(community_membership)) {
    stop(
      "`geneset_names` and `community_membership` must have the same length.",
      call. = FALSE
    )
  }

  # Pathway-specific stopwords removed before TF-IDF scoring
  pathway_stopwords <- c(
    "pathway", "pathways", "signaling", "signal", "regulation",
    "positive", "negative", "activity", "process", "involved",
    "via", "mediated", "induced", "related", "associated",
    "dependent", "independent"
  )

  unique_comms <- sort(unique(community_membership))

  superterms <- vapply(unique_comms, function(comm_id) {

    nodes <- geneset_names[community_membership == comm_id]

    if (length(nodes) == 0) return("Empty_Cluster")

    # Single gene set: split name into words directly
    if (length(nodes) == 1) {
      clean <- if (remove_prefix) sub("^[^_]+_", "", nodes) else nodes
      words <- strsplit(clean, "[_[:space:]]+")[[1]]
      return(paste(utils::head(words, 3), collapse = "_"))
    }

    # Multiple gene sets: TF-IDF pipeline
    cleaned <- if (remove_prefix) gsub("^[^_]+_", "", nodes) else nodes
    cleaned <- tolower(gsub("_", " ", cleaned))

    corpus <- tm::Corpus(tm::VectorSource(cleaned))
    corpus <- tm::tm_map(corpus, tm::content_transformer(tolower))
    corpus <- tm::tm_map(corpus, tm::removePunctuation)
    corpus <- tm::tm_map(corpus, tm::removeNumbers)
    corpus <- tm::tm_map(corpus, tm::removeWords,
                         c(tm::stopwords("english"), pathway_stopwords))
    corpus <- tm::tm_map(corpus, tm::stripWhitespace)

    dtm       <- tm::DocumentTermMatrix(corpus,
                                        control = list(weighting = tm::weightTfIdf))
    term_freq <- colSums(as.matrix(dtm))
    top_terms <- utils::head(sort(term_freq, decreasing = TRUE), n_terms)

    # Fallback: frequency-based when TF-IDF returns nothing
    if (length(top_terms) == 0) {
      all_words <- unlist(strsplit(cleaned, "\\s+"))
      all_words <- all_words[nchar(all_words) > 3]
      if (length(all_words) == 0) return("Cluster")
      word_freq      <- sort(table(all_words), decreasing = TRUE)
      top_term_names <- names(utils::head(word_freq, n_terms))
    } else {
      top_term_names <- names(top_terms)
    }

    # Capitalize first letter of each term
    top_term_names <- vapply(top_term_names, function(x) {
      paste0(toupper(substring(x, 1, 1)), substring(x, 2))
    }, character(1))

    paste(top_term_names, collapse = "/")

  }, character(1))

  names(superterms) <- as.character(unique_comms)

  mapping <- tibble::tibble(
    geneset   = geneset_names,
    community = as.integer(community_membership),
    superterm = superterms[as.character(community_membership)]
  ) %>%
    dplyr::arrange(.data[["community"]], .data[["geneset"]])

  summary_tbl <- mapping %>%
    dplyr::group_by(.data[["community"]], .data[["superterm"]]) %>%
    dplyr::summarise(n_genesets = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data[["n_genesets"]]))

  return(list(
    mapping = mapping,
    summary = summary_tbl
  ))
}


################################
# Function get_network_communities #
################################

#' Detect gene set communities and generate super-term labels
#'
#' Convenience wrapper that builds a binary adjacency network from a Jaccard
#' similarity matrix, runs a community-detection algorithm, and optionally
#' generates super-term labels for each community via [get_superterm()].
#' Designed to be the single step between [calc_jaccard()] and the network
#' plotting functions [network_clust()] / [network_clust_gg()].
#'
#' @param x A `JaccardResult` object (output of [calc_jaccard()]).
#' @param threshold Numeric between 0 and 1. Gene set pairs with a Jaccard
#'   similarity above this value are connected in the network. Default: `0.3`.
#' @param method Character. Community detection algorithm to use. One of:
#'   * `"louvain"` — [igraph::cluster_louvain()]: fast, recommended for most
#'     use cases. Default.
#'   * `"fast_greedy"` — [igraph::cluster_fast_greedy()]: optimizes modularity
#'     greedily, works well on mid-size networks.
#'   * `"walktrap"` — [igraph::cluster_walktrap()]: random-walk approach,
#'     tends to find smaller, tighter communities.
#' @param superterms Logical. If `TRUE`, calls [get_superterm()] and includes
#'   its output in `$superterms`. Default: `TRUE`.
#' @param n_terms Integer. Number of top TF-IDF terms per super-term label.
#'   Passed to [get_superterm()]. Default: `3`.
#' @param remove_prefix Logical. Remove database prefix before the first
#'   underscore (e.g., `"KEGG_"`). Passed to [get_superterm()]. Default: `TRUE`.
#' @param seed Integer. Random seed for reproducible community detection.
#'   Default: `174`.
#'
#' @return A named list with four elements:
#'   * `$communities`: The igraph communities object.
#'   * `$membership`: Named integer vector of community IDs, one per gene set.
#'   * `$adjacency_matrix`: Binary matrix (`1` if Jaccard > `threshold`).
#'   * `$superterms`: Output of [get_superterm()] with `$mapping` and
#'     `$summary`. `NULL` if `superterms = FALSE`.
#'
#' @examples
#' \dontrun{
#' gsl <- list_gmts("path/to/gmt_folder/")
#' res <- read.csv("path/to/results.csv")
#'
#' # Full workflow
#' jac   <- calc_jaccard(gsl, res, fdr_th = 0.05)
#' clust <- do_clust(jac)
#' net   <- get_network_communities(jac, threshold = 0.3, method = "louvain")
#'
#' net$membership            # community ID per gene set
#' net$superterms$mapping    # gene set -> superterm
#' net$superterms$summary    # community sizes and labels
#'
#' # Pass results to network plots
#' plots <- network_clust_gg(
#'   jac,
#'   clust_result   = clust,
#'   superterms     = TRUE,
#'   superterm_data = net$superterms
#' )
#' plots$combined
#' }
#'
#' @seealso [calc_jaccard()], [do_clust()], [get_superterm()],
#'   [network_clust()], [network_clust_gg()]
#' @importFrom magrittr %>%
#' @export

get_network_communities <- function(x,
                                    threshold     = 0.3,
                                    method        = "louvain",
                                    superterms    = TRUE,
                                    n_terms       = 3,
                                    remove_prefix = TRUE,
                                    seed          = 174) {

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package \"igraph\" must be installed to use this function.", call. = FALSE)
  }
  if (!inherits(x, "JaccardResult")) {
    stop(
      "`x` must be a `JaccardResult` object (output of `calc_jaccard()`).",
      call. = FALSE
    )
  }
  if (!is.numeric(threshold) || threshold <= 0 || threshold >= 1) {
    stop("`threshold` must be a numeric value strictly between 0 and 1.", call. = FALSE)
  }

  method <- match.arg(method, c("louvain", "fast_greedy", "walktrap"))

  # Build binary adjacency matrix
  adjacency_matrix        <- (x$jaccard_sim > threshold) * 1L
  diag(adjacency_matrix)  <- 0L

  g <- igraph::graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")

  if (igraph::vcount(g) == 0) {
    stop("The network has no nodes. Check that `x` contains gene sets.", call. = FALSE)
  }

  # Community detection
  set.seed(seed)
  communities <- switch(
    method,
    louvain     = igraph::cluster_louvain(g),
    fast_greedy = igraph::cluster_fast_greedy(g),
    walktrap    = igraph::cluster_walktrap(g)
  )

  membership <- igraph::membership(communities)

  # Super-terms
  st_result <- NULL
  if (superterms) {
    if (!requireNamespace("tm", quietly = TRUE)) {
      warning(
        "Package \"tm\" is required to generate super-terms. ",
        "Install it with install.packages(\"tm\") or set `superterms = FALSE`.",
        call. = FALSE
      )
    } else {
      st_result <- get_superterm(
        geneset_names        = names(membership),
        community_membership = membership,
        n_terms              = n_terms,
        remove_prefix        = remove_prefix
      )
    }
  }

  return(list(
    communities      = communities,
    membership       = membership,
    adjacency_matrix = adjacency_matrix,
    superterms       = st_result
  ))
}
