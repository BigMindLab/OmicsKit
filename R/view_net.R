###############################
# Plot enrichment clustering #
###############################

#' Plot a gene set network based on clustering output
#'
#' This function takes the output of \code{\link{build_net()}},
#' constructs a gene set similarity network using Jaccard similarity as edge
#' weights, and generates both a static igraph object and an interactive
#' \code{visNetwork} visualization.
#'
#' @param clustering_output A named list as returned by
#'   \code{\link{build_net}}, containing at least
#'   \code{clusters} (a tibble with columns \code{pathway} and \code{cluster})
#'   and \code{jaccard_matrix} (a symmetric numeric matrix).
#' @param edge_threshold A numeric scalar specifying the minimum Jaccard
#'   similarity required for an edge to be drawn between two gene sets.
#'   All similarities below this threshold are set to zero before network
#'   construction. Default is \code{0.25}.
#' @param save_html_path Optional character scalar. If non-\code{NULL}, the
#'   interactive \code{visNetwork} object will be saved as an HTML file at the
#'   specified path using \code{visNetwork::visSave}. If \code{NULL}, no HTML
#'   file is written. Default is \code{NULL}.
#'
#' @details
#' The static graph is built with \pkg{igraph} from a thresholded Jaccard
#' similarity matrix. Nodes represent gene sets (pathways), and edges represent
#' Jaccard similarity greater than or equal to \code{edge_threshold}.
#'
#' Cluster memberships are taken from \code{clustering_output$clusters} and
#' mapped to the node attribute \code{group} in the interactive network, which
#' can be used by \pkg{visNetwork} for coloring.
#'
#' This function does not perform any plotting side effects by default; it
#' returns objects that can be plotted or further customized by the caller.
#'
#' @return A named list with the following elements:
#' \itemize{
#'   \item \code{static_graph}: An \code{\link[igraph]{igraph}} object built
#'     from the thresholded Jaccard similarity matrix.
#'   \item \code{interactive_plot}: A \code{\link[visNetwork]{visNetworkProxy}}
#'     or \code{visNetwork} object representing the interactive network
#'     visualization.
#' }
#'
#' @examples
#' \dontrun{
#' clustering <- build_net(
#'   results_df   = res,
#'   gmt_path     = "gmt_dir/",
#'   fdr_threshold = 0.25,
#'   fdr_col       = "FDR_gsea",
#'   pathway_col   = "Pathway"
#' )
#'
#' plots <- view_net(
#'   clustering_output = clustering,
#'   edge_threshold    = 0.25,
#'   save_html_path    = "network.html"
#' )
#'
#' # Static plot
#' plot(plots$static_graph,
#'      vertex.label.cex = 0.7,
#'      edge.width       = igraph::E(plots$static_graph)$weight * 5)
#'
#' # Interactive plot (in RStudio viewer or browser)
#' plots$interactive_plot
#' }
#'
#' @importFrom igraph graph_from_adjacency_matrix V E as_data_frame
#' @importFrom visNetwork visNetwork visOptions visLayout visSave
#' @importFrom dplyr left_join rename
#' @importFrom tibble tibble
#' @export
view_net <- function(clustering_output,
                                 edge_threshold = 0.25,
                                 save_html_path = NULL) {
  # Basic input checks
  if (!is.list(clustering_output)) {
    stop("'clustering_output' must be a list as returned by 'calculate_geneset_clusters'.")
  }
  required_elements <- c("clusters", "jaccard_matrix")
  missing_elements <- setdiff(required_elements, names(clustering_output))
  if (length(missing_elements) > 0L) {
    stop(
      "'clustering_output' is missing required elements: ",
      paste(missing_elements, collapse = ", ")
    )
  }
  
  clusters_tbl <- clustering_output$clusters
  jaccard_matrix <- clustering_output$jaccard_matrix
  
  if (!is.matrix(jaccard_matrix) || !is.numeric(jaccard_matrix)) {
    stop("'clustering_output$jaccard_matrix' must be a numeric matrix.")
  }
  if (is.null(rownames(jaccard_matrix)) || is.null(colnames(jaccard_matrix))) {
    stop("The Jaccard matrix must have row and column names corresponding to pathways.")
  }
  if (!is.data.frame(clusters_tbl)) {
    stop("'clustering_output$clusters' must be a data.frame or tibble.")
  }
  if (!all(c("pathway", "cluster") %in% colnames(clusters_tbl))) {
    stop("'clusters' must contain 'pathway' and 'cluster' columns.")
  }
  if (!is.numeric(edge_threshold) || length(edge_threshold) != 1L || is.na(edge_threshold)) {
    stop("'edge_threshold' must be a single numeric value.")
  }
  
  # Threshold the Jaccard matrix to define edges
  adj <- jaccard_matrix
  adj[adj < edge_threshold] <- 0
  diag(adj) <- 0
  
  # Build igraph object
  g <- igraph::graph_from_adjacency_matrix(
    adjmatrix = adj,
    mode      = "undirected",
    weighted  = TRUE,
    diag      = FALSE
  )
  
  # Static graph (igraph object)
  static_graph <- g
  
  # Prepare nodes and edges for visNetwork
  pathways <- rownames(adj)
  
  nodes <- tibble::tibble(
    id    = pathways,
    label = pathways
  )
  
  # Merge cluster membership
  clusters_tbl <- dplyr::distinct(clusters_tbl, .data$pathway, .data$cluster)
  nodes <- dplyr::left_join(
    nodes,
    dplyr::rename(clusters_tbl, id = .data$pathway, group = .data$cluster),
    by = "id"
  )
  
  # Convert edges to data.frame
  edges_df <- igraph::as_data_frame(g, what = "edges")
  # visNetwork uses columns 'from' and 'to'
  colnames(edges_df)[1:2] <- c("from", "to")
  
  interactive_plot <- visNetwork::visNetwork(nodes, edges_df) |>
    visNetwork::visOptions(highlightNearest = TRUE) |>
    visNetwork::visLayout(randomSeed = 174)
  
  # Save HTML if requested
  if (!is.null(save_html_path)) {
    if (!is.character(save_html_path) || length(save_html_path) != 1L) {
      stop("'save_html_path' must be a single character string if provided.")
    }
    visNetwork::visSave(interactive_plot, file = save_html_path)
  }
  
  list(
    static_graph     = static_graph,
    interactive_plot = interactive_plot
  )
}