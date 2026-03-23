# Detect gene set communities and generate super-term labels

Convenience wrapper that builds a binary adjacency network from a
Jaccard similarity matrix, runs a community-detection algorithm, and
optionally generates super-term labels for each community via
[`get_superterm()`](https://danielgarbozo.github.io/OmicsKit/reference/get_superterm.md).
Designed to be the single step between
[`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md)
and the network plotting functions
[`network_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust.md)
/
[`network_clust_gg()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust_gg.md).

## Usage

``` r
get_network_communities(
  x,
  threshold = 0.3,
  method = "louvain",
  superterms = TRUE,
  n_terms = 3,
  remove_prefix = TRUE,
  seed = 174
)
```

## Arguments

- x:

  A `JaccardResult` object (output of
  [`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md)).

- threshold:

  Numeric between 0 and 1. Gene set pairs with a Jaccard similarity
  above this value are connected in the network. Default: `0.3`.

- method:

  Character. Community detection algorithm to use. One of:

  - `"louvain"` —
    [`igraph::cluster_louvain()`](https://r.igraph.org/reference/cluster_louvain.html):
    fast, recommended for most use cases. Default.

  - `"fast_greedy"` —
    [`igraph::cluster_fast_greedy()`](https://r.igraph.org/reference/cluster_fast_greedy.html):
    optimizes modularity greedily, works well on mid-size networks.

  - `"walktrap"` —
    [`igraph::cluster_walktrap()`](https://r.igraph.org/reference/cluster_walktrap.html):
    random-walk approach, tends to find smaller, tighter communities.

- superterms:

  Logical. If `TRUE`, calls
  [`get_superterm()`](https://danielgarbozo.github.io/OmicsKit/reference/get_superterm.md)
  and includes its output in `$superterms`. Default: `TRUE`.

- n_terms:

  Integer. Number of top TF-IDF terms per super-term label. Passed to
  [`get_superterm()`](https://danielgarbozo.github.io/OmicsKit/reference/get_superterm.md).
  Default: `3`.

- remove_prefix:

  Logical. Remove database prefix before the first underscore (e.g.,
  `"KEGG_"`). Passed to
  [`get_superterm()`](https://danielgarbozo.github.io/OmicsKit/reference/get_superterm.md).
  Default: `TRUE`.

- seed:

  Integer. Random seed for reproducible community detection. Default:
  `174`.

## Value

A named list with four elements:

- `$communities`: The igraph communities object.

- `$membership`: Named integer vector of community IDs, one per gene
  set.

- `$adjacency_matrix`: Binary matrix (`1` if Jaccard \> `threshold`).

- `$superterms`: Output of
  [`get_superterm()`](https://danielgarbozo.github.io/OmicsKit/reference/get_superterm.md)
  with `$mapping` and `$summary`. `NULL` if `superterms = FALSE`.

## See also

[`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md),
[`do_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/do_clust.md),
[`get_superterm()`](https://danielgarbozo.github.io/OmicsKit/reference/get_superterm.md),
[`network_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust.md),
[`network_clust_gg()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust_gg.md)

## Examples

``` r
if (FALSE) { # \dontrun{
gsl <- list_gmts("path/to/gmt_folder/")
res <- read.csv("path/to/results.csv")

# Full workflow
jac   <- geneset_similarity(gsl, res, fdr_th = 0.05)
clust <- do_clust(jac)
net   <- get_network_communities(jac, threshold = 0.3, method = "louvain")

net$membership            # community ID per gene set
net$superterms$mapping    # gene set -> superterm
net$superterms$summary    # community sizes and labels

# Pass results to network plots
plots <- network_clust_gg(
  jac,
  clust_result   = clust,
  superterms     = TRUE,
  superterm_data = net$superterms
)
plots$combined
} # }
```
