# Generate representative super-term labels for gene set communities

For each community in a gene set network, applies **TF-IDF** (Term
Frequency-Inverse Document Frequency) weighting to the words present in
gene set names to produce a short, representative label called a
*super-term*.

## Usage

``` r
get_superterm(
  geneset_names,
  community_membership,
  n_terms = 3,
  remove_prefix = TRUE
)
```

## Arguments

- geneset_names:

  Character vector of gene set names (nodes in the network).

- community_membership:

  A named numeric or integer vector mapping each gene set to its
  community ID. Typically the output of
  [`igraph::membership()`](https://r.igraph.org/reference/communities.html)
  applied to a community detection result (e.g.,
  [`igraph::cluster_louvain()`](https://r.igraph.org/reference/cluster_louvain.html)).
  Must have the same length as `geneset_names`. See
  [`get_network_communities()`](https://danielgarbozo.github.io/OmicsKit/reference/get_network_communities.md)
  for a simpler workflow.

- n_terms:

  Integer. Number of top TF-IDF terms to include in each label. Default:
  `3`.

- remove_prefix:

  Logical. If `TRUE`, removes the text before the first underscore in
  gene set names (e.g., strips the `"KEGG_"` prefix from
  `"KEGG_GLYCOLYSIS"`). Default: `TRUE`.

## Value

A named list with two elements:

- `$mapping`: A
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  with columns `geneset`, `community`, and `superterm` — one row per
  gene set, sorted by community.

- `$summary`: A
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  with columns `community`, `superterm`, and `n_genesets` — one row per
  community, sorted by decreasing size.

## Details

**How TF-IDF works here:** each gene set name is treated as a document
and each word as a term. TF-IDF upweights words that are frequent within
a community but rare across all communities, making the resulting label
specific to that cluster rather than generic. A frequency-based fallback
is used when TF-IDF returns no terms (e.g., very small communities).

Common pathway words (`"pathway"`, `"signaling"`, `"regulation"`, etc.)
and standard English stopwords are removed before scoring.

**Note:** this function is most easily used through
[`get_network_communities()`](https://danielgarbozo.github.io/OmicsKit/reference/get_network_communities.md),
which handles community detection and calls `get_superterm()`
internally. If you prefer to call it directly, you need a community
membership vector:

    adj        <- (jac$jaccard_sim > 0.3) * 1
    g          <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
    comm       <- igraph::cluster_louvain(g)
    membership <- igraph::membership(comm)

    st <- get_superterm(
      geneset_names        = names(membership),
      community_membership = membership
    )

## See also

[`get_network_communities()`](https://danielgarbozo.github.io/OmicsKit/reference/get_network_communities.md),
[`network_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust.md),
[`network_clust_gg()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust_gg.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Recommended: use get_network_communities() which calls this internally
net <- get_network_communities(jac, threshold = 0.3)
net$superterms$mapping
net$superterms$summary

# Direct usage with a pre-computed membership vector
adj        <- (jac$jaccard_sim > 0.3) * 1
g          <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
comm       <- igraph::cluster_louvain(g)
membership <- igraph::membership(comm)

st <- get_superterm(
  geneset_names        = names(membership),
  community_membership = membership,
  n_terms              = 3,
  remove_prefix        = TRUE
)

st$mapping   # per-gene-set labels
st$summary   # per-community summary
} # }
```
