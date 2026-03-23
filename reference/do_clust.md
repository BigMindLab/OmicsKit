# Hierarchical clustering of gene sets with silhouette-based optimization

Performs hierarchical clustering on a Jaccard distance matrix, selects
the optimal number of clusters by maximizing average silhouette width,
and returns cluster assignments, a silhouette ggplot2 object, and a
ComplexHeatmap with dendrogram.

## Usage

``` r
do_clust(x, method = "ward.D2", max_k = NULL)
```

## Arguments

- x:

  A `JaccardResult` object (output of
  [`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md))
  or an object of class `dist`.

- method:

  Agglomeration method passed to
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html). Default:
  `"ward.D2"`.

- max_k:

  Maximum number of clusters to evaluate in silhouette analysis.
  Default: `NULL`, which sets it automatically to
  `max(1, floor(n / 2))`.

## Value

A named list with five elements:

- `$hclust`: The
  [`stats::hclust()`](https://rdrr.io/r/stats/hclust.html) object.

- `$cluster_assignments`: A
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  with columns `NAME` and `cluster`.

- `$optimal_k`: Integer. The optimal number of clusters.

- `$silhouette_plot`: A ggplot2 object of average silhouette width vs.
  k.

- `$heatmap`: A `ComplexHeatmap` object. Display with
  `ComplexHeatmap::draw(result$heatmap)`.

## See also

[`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md),
[`get_network_communities()`](https://danielgarbozo.github.io/OmicsKit/reference/get_network_communities.md),
[`network_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust.md),
[`network_clust_gg()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust_gg.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires ComplexHeatmap and cluster packages
geneset_list <- list(
  KEGG_APOPTOSIS      = c("TP53", "BCL2", "CASP3", "BAX"),
  KEGG_CELL_CYCLE     = c("CDK2", "CCND1", "TP53", "RB1"),
  HALLMARK_HYPOXIA    = c("HIF1A", "VEGFA", "LDHA", "BNIP3"),
  HALLMARK_GLYCOLYSIS = c("LDHA", "ENO1", "PKM", "HIF1A"),
  KEGG_P53_PATHWAY    = c("TP53", "MDM2", "CDKN1A", "BAX")
)

results <- data.frame(
  GeneSet = names(geneset_list),
  FDR     = c(0.01, 0.02, 0.03, 0.04, 0.01)
)

jac   <- geneset_similarity(geneset_list, results)
clust <- do_clust(jac)

clust$silhouette_plot               # ggplot2 silhouette curve
ComplexHeatmap::draw(clust$heatmap) # Jaccard heatmap with dendrogram
clust$optimal_k                     # selected number of clusters
clust$cluster_assignments           # tibble: NAME | cluster
} # }
```
