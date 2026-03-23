# Compute Jaccard similarity and distance matrices for gene sets

Filters a named list of gene sets by a significance threshold and
computes pairwise Jaccard similarity and distance matrices for the
retained sets. The output object can be passed directly to
[`do_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/do_clust.md),
[`get_network_communities()`](https://danielgarbozo.github.io/OmicsKit/reference/get_network_communities.md),
[`network_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust.md),
or
[`network_clust_gg()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust_gg.md),
or its individual slots can be used independently (e.g., `$dist_mat` for
UMAP, `$jaccard_sim` for custom visualizations).

## Usage

``` r
geneset_similarity(geneset_list, results, fdr_th = 0.05)
```

## Arguments

- geneset_list:

  A named list where each element is a character vector of gene symbols
  belonging to that gene set. Typically the output of
  [`list_gmts()`](https://danielgarbozo.github.io/OmicsKit/reference/list_gmts.md).

- results:

  A data frame with at least two columns: `GeneSet` (gene set names) and
  `FDR` (adjusted p-values).

- fdr_th:

  Numeric. FDR cutoff to retain significant gene sets. Default: `0.05`.

## Value

An object of class `JaccardResult`, a named list with three slots:

- `$jaccard_sim`: Numeric matrix of pairwise Jaccard similarities.

- `$dist_mat`: A `dist` object of 1 - Jaccard similarity, suitable for
  clustering or UMAP.

- `$geneset_list`: Named list of gene sets retained after FDR filtering.

## See also

[`list_gmts()`](https://danielgarbozo.github.io/OmicsKit/reference/list_gmts.md),
[`do_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/do_clust.md),
[`get_network_communities()`](https://danielgarbozo.github.io/OmicsKit/reference/get_network_communities.md),
[`network_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust.md),
[`network_clust_gg()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust_gg.md)

## Examples

``` r
geneset_list <- list(
  KEGG_APOPTOSIS      = c("TP53", "BCL2", "CASP3", "BAX"),
  KEGG_CELL_CYCLE     = c("CDK2", "CCND1", "TP53", "RB1"),
  HALLMARK_HYPOXIA    = c("HIF1A", "VEGFA", "LDHA", "BNIP3"),
  HALLMARK_GLYCOLYSIS = c("LDHA", "ENO1", "PKM", "HIF1A")
)

results <- data.frame(
  GeneSet = names(geneset_list),
  FDR     = c(0.01, 0.03, 0.04, 0.20)
)

# Only the first three gene sets pass the FDR threshold
jac <- geneset_similarity(geneset_list, results, fdr_th = 0.05)

jac$jaccard_sim   # similarity matrix
#>                  KEGG_APOPTOSIS KEGG_CELL_CYCLE HALLMARK_HYPOXIA
#> KEGG_APOPTOSIS        1.0000000       0.1428571                0
#> KEGG_CELL_CYCLE       0.1428571       1.0000000                0
#> HALLMARK_HYPOXIA      0.0000000       0.0000000                1
jac$dist_mat      # distance object (usable in UMAP, clustering, etc.)
#>                  KEGG_APOPTOSIS KEGG_CELL_CYCLE
#> KEGG_CELL_CYCLE       0.8571429                
#> HALLMARK_HYPOXIA      1.0000000       1.0000000
jac$geneset_list  # filtered gene sets
#> $KEGG_APOPTOSIS
#> [1] "TP53"  "BCL2"  "CASP3" "BAX"  
#> 
#> $KEGG_CELL_CYCLE
#> [1] "CDK2"  "CCND1" "TP53"  "RB1"  
#> 
#> $HALLMARK_HYPOXIA
#> [1] "HIF1A" "VEGFA" "LDHA"  "BNIP3"
#> 
```
