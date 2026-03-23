# Example CAMERA enrichment results for pathway analysis clustering

A data frame simulating the output of a CAMERA differential expression
analysis, containing significance values for the 40 gene sets in
[geneset_list](https://danielgarbozo.github.io/OmicsKit/reference/geneset_list.md).
Approximately 60% of gene sets have FDR \< 0.05, providing enough
significant sets for meaningful clustering. Designed to be used
alongside
[geneset_list](https://danielgarbozo.github.io/OmicsKit/reference/geneset_list.md)
as input to
[`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md).

## Usage

``` r
camera_results
```

## Format

A data frame with 40 rows and 4 columns:

- GeneSet:

  Character. Gene set name, matching the names in
  [geneset_list](https://danielgarbozo.github.io/OmicsKit/reference/geneset_list.md).

- Direction:

  Character. Enrichment direction: `"Up"` or `"Down"`.

- PValue:

  Numeric. Raw p-value from the simulated CAMERA test.

- FDR:

  Numeric. Benjamini-Hochberg adjusted p-value.

## Source

Simulated with `set.seed(1905)` in `data-raw/example_PA.R` for OmicsKit
examples.

## See also

[`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md),
[geneset_list](https://danielgarbozo.github.io/OmicsKit/reference/geneset_list.md)

## Examples

``` r
data(camera_results)

# Overview
head(camera_results)
#>                  GeneSet Direction      PValue        FDR
#> 1         KEGG_APOPTOSIS        Up 0.017401305 0.04321579
#> 2     HALLMARK_APOPTOSIS        Up 0.012585092 0.04321579
#> 3 GO_INTRINSIC_APOPTOSIS      Down 0.010449980 0.04321579
#> 4 GO_EXTRINSIC_APOPTOSIS        Up 0.007223143 0.04321579
#> 5   HALLMARK_P53_PATHWAY        Up 0.016167120 0.04321579
#> 6     KEGG_P53_SIGNALING        Up 0.006856107 0.04321579

# How many gene sets are significant?
sum(camera_results$FDR < 0.05)
#> [1] 24

# Use with geneset_similarity()
data(geneset_list)
jac <- geneset_similarity(geneset_list, camera_results, fdr_th = 0.05)
```
