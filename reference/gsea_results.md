# Simulated GSEA pathway analysis results for TCGA-LUAD

A simulated data frame representing the output of
[`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md)
for three pairwise comparisons of TCGA-LUAD samples across 40 gene sets
from three MSigDB collections (HALLMARK, KEGG, GO). Gene sets and gene
memberships are derived from
[geneset_list](https://danielgarbozo.github.io/OmicsKit/reference/geneset_list.md).
NES values and FDR are simulated with `set.seed(174)` to produce
realistic enrichment patterns, where ~60% of gene sets per comparison
are significant (FDR \< 0.05).

## Usage

``` r
gsea_results
```

## Format

A data frame with 120 rows (40 gene sets x 3 comparisons) and 15
columns:

- NAME:

  Character. Gene set name, matching the names in
  [geneset_list](https://danielgarbozo.github.io/OmicsKit/reference/geneset_list.md).

- SIZE:

  Integer. Number of genes in the gene set.

- ES:

  Numeric. Enrichment score.

- NES:

  Numeric. Normalized enrichment score.

- NOM p-val:

  Numeric. Nominal p-value.

- FDR:

  Numeric. False discovery rate. Approximately 60% of gene sets per
  comparison have FDR \< 0.05.

- FWER p-val:

  Numeric. Family-wise error rate.

- RANK AT MAX:

  Integer. Gene rank at maximum enrichment score.

- Log10FDR:

  Numeric. `-log10(FDR)`.

- tags:

  Numeric. Fraction of gene set in the leading edge (0-1).

- list:

  Numeric. Fraction of the ranked list used (0-1).

- signal:

  Numeric. Enrichment signal strength (0-1).

- LEADING EDGE:

  Character. Leading edge string in GSEA format (e.g.,
  `"tags=20%, list=35%, signal=15%"`).

- COLLECTION:

  Character. MSigDB collection name: `"HALLMARK"`, `"KEGG"`, or `"GO"`.

- COMPARISON:

  Character. Comparison name: `"TumorVsNormal"`, `"MetastasisVsNormal"`,
  or `"MetastasisVsTumor"`.

## Source

Simulated with `set.seed(174)` in `data-raw/gsea_results.R`. Gene set
names and memberships derived from
[geneset_list](https://danielgarbozo.github.io/OmicsKit/reference/geneset_list.md).
NES values and significance are simulated to reflect realistic GSEA
output patterns.

## Details

This dataset is designed to demonstrate
[`splot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/splot_PA.md),
[`multiplot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/multiplot_PA.md),
[`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md),
[`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md),
and
[`heatmap_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_PA.md)
without requiring external GSEA output files.

## See also

[`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md)
which produces this format from real GSEA output;
[`splot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/splot_PA.md),
[`multiplot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/multiplot_PA.md)
for visualization;
[`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md),
[`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md)
for gene-level annotation;
[geneset_list](https://danielgarbozo.github.io/OmicsKit/reference/geneset_list.md)
for the gene set memberships used here.

## Examples

``` r
data(gsea_results)

# Overview
dim(gsea_results)
#> [1] 120  15
table(gsea_results$COMPARISON)
#> 
#> MetastasisVsNormal  MetastasisVsTumor      TumorVsNormal 
#>                 40                 40                 40 
table(gsea_results$COLLECTION)
#> 
#>       GO HALLMARK     KEGG 
#>       33       42       45 

# How many gene sets are significant per comparison?
tapply(gsea_results$FDR < 0.05, gsea_results$COMPARISON, sum)
#> MetastasisVsNormal  MetastasisVsTumor      TumorVsNormal 
#>                 24                 24                 24 

# Single comparison plot
single <- gsea_results[gsea_results$COMPARISON == "TumorVsNormal", ]
if (FALSE) { # \dontrun{
splot_PA(
  data           = single,
  geneset_col    = "NAME",
  collection_col = "COLLECTION",
  nes_col        = "NES",
  fdr_col        = "FDR"
)
} # }

# Multi-comparison plot
if (FALSE) { # \dontrun{
multiplot_PA(
  data             = gsea_results,
  comparison_col   = "COMPARISON",
  facet_col        = "NAME",
  fdr_col          = "FDR",
  comparison_order = c("TumorVsNormal", "MetastasisVsNormal",
                       "MetastasisVsTumor")
)
} # }
```
