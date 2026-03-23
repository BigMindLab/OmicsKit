# Example gene set list for pathway analysis clustering

A named list of 40 curated gene sets spanning four biological themes:
apoptosis & cell death, cell cycle & DNA damage, immune response &
inflammation, and metabolism. Gene set names follow standard database
conventions (`KEGG_`, `HALLMARK_`, `GO_`) and gene symbols are real
human genes. Designed to be used as input to
[`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md).

## Usage

``` r
geneset_list
```

## Format

A named list of 40 elements. Each element is a character vector of human
gene symbols (HGNC) belonging to that gene set. Gene set sizes range
from 11 to 20 genes.

## Source

Curated manually for OmicsKit examples, based on KEGG, MSigDB Hallmark,
and Gene Ontology gene set collections.

## See also

[`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md),
[camera_results](https://danielgarbozo.github.io/OmicsKit/reference/camera_results.md)

## Examples

``` r
data(geneset_list)

# How many gene sets?
length(geneset_list)
#> [1] 40

# Inspect one gene set
geneset_list[["KEGG_APOPTOSIS"]]
#>  [1] "TP53"      "BCL2"      "BCL2L1"    "BAX"       "BAD"       "BID"      
#>  [7] "CASP3"     "CASP8"     "CASP9"     "CYCS"      "APAF1"     "FADD"     
#> [13] "FAS"       "TNFRSF10A" "TNFRSF10B" "MCL1"      "XIAP"      "DIABLO"   

# Use with geneset_similarity()
data(camera_results)
jac <- geneset_similarity(geneset_list, camera_results, fdr_th = 0.05)
```
