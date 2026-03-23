# Obtain annotations to display significance.

This function will create asteriscs (\*) from DESeq2 results objects to
represent the significance of comparisons.

## Usage

``` r
get_stars(geneID, object, thresholds = c(0.001, 0.01, 0.1, 0.25))
```

## Arguments

- geneID:

  Ensembl ID of the gene of interest.

- object:

  DESeq2 results object of a comparison.

- thresholds:

  Vector with 4 values of significance. Default c(0.001, 0.01, 0.1,
  0.25).

## Value

A single character string: `"****"`, `"***"`, `"**"`, `"*"`, `"ns"` (not
significant), or `"Gene ID not found"` if the gene is absent from
`object`.

## See also

[`detect_filter()`](https://danielgarbozo.github.io/OmicsKit/reference/detect_filter.md)
to identify detectable genes before annotating;
[`nice_VSB()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_VSB.md)
where significance stars can be added to plots;
[deseq2_results](https://danielgarbozo.github.io/OmicsKit/reference/deseq2_results.md)
for an example input.

## Examples

``` r
data(deseq2_results)

# get_stars expects a column named "ensembl"
res <- deseq2_results
colnames(res)[colnames(res) == "gene_id"] <- "ensembl"

# Get significance stars for the most significant gene
get_stars(
  geneID = res$ensembl[1],
  object = res
)
#> [1] "****"

# Custom thresholds
get_stars(
  geneID     = res$ensembl[1],
  object     = res,
  thresholds = c(0.001, 0.01, 0.05, 0.10)
)
#> [1] "****"

# Non-significant gene
get_stars(
  geneID = res$ensembl[nrow(res)],
  object = res
)
#> [1] "ns"
```
