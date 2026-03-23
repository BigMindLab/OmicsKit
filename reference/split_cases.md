# Obtain 10 exclusive cases from 3 comparisons.

When performing differential expression analysis of a study with 3
phenotypes, including the baseline, there are 10 mutually exclusive
cases where genes can fall into. This function allows us to obtain these
10 cases and saves them into a list.

## Usage

``` r
split_cases(
  df.BvsA = NULL,
  df.CvsA = NULL,
  df.BvsC = NULL,
  unique_id = "ensembl",
  significance_var = "padj",
  significance_cutoff = 0.25,
  change_var = "log2FoldChange",
  change_cutoff = 0
)
```

## Arguments

- df.BvsA:

  Data frame comparing the first condition to the baseline.

- df.CvsA:

  Data frame comparing the second condition to the baseline.

- df.BvsC:

  Data frame comparing the two conditions of the study.

- unique_id:

  Column name with the unique identifiers in the data tables (i.e.
  ensembl for DESeq2, set NAME for GSEA).

- significance_var:

  Variable of significance to filter by (padj for DESeq2 or FDR for
  GSEA).

- significance_cutoff:

  Cut-off of the significance variable.

- change_var:

  Variable that indicates the direction of the change (i.e.
  log2FoldChange in DESeq2, NES in GSEA).

- change_cutoff:

  The values of the change variable will be filtered by \|change_var\|
  \> change_cutoff. Default: 0.

## Value

A named list of 10 data frames (`$Case1` through `$Case10`), each
containing the genes belonging to that mutually exclusive expression
pattern. Cases 1–6 contain a `trend` column (`"up"` or `"dn"`). Case 10
contains genes not significant in any comparison.

## Details

The 10 cases are:

- **Case 1** : Ladder: significant in all 3, same direction.

- **Case 2** : Stronger in condition 1: significant in all 3, direction
  reverses between conditions.

- **Case 3** : Stronger in condition 2.

- **Case 4** : Marker of condition 2: significant in CvsA and BvsC only.

- **Case 5** : Marker of condition 1: significant in BvsA and BvsC only.

- **Case 6** : Shared: significant in BvsA and CvsA only.

- **Cases 7-9** : Significant in only one comparison.

- **Case 10** : Not significant in any comparison.

## See also

[`detect_filter()`](https://danielgarbozo.github.io/OmicsKit/reference/detect_filter.md)
to pre-filter genes before splitting;
[`nice_Volcano()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_Volcano.md)
to visualize individual comparison results.

## Examples

``` r
if (FALSE) { # \dontrun{
# split_cases requires three DESeq2 comparisons.
# Simulate a 3-phenotype study: Normal (A), Primary (B), Metastasis (C)
set.seed(174)
n_genes <- 500

make_res <- function(seed) {
  set.seed(seed)
  data.frame(
    ensembl        = paste0("ENSG", sprintf("%011d", seq_len(n_genes))),
    log2FoldChange = rnorm(n_genes, 0, 1.5),
    padj           = runif(n_genes, 0, 0.5),
    stringsAsFactors = FALSE
  )
}

df_BvsA <- make_res(1)
df_CvsA <- make_res(2)
df_BvsC <- make_res(3)

cases <- split_cases(
  df.BvsA              = df_BvsA,
  df.CvsA              = df_CvsA,
  df.BvsC              = df_BvsC,
  unique_id            = "ensembl",
  significance_var     = "padj",
  significance_cutoff  = 0.25,
  change_var           = "log2FoldChange",
  change_cutoff        = 0
)

# Number of genes per case
sapply(cases, nrow)

# Inspect Case 1 (ladder genes: significant in all 3 comparisons)
head(cases$Case1)
} # }
```
