# Save the results of a Differential Expression analysis

This function takes as input the output of the function "results()" of
DEseq2. And will save 3 tables:

- A table with all genes

- A table including only the over-expressed genes

- A table including only the under-expressed genes

## Usage

``` r
save_results(df, name, l2fc = 0, cutoff_alpha = 0.25)
```

## Arguments

- df:

  A dataframe with the results of a Differential Expression analysis.

- name:

  The name to be used to save the tables, without file extension.

- l2fc:

  The cut-off of Log2(Fold Change) for the over- and under-expressed
  tables. Default = 0.

- cutoff_alpha:

  The cut-off of the False Discovery Rate (FDR o padj). Default = 0.25.

## Value

Invisibly returns `NULL`. Saves three `.xlsx` files to the working
directory:

- `<name>_full.xlsx` : all genes.

- `<name>_up_log2FC><l2fc>_FDR<<cutoff_alpha>.xlsx` : upregulated genes.

- `<name>_down_log2FC<<l2fc>_FDR<<cutoff_alpha>.xlsx` : downregulated
  genes.

## See also

[`detect_filter()`](https://danielgarbozo.github.io/OmicsKit/reference/detect_filter.md)
to further filter saved results;
[deseq2_results](https://danielgarbozo.github.io/OmicsKit/reference/deseq2_results.md)
for an example input.

## Examples

``` r
if (FALSE) { # \dontrun{
data(deseq2_results)

# Save full results + over/under-expressed tables as .xlsx files
save_results(
  df           = deseq2_results,
  name         = "TCGA_LUAD_TumorVsNormal",
  l2fc         = 1,
  cutoff_alpha = 0.05
)

# Creates:
#   TCGA_LUAD_TumorVsNormal_full.xlsx
#   TCGA_LUAD_TumorVsNormal_up_log2FC>1_FDR<0.05.xlsx
#   TCGA_LUAD_TumorVsNormal_down_log2FC<1_FDR<0.05.xlsx
} # }
```
