# Merge GSEA result files into a single data frame

Reads all `.tsv` files produced by `GSEA_merge.sh` (from the GSEA.sh
pipeline) from a directory, standardizes numeric columns, parses the
leading edge string, computes `-log10(FDR)`, and returns a single merged
data frame ready for downstream use with
[`splot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/splot_PA.md),
[`multiplot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/multiplot_PA.md)
,
[`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md),
and
[`heatmap_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_PA.md).

## Usage

``` r
merge_PA(input_directory, fdr_replace = 0.001)
```

## Arguments

- input_directory:

  Character. Path to the directory containing one or more GSEA
  collection result files in `.tsv` format (e.g., output of
  `GSEA_merge.sh`). Each file must end in `.tsv`.

- fdr_replace:

  Numeric. Value used to replace `FDR q-val = 0`. This occurs when no
  permutation produced an NES as extreme as the observed one, meaning
  the true FDR is below `1 / n_permutations`. With the standard 1,000
  permutations, the recommended value is `0.001`. Adjust to
  `1 / n_permutations` if a different number of permutations was used.
  Default: `0.001`.

## Value

A data frame (`gsea_data`) containing all merged and processed GSEA
results with standardized column names.

## Details

**Input file format:** Each `.tsv` file corresponds to one MSigDB
collection (e.g., `H.tsv`, `C2.tsv`) and must follow the standard GSEA
output format with the following columns:

- `NAME`: gene set name.

- `SIZE`: number of genes in the gene set.

- `ES`: enrichment score.

- `NES`: normalized enrichment score.

- `NOM p-val`: nominal p-value.

- `FDR q-val`: false discovery rate. Values of exactly `0` indicate that
  no permutation produced an equally extreme NES (i.e., the true FDR is
  below the permutation resolution `1 / n_permutations`). These are
  replaced by `fdr_replace` to avoid `-Inf` in log-transforms.

- `FWER p-val`: family-wise error rate.

- `RANK AT MAX`: gene rank at maximum enrichment score.

- `LEADING EDGE`: string encoding the leading edge statistics in the
  format `"tags=XX%, list=XX%, signal=XX%"`. Parsed into three numeric
  columns: `tags` (fraction of gene set in leading edge), `list`
  (fraction of ranked list used), and `signal` (enrichment signal
  strength).

- `Comparison`: name of the comparison (e.g., `"TumorVsNormal"`).
  Renamed to `COMPARISON` in the output. Required for visualization with
  [`splot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/splot_PA.md)
  or
  [`multiplot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/multiplot_PA.md)
  .

- `GS<br> follow link to MSigDB` and `GS DETAILS`: removed
  automatically.

**Output columns:** All input columns (minus the two removed above)
plus:

- `COLLECTION`: name of the MSigDB collection, derived from the file
  name by removing the `.tsv` suffix.

- `tags`, `list`, `signal`: numeric leading edge components (0-1 scale).

- `Log10FDR`: `-log10(FDR)` computed after applying `fdr_replace`.

- `FDR`: renamed from `FDR q-val`.

- `COMPARISON`: renamed from `Comparison`.

## Note

The input `.tsv` files must contain a `Comparison` column identifying
each comparison (e.g., `"TumorVsNormal"`). This column is renamed to
`COMPARISON` in the output and is required by
[`splot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/splot_PA.md)
or
[`multiplot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/multiplot_PA.md)
to operate in multi-comparison mode. If your files come from a single
comparison and do not have this column, add it manually to each file
before merging: `your_data$Comparison <- "YourComparisonName"`.

## See also

[`splot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/splot_PA.md)
or
[`multiplot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/multiplot_PA.md)
for visualization of merged results;
[`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md)
and
[`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md)
for gene-level annotation;
[`heatmap_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_PA.md)
for leading edge heatmaps;
[`list_gmts()`](https://danielgarbozo.github.io/OmicsKit/reference/list_gmts.md)
to load gene sets.

## Examples

``` r
if (FALSE) { # \dontrun{
# Merge all GSEA collection TSV files from a directory
gsea_data <- merge_PA(
  input_directory = "path/to/gsea_results/",
  fdr_replace     = 0.001   # standard for 1000 permutations
)

# Inspect result
head(gsea_data)
colnames(gsea_data)

# Use directly in downstream functions
gsl        <- list_gmts("path/to/gmt_folder/")
ranked     <- deseq2_results$gene_id[order(deseq2_results$stat,
                                           decreasing = TRUE)]
gene_lists <- getgenesPA(gsea_data, gsl, ranked)
pa_annot   <- addgenesPA(gsea_data, gene_lists)

plot_PA(gsea_data, comparison_col = "COMPARISON")
} # }
```
