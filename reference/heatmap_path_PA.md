# Plot leading edge heatmaps from GSEA analysis results using file paths

Generates one heatmap per gene set from GSEA/CAMERA/PADOG output by
reading all required inputs from file paths. For each gene set, the
leading edge genes are extracted, ordered by their rank in the ranked
gene list, and plotted as a scaled row heatmap against the expression
matrix.

## Usage

``` r
heatmap_path_PA(
  main_dir = NULL,
  expression_file,
  metadata_file,
  gmt_file,
  ranked_genes_file,
  gsea_file,
  output_dir = "leading_edge_heatmaps",
  sample_col = "Sample",
  group_col = "group",
  save_dataframe = FALSE
)
```

## Arguments

- main_dir:

  Character or `NULL`. Optional base directory prepended to all relative
  file paths. If `NULL` (default), all paths are used as-is.

- expression_file:

  Character. Path to a tab-delimited expression data file. Rows are
  genes (first column or a column named `NAME` used as row names),
  columns are sample IDs. Recommended input: VST-transformed counts.

- metadata_file:

  Character. Path to an Excel (`.xlsx`) metadata file. Must contain a
  column matching `sample_col` (sample IDs) and a column matching
  `group_col` (condition labels, e.g., `"Control"`, `"Treatment"`).

- gmt_file:

  Character. Path to a `.gmt` file defining gene sets. Each row
  contains: gene set name (column 1), description (column 2, ignored),
  and gene symbols (columns 3+).

- ranked_genes_file:

  Character. Path to a tab-delimited file where the first column
  contains gene symbols ordered by their ranking metric (e.g., log2FC or
  signal-to-noise ratio), from most positive to most negative. Used to
  order leading edge genes within each heatmap.

- gsea_file:

  Character. Path to a GSEA results `.tsv` file containing at least the
  columns `NAME`, `SIZE`, and `tags` (from the `LEADING EDGE` column
  parsed by
  [`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md)).

- output_dir:

  Character. Directory where heatmap files are saved. Created
  automatically if it does not exist. Default:
  `"leading_edge_heatmaps"`.

- sample_col:

  Name of the sample ID column in the metadata file. Default:
  `"Sample"`.

- group_col:

  Name of the condition/group column in the metadata file (e.g.,
  `"Control"` vs `"Treatment"`). Used for heatmap column annotations.
  Default: `"group"`.

- save_dataframe:

  Logical. If `TRUE`, saves the intermediate data frame (gene sets with
  computed leading edge genes) as a `.tsv` file in `output_dir` before
  plotting. Useful for inspection or reuse. Default: `FALSE`.

## Value

Invisibly returns `TRUE` upon completion. Saves two files per gene set
in `output_dir`:

- `<geneset_name>_heatmap.pdf`

- `<geneset_name>_heatmap.jpg`

If `save_dataframe = TRUE`, also saves
`<output_dir>/leading_edge_genes_df.tsv`.

## Details

This function is the file-path-based alternative to
[`heatmap_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_PA.md),
which accepts R objects directly. Use this version when working from raw
output files on disk (e.g., directly after running `GSEA_merge.sh`).

## Note

For a more flexible workflow that accepts R objects directly (avoiding
repeated file reads), use
[`heatmap_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_PA.md)
instead, which takes `expression_data`, `metadata`, and `pa_data_annot`
as R objects and integrates with
[`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md)
and
[`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md).

## See also

[`heatmap_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_PA.md)
for the R-object-based alternative;
[`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md)
and
[`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md)
for extracting leading edge genes from R objects;
[`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md)
to generate the GSEA results input;
[`list_gmts()`](https://danielgarbozo.github.io/OmicsKit/reference/list_gmts.md)
to load GMT files as R objects.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run with all files in a single base directory
heatmap_path_PA(
  main_dir          = "path/to/analysis/",
  expression_file   = "vst_expression.tsv",
  metadata_file     = "metadata.xlsx",
  gmt_file          = "genesets.gmt",
  ranked_genes_file = "ranked_genes.tsv",
  gsea_file         = "gsea_results.tsv",
  output_dir        = "leading_edge_heatmaps",
  sample_col        = "Sample",
  group_col         = "group",
  save_dataframe    = TRUE
)
# Saves:
#   leading_edge_heatmaps/<geneset>_heatmap.pdf
#   leading_edge_heatmaps/<geneset>_heatmap.jpg
#   leading_edge_heatmaps/leading_edge_genes_df.tsv  (if save_dataframe = TRUE)

# Run with absolute paths (no main_dir)
heatmap_path_PA(
  expression_file   = "/data/vst_counts.tsv",
  metadata_file     = "/data/metadata.xlsx",
  gmt_file          = "/data/h.all.v2023.gmt",
  ranked_genes_file = "/data/ranked_genes.tsv",
  gsea_file         = "/data/gsea_results.tsv"
)
} # }
```
