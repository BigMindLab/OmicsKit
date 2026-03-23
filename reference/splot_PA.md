# Pathway analysis visualization for a single comparison

Generates a publication-quality multi-panel pathway enrichment plot for
a single comparison using patchwork. Gene sets appear on the y-axis
grouped by MSigDB collection, NES on the x-axis, and -log10(FDR) as fill
color. Six panels are assembled side by side: a "Pathways" label, gene
set names, the NES bar chart, collection labels, a "MSigDB" label, and
the color legend.

## Usage

``` r
splot_PA(
  data,
  geneset_col = "NAME",
  collection_col = "COLLECTION",
  nes_col = "NES",
  fdr_col = "FDR",
  order = "desc",
  fill_limits = NULL,
  fill_palette = c("white", "red"),
  theme_params = list()
)
```

## Arguments

- data:

  A data frame of pathway analysis results for a single comparison.
  Typically the output of
  [`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md)
  filtered to one value of the `COMPARISON` column, or results from a
  single CAMERA/GSEA run. Must contain the columns specified by
  `geneset_col`, `collection_col`, `nes_col`, and `fdr_col`.

- geneset_col:

  Name of the column containing gene set labels shown on the y-axis.
  Default: `"NAME"`.

- collection_col:

  Name of the column containing MSigDB collection labels used to group
  gene sets (e.g., `"KEGG"`, `"HALLMARK"`, `"GO"`). Default:
  `"COLLECTION"`.

- nes_col:

  Name of the column containing NES values (x-axis). Default: `"NES"`.

- fdr_col:

  Name of the column containing FDR values. `-log10(FDR)` is computed
  internally and used as the fill color. Default: `"FDR"`.

- order:

  One of `"desc"` or `"asc"`. Sort order for NES values on the y-axis.
  Default: `"desc"`.

- fill_limits:

  Numeric vector of length 2 setting the color scale range for
  `-log10(FDR)`. Values outside this range are clamped to the nearest
  limit. For example, `fill_limits = c(0, 5)` maps all gene sets with
  `-log10(FDR) >= 5` (i.e., FDR \<= 0.00001) to the maximum color (red),
  and any value below 0 to the minimum color (white). Useful when a few
  gene sets have extreme significance that washes out color variation in
  the rest. Default: `NULL` (auto uses the actual data range).

- fill_palette:

  Character vector of two colors for the fill gradient (low to high
  -log10(FDR)). Default: `c("white", "red")`.

- theme_params:

  Named list to override default theme parameters. See Details.

## Value

A `patchwork` object combining six ggplot2 panels.

## Details

For visualizing enrichment across multiple comparisons, use
[`multiplot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/multiplot_PA.md)
instead.

`theme_params` accepts any of the following named elements:

- `side_label_size`:

  Size for "Pathways" and "MSigDB" labels. Default: `35`.

- `geneset_text_size`:

  Text size for gene set labels. Default: `5`.

- `collection_text_size`:

  Text size for collection labels. Default: `5`.

- `panel_widths`:

  Patchwork relative widths for the 6 panels. Default:
  `c(4, 25, 15, 3, 10, 3)`.

- `col_size`:

  Border linewidth for `geom_col`. Default: `1`.

- `axis_title_size`:

  Font size for axis titles. Default: `45`.

- `axis_text_size_x`:

  Font size for x-axis labels. Default: `30`.

- `tick_size`:

  Linewidth for axis ticks. Default: `1.5`.

- `tick_length`:

  Length of axis ticks in cm. Default: `0.3`.

- `panel_spacing_single`:

  Spacing between facets. Default: `4`.

## See also

[`multiplot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/multiplot_PA.md)
for multi-comparison faceted barplots;
[`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md)
to generate the input data frame;
[camera_results](https://danielgarbozo.github.io/OmicsKit/reference/camera_results.md)
for a minimal example dataset.

## Examples

``` r
if (FALSE) { # \dontrun{
gsea_results <- merge_PA("path/to/gsea_results/")

# Filter to one comparison
single <- gsea_results[gsea_results$COMPARISON == "TumorVsNormal", ]

splot_PA(
  data           = single,
  geneset_col    = "NAME",
  collection_col = "COLLECTION",
  nes_col        = "NES",
  fdr_col        = "FDR"
)

# Cap color scale at -log10(FDR) = 5 so subtle differences are visible
# (gene sets with FDR <= 0.00001 all get the same max red color)
splot_PA(
  data        = single,
  geneset_col = "NAME", collection_col = "COLLECTION",
  nes_col     = "NES",  fdr_col        = "FDR",
  fill_limits = c(0, 5)
)
} # }
```
