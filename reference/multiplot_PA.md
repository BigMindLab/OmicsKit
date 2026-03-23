# Pathway analysis visualization across multiple comparisons

Generates a faceted barplot showing NES values across multiple
comparisons for a set of gene sets. Each facet represents one gene set
and bars represent the NES per comparison, colored by -log10(FDR). This
layout makes it easy to compare how enrichment of gene sets changes
across conditions (e.g., TumorVsNormal, MetastasisVsNormal).

## Usage

``` r
multiplot_PA(
  data,
  comparison_col = "COMPARISON",
  facet_col = "NAME",
  axis_y = "NES",
  fdr_col = "FDR",
  comparison_order = NULL,
  custom_labels = NULL,
  ncol_wrap = 2,
  free_y = TRUE,
  fill_limits = NULL,
  fill_palette = c("white", "red"),
  theme_params = list()
)
```

## Arguments

- data:

  A data frame of pathway analysis results containing two or more
  comparisons. Typically the output of
  [`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md).

- comparison_col:

  Name of the column identifying each comparison. Appears on the x-axis
  of each facet. Default: `"COMPARISON"`.

- facet_col:

  Name of the column used to define facets one facet per unique value.
  Can be the original gene set name column (e.g., `"NAME"`) or a
  manually curated column with cleaner or shorter labels (e.g.,
  `"clean_name"`). Default: `"NAME"`.

- axis_y:

  Name of the column to use for the y-axis. Default: `"NES"`.

- fdr_col:

  Name of the column containing FDR values. `-log10(FDR)` is computed
  internally and used as the fill color. Default: `"FDR"`.

- comparison_order:

  Character vector specifying the left-to-right order of comparisons on
  the x-axis of each facet. For example,
  `comparison_order = c("BvsA", "CvsA")` places `BvsA` on the left and
  `CvsA` on the right. If `NULL` (default), the order follows the factor
  levels of `comparison_col` as they appear in `data`.

- custom_labels:

  Named character vector of x-axis tick labels. Useful for shortening
  comparison names on the axis. For example,
  `custom_labels = c(TumorVsNormal = "Tumor", MetastasisVsNormal = "Mets")`.
  Default: `NULL`.

- ncol_wrap:

  Integer. Number of columns in `facet_wrap`. Default: `2`.

- free_y:

  Logical. If `TRUE`, each facet uses its own y-axis scale. Default:
  `TRUE`.

- fill_limits:

  Numeric vector of length 2 setting the color scale range for
  `-log10(FDR)`. Values outside this range are clamped to the nearest
  limit. For example, `fill_limits = c(0, 5)` maps all gene sets with
  `-log10(FDR) >= 5` (FDR \<= 0.00001) to maximum red, and any value
  below 0 to white. Useful when one gene set has extreme significance
  that makes the rest appear uniform. Default: `NULL` (auto).

- fill_palette:

  Character vector of two colors for the fill gradient (low to high
  -log10(FDR)). Default: `c("white", "red")`.

- theme_params:

  Named list to override default theme parameters. See Details.

## Value

A ggplot2 object.

## Details

All comparisons must be combined in a single data frame with a column
identifying each comparison as produced by
[`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md).

For visualizing a single comparison with full collection grouping, use
[`splot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/splot_PA.md)
instead.

`theme_params` accepts any of the following named elements:

- `bar_col`:

  Bar border color. Default: `"black"`.

- `bar_size`:

  Bar border linewidth. Default: `0.5`.

- `bar_width`:

  Bar width. Default: `0.6`.

- `hline_size`:

  Linewidth for horizontal line at y = 0. Default: `2`.

- `axis_title_size`:

  Font size for axis titles. Default: `45`.

- `axis_text_size_x`:

  Font size for x-axis labels. Default: `30`.

- `axis_text_size_y`:

  Font size for y-axis labels. Default: `50`.

- `tick_size`:

  Linewidth for axis ticks. Default: `1.5`.

- `tick_length`:

  Length of axis ticks in cm. Default: `0.3`.

- `strip_text_size`:

  Font size for facet strip labels. Default: `50`.

- `panel_spacing_multi`:

  Spacing between facets. Default: `0.6`.

## See also

[`splot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/splot_PA.md)
for single-comparison patchwork plots;
[`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md)
to generate the input data frame;
[camera_results](https://danielgarbozo.github.io/OmicsKit/reference/camera_results.md)
for a minimal example dataset.

## Examples

``` r
if (FALSE) { # \dontrun{
gsea_results <- merge_PA("path/to/gsea_results/")

# Basic multi-comparison plot
multiplot_PA(
  data           = gsea_results,
  comparison_col = "COMPARISON",
  facet_col      = "NAME",
  fdr_col        = "FDR",
  ncol_wrap      = 3
)

# Control left-to-right order of comparisons on the x-axis
multiplot_PA(
  data             = gsea_results,
  comparison_col   = "COMPARISON",
  facet_col        = "NAME",
  fdr_col          = "FDR",
  comparison_order = c("BvsA", "CvsA")   # BvsA on the left, CvsA on the right
)

# Use cleaner facet labels and shorten x-axis tick names
gsea_results$clean_name <- gsub("_", " ", gsea_results$NAME)

multiplot_PA(
  data             = gsea_results,
  comparison_col   = "COMPARISON",
  facet_col        = "clean_name",
  fdr_col          = "FDR",
  comparison_order = c("BvsA", "CvsA"),
  custom_labels    = c(BvsA = "Tumor", CvsA = "Metastasis")
)
} # }
```
