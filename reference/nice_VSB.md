# Function to make Violin-Scatter-Box plots from data frames.

This function will make a Boxplot, using a DEseq object. It will show
the data points on top with a small deviation (jitter) for a better
visualization.

## Usage

``` r
nice_VSB(
  object = NULL,
  annotations,
  variables = c(fill = "VarFill", shape = "VarShape"),
  genename = NULL,
  symbol = NULL,
  labels = c("N", "P", "R", "M"),
  categories = c("normal", "primary", "recurrence", "metastasis"),
  colors = NULL,
  shapes = NULL,
  markersize = NULL,
  alpha = 0.8,
  jitter = 0.2,
  title_size = c(axis = 20, fig = 24),
  label_size = c(x = 20, y = 16),
  legend_size = c(title = 14, elements = 12)
)
```

## Arguments

- object:

  A data frame object with normalized counts genes(in rows) across
  samples(in columns).

- annotations:

  Data frame with annotations.

- variables:

  To indicate the variables to be used as Shape and Fill of the markers.

- genename:

  The gene name to be used for the plot.

- symbol:

  The gene symbol to display in the plot title. To obtain gene symbols
  from Ensembl IDs, use
  [`get_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/get_annotations.md).

- labels:

  A vector containing the x-labels of the box-plot. Default: c("N", "P",
  "R", "M").

- categories:

  A vector containing the labels for the legend. Default: c("normal",
  "primary", "recurrence", "metastasis").

- colors:

  Vector of colors to be used for the categories of the variable
  assigned as Marker Fill.

- shapes:

  Vector of shapes to be used for the categories of the variable
  assigned as Marker Shape.

- markersize:

  Size of the marker.

- alpha:

  Transparency of the marker, which goes from 0 (transparent) to 1 (no
  transparent). Default: 0.8.

- jitter:

  Random deviation added to the dots. Default: 0.2.

- title_size:

  Font of the title and axis names. Default: c(axis = 20, fig = 24).

- label_size:

  Font of the labels (x-axis) and numbers (y-axis). Default: c(x = 20, y
  = 16).

- legend_size:

  Font of the title and elements of the legend. Default: c(title = 14,
  elements = 12).

## Value

A ggplot2 object.

## See also

[`nice_Volcano()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_Volcano.md)
for genome-wide visualization;
[`detect_filter()`](https://danielgarbozo.github.io/OmicsKit/reference/detect_filter.md)
to identify reliably expressed genes;
[`get_stars()`](https://danielgarbozo.github.io/OmicsKit/reference/get_stars.md)
to add significance annotations;
[norm_counts](https://danielgarbozo.github.io/OmicsKit/reference/norm_counts.md)
for an example input matrix.

## Examples

``` r
data(norm_counts)
data(sampledata)

nice_VSB(
  object      = norm_counts,
  annotations = sampledata,
  variables   = c(fill = "sample_type"),
  genename    = rownames(norm_counts)[1],
  categories  = c("normal", "tumor"),
  labels      = c("Normal", "Tumor"),
  colors      = c("steelblue", "firebrick"),
  shapes      = 21,
  markersize  = 3
)
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the OmicsKit package.
#>   Please report the issue to the authors.

```
