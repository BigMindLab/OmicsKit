# Function to make Box-Scatter-Violin plots from DEseq2 output directly.

This function will make a Boxplot, using a DEseq object. It will show
the data points on top with a small deviation (jitter) for a better
visualization.

## Usage

``` r
nice_VSB_DEseq2(
  object = NULL,
  variables = c(fill = "VarFill", shape = "VarShape"),
  genename = NULL,
  symbol = NULL,
  labels = c("N", "P", "R", "M"),
  categories = c("normal", "primary", "recurrence", "metastasis"),
  colors = NULL,
  shapes = NULL,
  markersize = NULL,
  alpha = 0.8,
  width = NULL,
  height = NULL,
  jitter = 0.2,
  dpi = 150,
  save = FALSE,
  title_size = c(axis = 20, fig = 24),
  label_size = c(x = 20, y = 16),
  legend_size = c(title = 14, elements = 12)
)
```

## Arguments

- object:

  A DEseq object already transformed with the variance stabilizing or
  rlog transformations.

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

- width:

  Width of the plot.

- height:

  Height of the plot.

- jitter:

  Random deviation added to the dots. Default: 0.2.

- dpi:

  DPI of the plot. Default: 150.

- save:

  To save the plot. Default: FALSE.

- title_size:

  Font of the title and axis names. Default: c(axis = 20, fig = 24).

- label_size:

  Font of the labels (x-axis) and numbers (y-axis). Default: c(x = 20, y
  = 16).

- legend_size:

  Font of the title and elements of the legend. Default: c(title = 14,
  elements = 12).

## Examples

``` r
if (FALSE) { # \dontrun{
# requires a DESeq2 object

data(sampledata)

nice_VSB_DEseq2(
  object      = vst,
  annotations = sampledata,
  variables   = c(fill = "sample_type"),
  genename    = rownames(norm_counts)[1],
  categories  = c("normal", "tumor"),
  labels      = c("Normal", "Tumor"),
  colors      = c("steelblue", "firebrick"),
  shapes      = 21,
  markersize  = 3
)
} # }
```
