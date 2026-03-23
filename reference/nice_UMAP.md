# Function to make UMAP plots.

Function to make UMAP plots.

## Usage

``` r
nice_UMAP(
  object,
  annotations = NULL,
  neighbors = 5,
  components = 2,
  epochs = 10000,
  seed = 0,
  variables = c(fill = "VarFill", shape = "VarShape"),
  legend_names = c(fill = "Label Fill", shape = "Label Shape"),
  size = 5,
  alpha = 1,
  colors = NULL,
  shapes = NULL,
  title = NULL,
  legend_title = 16,
  legend_elements = 14,
  legend_pos = NULL,
  labels = NULL,
  name_tags = NULL,
  cluster_data = FALSE,
  min_points = 7,
  transform = FALSE,
  returnData = FALSE
)
```

## Arguments

- object:

  A matrix of counts with genes as rows and sample ids as columns.

- annotations:

  A data frame of annotation, including sample ids and variables to
  plot. Default: NULL.

- neighbors:

  Number of nearest neighbors to consider. Default: 5.

- components:

  Number of components to be consider for the dimensional reduction.
  Default: 3.

- epochs:

  Number of iterations. Default: 10000.

- seed:

  Set a random seed state. Default: 0.

- variables:

  To indicate the variables to be used as Shape and Fill of the markers.

- legend_names:

  The names to be used for the legend of the Shape and Fill.

- size:

  Size of the marker. Default: 5.

- alpha:

  Transparency of the marker, which goes from 0 (transparent) to 1 (no
  transparent). Default: 1.

- colors:

  Vector of colors to be used for the categories of the variable
  assigned as VarFill.

- shapes:

  Vector of shapes to be used for the categories of the variable
  assigned as VarShape.

- title:

  Plot title. Default: NULL.

- legend_title:

  Font of the legend title. Default: 16.

- legend_elements:

  Font of the elements of the legend Default: 14.

- legend_pos:

  Position of the legend inside the plot. Example: c(0.80, 0.80).
  Default: NULL.

- labels:

  A vector containing the variable to be used as labels (name inside the
  marker), and the label size. Example: c(var = "patient", size = 2).
  Default: NULL (no labels).

- name_tags:

  A vector containing the variable to be used as name tags (name outside
  the marker), tag size, minimum distance in order to add an arrow
  connecting the tag and the marker, and minimum distance from the tag
  and the center of the marker. Example: c(var = "label", size = 3,
  minlen = 2, box = 0.5). Default: NULL (no name tags).

- cluster_data:

  Indicates if the function generates the clusters (TRUE) or not
  (FALSE). This new cluster variable can be used as fill or shape.
  Default: FALSE.

- min_points:

  Minimum number of neighbors to form a cluster. Default: 7.

- transform:

  Logical. Indicates whether to log2 transform the input `object` or
  not. Default: FALSE.

- returnData:

  Indicates if the function should return the data (TRUE) or the plot
  (FALSE). Default: FALSE.

## Value

A ggplot2 object if `returnData = FALSE` (default). If
`returnData = TRUE`, a data frame with UMAP coordinates and sample
annotations.

## References

McInnes, L., Healy, J., & Melville, J. (2018). Umap: Uniform Manifold
Approximation and Projection for Dimension Reduction. *arXiv preprint
arXiv:1802.03426*. <https://arxiv.org/abs/1802.03426>

## See also

[`nice_PCA()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_PCA.md),
[`nice_tSNE()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_tSNE.md)
for alternative dimensionality reduction methods;
[vst_counts](https://danielgarbozo.github.io/OmicsKit/reference/vst_counts.md)
for the recommended input matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
data(vst_counts)
data(sampledata)

sampledata_u <- sampledata
colnames(sampledata_u)[colnames(sampledata_u) == "patient_id"] <- "id"

nice_UMAP(
  object       = vst_counts,
  annotations  = sampledata_u,
  variables    = c(fill = "sample_type"),
  legend_names = c(fill = "Sample Type"),
  colors       = c("steelblue", "firebrick"),
  shapes       = c(21, 21),
  title        = "TCGA-LUAD UMAP",
  neighbors    = 5,
  epochs       = 1000,
  seed         = 1905
)
} # }
```
