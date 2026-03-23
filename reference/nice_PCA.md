# Function to make nice PCA plots

This was inspired on the plotPCA function from DESeq2, made by Wolfgang
Huber But including some improvements made by David Requena. Now it
allows:

- To choose which PCs to plot.

- To use one or two features to represent as the fill or shape of the
  markers.

- To provide the colors, shapes and fonts.

## Usage

``` r
nice_PCA(
  object,
  annotations = NULL,
  PCs = c(1, 2),
  ntop = NULL,
  variables = c(fill = "VarFill", shape = "VarShape"),
  legend_names = c(fill = "Sample Type", shape = "Library"),
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
  scale = FALSE,
  n_clusters = 3,
  transform = FALSE,
  outPCs = 50,
  returnData = FALSE
)
```

## Arguments

- object:

  A matrix of counts with genes as rows and sample ids as columns.

- annotations:

  A data frame of annotation, including sample ids and variables to
  plot. Default: NULL.

- PCs:

  A vector indicating the two Principal Components to plot. Default:
  c(1,2).

- ntop:

  Number of top genes to use for principal components, selected by
  highest row variance. Default: NULL.

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
  assigned as Marker Fill.

- shapes:

  Vector of shapes to be used for the categories of the variable
  assigned as Marker Shape.

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

- scale:

  Logical. Indicates whether to scale the data to have unit variances or
  not. Default: FALSE.

- n_clusters:

  Number of cluster categories. Default: 3.

- transform:

  Logical. Indicates whether to log2 transform the input `object` or
  not. Default: FALSE.

- outPCs:

  Number of Principal Components to keep if `returnData` is TRUE.
  Default: 50.

- returnData:

  Indicates if the function should return the data (TRUE) or the plot
  (FALSE). Default: FALSE.

## Value

A ggplot2 object if `returnData = FALSE` (default). If
`returnData = TRUE`, a numeric matrix of PCA coordinates with dimensions
samples × `outPCs`, with a `percentVar` attribute containing the
proportion of variance explained per component.

## See also

[`nice_UMAP()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_UMAP.md),
[`nice_tSNE()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_tSNE.md)
for other alternatives;
[vst_counts](https://danielgarbozo.github.io/OmicsKit/reference/vst_counts.md)
for the recommended input matrix.

## Examples

``` r
data(vst_counts)
data(sampledata)

# nice_PCA joins by a column named "id" in annotations
sampledata_pca <- sampledata
colnames(sampledata_pca)[colnames(sampledata_pca) == "patient_id"] <- "id"

nice_PCA(
  object       = vst_counts,
  annotations  = sampledata_pca,
  variables    = c(fill = "sample_type"),
  legend_names = c(fill = "Sample Type"),
  colors       = c("steelblue", "firebrick"),
  shapes       = c(21, 21),
  title        = "TCGA-LUAD PCA"
)


# Return PCA coordinates instead of plot
pca_data <- nice_PCA(
  object       = vst_counts,
  annotations  = sampledata_pca,
  variables    = c(fill = "sample_type"),
  legend_names = c(fill = "Sample Type"),
  colors       = c("steelblue", "firebrick"),
  shapes       = c(21, 21),
  returnData   = TRUE
)
head(pca_data)
#>                        PC1       PC2        PC3       PC4        PC5        PC6
#> TCGA.38.4627.11A -78.79063 59.286255 -12.535203  45.69875  -8.521804 31.5738405
#> TCGA.44.2661.11A -72.45790  5.943699 -20.779773   6.92918  16.160545  0.4116424
#> TCGA.44.2662.11A -96.38040 10.479110 -17.724971  44.79115  30.191426  5.4811386
#> TCGA.44.2665.11A -75.59222  5.412153  -4.878161  18.38946  26.132611  1.4414801
#> TCGA.44.3396.11A -83.09770  1.279103   4.235361 -13.36483   2.182167 -2.3089167
#> TCGA.49.4490.11A -79.21172 48.097839   4.554704  38.25393 -29.732007 27.7126528
#>                          PC7          PC8         PC9       PC10      PC11
#> TCGA.38.4627.11A -24.8811054  25.17386382 -18.8116041   5.287110  8.995741
#> TCGA.44.2661.11A   2.3514439  -0.09795269  27.7216192  -7.955043 -6.288688
#> TCGA.44.2662.11A   0.6071834 -24.82756200  23.2520922 -24.565973 -5.875162
#> TCGA.44.2665.11A  -3.4136739  -4.92306119  22.4201811 -15.787971 -4.662476
#> TCGA.44.3396.11A  23.1974629  -1.36924768   0.4289105   8.639563  2.303862
#> TCGA.49.4490.11A -18.4567315  44.37886818  -9.1173455  32.866384  2.448545
#>                        PC12       PC13      PC14       PC15       PC16
#> TCGA.38.4627.11A -27.247461 -23.282455 13.627227  10.033671  2.1296466
#> TCGA.44.2661.11A   1.737105  -3.823004  2.451978   3.167976 -0.4638886
#> TCGA.44.2662.11A  12.768091  -6.553733  3.264268   9.912506  4.0648559
#> TCGA.44.2665.11A   2.590769   4.443403  6.650630   6.212390 -1.7501553
#> TCGA.44.3396.11A -15.512437   0.298238  1.815369   1.570493 -5.7752545
#> TCGA.49.4490.11A -31.561690 -25.225557 12.121584 -16.891459  4.6797697
#>                        PC17       PC18       PC19      PC20       PC21
#> TCGA.38.4627.11A   6.616933  -2.141496  -6.080581 -5.543456 -10.666033
#> TCGA.44.2661.11A   6.726650  -4.154546 -13.310570 21.211995 -25.943320
#> TCGA.44.2662.11A -22.862328  -3.669719  -5.370658 -1.449050  27.634535
#> TCGA.44.2665.11A  -2.470355 -19.226250 -29.554969  6.935876 -18.848270
#> TCGA.44.3396.11A -12.633249 -12.920533  -7.301930 -2.138942   6.148496
#> TCGA.49.4490.11A  23.403437  11.491121   7.785900 -4.927255   9.185427
#>                        PC22       PC23       PC24       PC25         PC26
#> TCGA.38.4627.11A -12.757711  31.969389 -29.012472 -16.640074 -11.73512008
#> TCGA.44.2661.11A -17.314486  -7.799813 -10.727717  -8.499543  -5.61731976
#> TCGA.44.2662.11A  11.903944  19.740749 -12.970173  29.226850  -0.05879115
#> TCGA.44.2665.11A  -5.227457 -13.646929  25.116920 -11.061594  -9.27940256
#> TCGA.44.3396.11A  -1.700231   9.263947  26.591122 -25.708918  -3.29802557
#> TCGA.49.4490.11A  11.004490 -37.102681   8.257585  17.891726  -0.47593621
#>                        PC27       PC28       PC29        PC30       PC31
#> TCGA.38.4627.11A -13.085979   5.981191   7.908326   6.5438103 -3.0156000
#> TCGA.44.2661.11A  13.901876   1.180475 -28.600207 -28.0720520  6.9040333
#> TCGA.44.2662.11A  25.640151  -3.361286  -6.938338  11.9338064  0.4268999
#> TCGA.44.2665.11A  16.465812   2.770048  34.558252   2.9641156 -4.0022625
#> TCGA.44.3396.11A  -7.296370 -33.642971 -19.531914  16.8185457  2.3830728
#> TCGA.49.4490.11A   9.459282  -8.336114  -2.189187   0.6766476  2.2589409
#>                           PC32
#> TCGA.38.4627.11A -3.614341e-13
#> TCGA.44.2661.11A -3.040833e-14
#> TCGA.44.2662.11A -2.852745e-14
#> TCGA.44.2665.11A -1.366383e-13
#> TCGA.44.3396.11A -1.641928e-13
#> TCGA.49.4490.11A -1.264051e-13
```
