
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OmicsKit

<!-- badges: start -->
<!-- badges: end -->

The goal of OmicsKit is to help in manipulating tables and generating
plots for multi-omics analyses including genomics, transcriptomics,
proteomics, methylomics and immunoinformatics.

## Installation

You can install the development version of OmicsKit from
[GitHub](https://github.com/) with:

``` r
# Install the remotes package if needed
install.packages("remotes")

# Install from GitHub
remotes::install_github("BigMindLab/OmicsKit")
```

## Key features

- **Gene Annotation**: Retrieve information from Ensembl and BioMart to
  annotate gene counts tables, including transcript and gene names,
  genomic coordinates, and cross-references from various annotation
  databases.

``` r
# Example on ge

tx2gene <- get_annotations(rownames(txi$counts), version = "103", format = "xlsx")
```

- **Dimensionality Reduction**: Generate a range of visually appealing
  plots to visualize high-dimensional data. Includes unsupervised
  clustering methods as well.
  - PCA (Principal Component Analysis)

``` r
nice_PCA(transf.data,
         PCs = c(1, 2),
         ntop =  nrow(assay(transf.data)),
         variables = c("group", "sex"),
         legend_names = c("Group", "Sex"),
         size = 9, alpha = 1,
         colors = my_colors,
         shapes = 21:22,
         legend_title = 10,
         legend_elements = 8,
         legend_pos = NULL,
         labels = c(var = "id", size = 3))
```

    + tSNE (t-distributed Stochastic Neighbor Embedding)

    + UMAP (Uniform Manifold Approximation and Projection)

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(OmicsKit)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
