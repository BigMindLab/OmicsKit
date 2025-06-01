
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OmicsKit

<!-- badges: start -->

[![R-CMD-check](https://github.com/BigMindLab/OmicsKit/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BigMindLab/OmicsKit/actions/workflows/R-CMD-check.yaml)
[![License: CC
BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of `OmicsKit` is to help in manipulating tables and generating
plots for multi-omics analyses including genomics, transcriptomics,
proteomics, methylomics and immunoinformatics.

-----

## Developed by

  - David R. Requena Anicama, Ph.D.
    
      - Author’s name: David Requena
      - [Google
        Scholar](https://scholar.google.com/citations?user=uI01iS4AAAAJ&hl=en)
      - [ORCID: 0000-0002-5968-1133](https://orcid.org/0000-0002-5968-1133)

  - Daniel F. Guevara Díaz, B.Sc.
    
      - Author’s name: Daniel F. Guevara-Díaz
      - [Google
        Scholar](https://scholar.google.com/citations?hl=en&user=tqT7vr8AAAAJ)
      - [ORCID: 0009-0001-2786-8729](https://orcid.org/0009-0001-2786-8729)

## License

CC BY-NC-SA 4.0

## Contact

<david.requena@nyulangone.org>

## Installation

Install the dependencies with the following:

``` r
# Install CRAN packages
install.packages("dplyr")
install.packages("matrixStats")
install.packages("ggplot2")
install.packages("tsne")
install.packages("magrittr")
install.packages("rlang")
install.packages("stats")
install.packages("tibble")
install.packages("umap")
install.packages("lifecycle")
install.packages("ggrepel")
install.packages("openxlsx")
install.packages("remotes")

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
BiocManager::install("SummarizedExperiment")
```

You can install the development version of `OmicsKit` from
[GitHub](https://github.com/) with:

``` r
# Install from GitHub
remotes::install_github("BigMindLab/OmicsKit")

# Call library for usage
library("OmicsKit")
```
For a more detailed workflow on Differential Expression Analysis with
the application of the `OmicsKit` suit please check the custom
[BigMind](https://github.com/BigMindLab)’s pipeline for
[DESeq2](https://github.com/BigMindLab/DESeq2).
