# OmicsKit

`OmicsKit` is an R package that streamlines multi-omics analyses in
genomics, transcriptomics, proteomics, and methylomics by simplifying
data handling, generating publication-quality visualizations, and
offering flexible utilities to accelerate end-to-end analysis.

------------------------------------------------------------------------

## Developed by

- David R. Requena Anicama, PhD

  - Author’s name: David Requena
  - [Google
    Scholar](https://scholar.google.com/citations?user=uI01iS4AAAAJ&hl=en)
  - [ORCID: 0000-0002-5968-1133](https://orcid.org/0000-0002-5968-1133)

- Daniel F. Guevara Díaz, BSc

  - Author’s name: Daniel F. Guevara-Diaz
  - [Google
    Scholar](https://scholar.google.com/citations?hl=en&user=tqT7vr8AAAAJ)
  - [ORCID: 0009-0001-2786-8729](https://orcid.org/0009-0001-2786-8729)

- Daniel E. Garbozo Santillan, BSc

  - Author’s name: Daniel E. Garbozo
  - [Google Scholar](https://danielgarbozo.github.io/OmicsKit/)
  - [ORCID: 0009-0003-2495-6568](https://orcid.org/0009-0003-2495-6568)

- Angela D. C. Alarcon Guerrero, BSc(s)

  - Author’s name: Angela D. C. Alarcon Guerrero
  - [Google Scholar](https://danielgarbozo.github.io/OmicsKit/)
  - [ORCID: 0000-0003-0293-5603](https://orcid.org/0000-0003-0293-5603)

## License

CC BY-NC-SA 4.0

## Contact

<david.requena@nyulangone.org>

## Installation

Install the dependencies with the following:

``` r
# Install CRAN packages
install.packages("ggplot2")
install.packages("dplyr")
install.packages("stats")
install.packages("matrixStats")
install.packages("umap")
install.packages("tsne")
install.packages("magrittr")
install.packages("rlang")
install.packages("tibble")
install.packages("remotes")
```

You can install the development version of `OmicsKit` from
[GitHub](https://github.com/) with:

``` r
# Install from GitHub
remotes::install_github("BigMindLab/OmicsKit")

# Call library for usage
library("OmicsKit")
```

For a more detailed Differential Expression Analysis workflow using the
`OmicsKit` suite, see [BigMind](https://github.com/BigMindLab)’s GitHub
organization and its custom [DESeq2
pipeline](https://github.com/BigMindLab/DESeq2).
