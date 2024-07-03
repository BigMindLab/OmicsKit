# OmicsKit
A bioinformatics library for the analysis of omics data.

---

### Developed by:

David R. Requena Anicama, Ph.D.

- Author's name: David Requena / [Google Scholar](https://scholar.google.com/citations?user=uI01iS4AAAAJ&hl=en) / [ORCID: 0000-0002-5968-1133](https://orcid.org/0000-0002-5968-1133)

Daniel F. Guevara Díaz, B.Sc.(s)

- Author's name: Daniel F. Guevara-Díaz / [Google Scholar](https://scholar.google.com/citations?hl=en&user=tqT7vr8AAAAJ) / [ORCID: 0009-0001-2786-8729](https://orcid.org/0009-0001-2786-8729)

---

### Description

This library contains functions that help in manipulating tables and generating plots for transcriptomics and gene-set enrichment analysis.

It contains:
- Get information from Ensembl and BioMart to annotate gene counts tables. This information includes transcript and gene names, coordinates, and identifiers from different annotation databases.
- Functions to make visually appealing plots: PCA, tSNE, UMAP, Volcano, Balloon, and BSV (Box-Scatter-Violin) plots.
- Calculate and extract normalized counts, such as TPMs, RPKMs, FPKMs, and the normalized counts from DESeq2.
- Filter and export differential expression analysis results in MS Excel or CSV formal. Filtering can be done by expression change (log2FC), significance (FDR), or detectability (Requena et al. 2024. Manuscript under review).
- Organize the results from 3 pairwise differential expression analyses (i.e. B vs A, C vs A, C vs B) into 10 expression cases (BigMind. 2024. Manuscript in preparation.)
- Organize the results from 3 pairwise gene set enrichment analyses (i.e. B vs A, C vs A, C vs B) into 10 expression cases (BigMind. 2024. Manuscript in preparation.)

### License

CC BY-NC-SA 4.0

### Contact

david[dot]requena[at]nyulangone[dot]org
