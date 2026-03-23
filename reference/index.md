# Package index

## Differential Expression Analysis (DEA)

Functions for data quality control, normalization, dimensionality
reduction, annotation, and differential expression visualization.

- [`power_analysis()`](https://danielgarbozo.github.io/OmicsKit/reference/power_analysis.md)
  : Power analysis for RNA-seq differential expression with optional
  plotting
- [`tpm()`](https://danielgarbozo.github.io/OmicsKit/reference/tpm.md) :
  Function to calculate the TPMs from a table of raw gene counts.
- [`nice_PCA()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_PCA.md)
  : Function to make nice PCA plots
- [`nice_UMAP()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_UMAP.md)
  : Function to make UMAP plots.
- [`nice_tSNE()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_tSNE.md)
  : Function to make tSNE plots.
- [`get_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/get_annotations.md)
  : Get annotations from Ensembl.
- [`add_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/add_annotations.md)
  : A function to add annotations to a table of gene counts.
- [`save_results()`](https://danielgarbozo.github.io/OmicsKit/reference/save_results.md)
  : Save the results of a Differential Expression analysis
- [`split_cases()`](https://danielgarbozo.github.io/OmicsKit/reference/split_cases.md)
  : Obtain 10 exclusive cases from 3 comparisons.
- [`nice_Volcano()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_Volcano.md)
  : Function to draw Volcano plots.
- [`detect_filter()`](https://danielgarbozo.github.io/OmicsKit/reference/detect_filter.md)
  : Find detectable genes across comparisons.
- [`get_stars()`](https://danielgarbozo.github.io/OmicsKit/reference/get_stars.md)
  : Obtain annotations to display significance.
- [`nice_VSB()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_VSB.md)
  : Function to make Violin-Scatter-Box plots from data frames.
- [`nice_VSB_DEseq2()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_VSB_DEseq2.md)
  : Function to make Box-Scatter-Violin plots from DEseq2 output
  directly.

## Genomics

Survival analysis and genomics visualization utilities.

- [`nice_KM()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_KM.md)
  : Function to generate Kaplan Meier survival plots for a given binary
  gene variable

## Pathway Analysis (PA)

Tools for loading, merging, and visualizing GSEA / pathway analysis
results, including single- and multi-comparison plots and heatmaps.

- [`list_gmts()`](https://danielgarbozo.github.io/OmicsKit/reference/list_gmts.md)
  : Read GMT files from a directory into a named gene set list
- [`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md)
  : Merge GSEA result files into a single data frame
- [`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md)
  : Extract gene members from pathway analysis results
- [`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md)
  : Add gene columns to pathway analysis results
- [`multiplot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/multiplot_PA.md)
  : Pathway analysis visualization across multiple comparisons
- [`splot_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/splot_PA.md)
  : Pathway analysis visualization for a single comparison
- [`heatmap_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_PA.md)
  : Plot leading edge heatmaps from pathway analysis results
- [`heatmap_path_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_path_PA.md)
  : Plot leading edge heatmaps from GSEA analysis results using file
  paths

## PA Clustering

Pathway clustering via Jaccard similarity, hierarchical clustering, and
network-based community detection.

- [`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md)
  : Compute Jaccard similarity and distance matrices for gene sets
- [`do_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/do_clust.md)
  : Hierarchical clustering of gene sets with silhouette-based
  optimization
- [`get_superterm()`](https://danielgarbozo.github.io/OmicsKit/reference/get_superterm.md)
  : Generate representative super-term labels for gene set communities
- [`get_network_communities()`](https://danielgarbozo.github.io/OmicsKit/reference/get_network_communities.md)
  : Detect gene set communities and generate super-term labels
- [`network_clust()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust.md)
  : Gene set network clustering with igraph base R graphics
- [`network_clust_gg()`](https://danielgarbozo.github.io/OmicsKit/reference/network_clust_gg.md)
  : Gene set network clustering with ggraph graphics

## Example datasets

Built-in datasets for reproducible examples and vignettes.

- [`camera_results`](https://danielgarbozo.github.io/OmicsKit/reference/camera_results.md)
  : Example CAMERA enrichment results for pathway analysis clustering
- [`deseq2_results`](https://danielgarbozo.github.io/OmicsKit/reference/deseq2_results.md)
  : DESeq2 differential expression results for TCGA-LUAD
- [`geneset_list`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_list.md)
  : Example gene set list for pathway analysis clustering
- [`gsea_results`](https://danielgarbozo.github.io/OmicsKit/reference/gsea_results.md)
  : Simulated GSEA pathway analysis results for TCGA-LUAD
- [`norm_counts`](https://danielgarbozo.github.io/OmicsKit/reference/norm_counts.md)
  : Normalized counts matrix for TCGA-LUAD
- [`raw_counts`](https://danielgarbozo.github.io/OmicsKit/reference/raw_counts.md)
  : GDC TCGA Lung Adenocarcinoma (LUAD) - Raw STAR counts
- [`sampledata`](https://danielgarbozo.github.io/OmicsKit/reference/sampledata.md)
  : GDC TCGA Lung Adenocarcinoma (LUAD) - Metadata
- [`vst_counts`](https://danielgarbozo.github.io/OmicsKit/reference/vst_counts.md)
  : Variance-stabilized counts matrix for TCGA-LUAD
