# Plot leading edge heatmaps from pathway analysis results

Generates heatmaps of gene expression for each gene set in
`pa_data_annot`, using the `all_genes`, `le_genes` (GSEA output only),
and/or `top_genes` columns produced by
[`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md).
Genes within each heatmap are ordered by their position in
`ranked_genes`.

## Usage

``` r
heatmap_PA(
  expression_data,
  metadata,
  pa_data_annot,
  ranked_genes,
  plot_genes = c("all_genes", "le_genes"),
  sample_col = "Sample",
  group_col = "group",
  out_dir = "heatmaps_PA",
  pdf = TRUE,
  jpg = TRUE
)
```

## Arguments

- expression_data:

  A numeric matrix or data frame of expression values with gene symbols
  or Ensembl IDs as row names and sample IDs as column names.
  Recommended input: VST-transformed counts from
  [vst_counts](https://danielgarbozo.github.io/OmicsKit/reference/vst_counts.md)
  or normalized coutns
  [norm_counts](https://danielgarbozo.github.io/OmicsKit/reference/norm_counts.md).

- metadata:

  A data frame of sample annotations. Must contain a column matching
  `sample_col` (sample IDs) and a column matching `group_col` (condition
  labels, e.g., `"Control"`, `"Treatment"`).

- pa_data_annot:

  A data frame of pathway analysis results enriched with gene columns.
  Must contain the column `NAME` and at least one of `all_genes`,
  `le_genes`, or `top_genes` (comma-separated gene symbols per gene
  set). Typically the output of
  [`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md).

- ranked_genes:

  A character vector of gene symbols ordered by their ranking metric
  (e.g., stat, log2FC or signal-to-noise ratio), used to sort genes
  within each heatmap row.

- plot_genes:

  Character vector specifying which gene columns to plot. One or both of
  `"all_genes"` and `"le_genes"`, and `"top_genes"`. Each selection
  produces its own set of output files in a dedicated subfolder.
  Default: `c("all_genes", "le_genes")`.

- sample_col:

  Name of the sample ID column in `metadata`. Default: `"Sample"`.

- group_col:

  Name of the condition/group column in `metadata` (e.g., `"Control"` vs
  `"Treatment"`). Used for heatmap column annotations. Default:
  `"group"`.

- out_dir:

  Character. Path to the output directory. Subdirectories are created
  automatically based on `pdf`, `jpg`, and `plot_genes`:

  - `<out_dir>/pdf/all_genes/`

  - `<out_dir>/pdf/le_genes/`

  - `<out_dir>/pdf/top_genes/`

  - `<out_dir>/jpg/all_genes/`

  - `<out_dir>/jpg/le_genes/`

  - `<out_dir>/jpg/top_genes/` Default: `"heatmaps_PA"`.

- pdf:

  Logical. If `TRUE`, saves PDF heatmaps. Default: `TRUE`.

- jpg:

  Logical. If `TRUE`, saves JPG heatmaps. Default: `TRUE`.

## Value

Invisibly returns `TRUE` upon completion. Saves heatmap files to the
corresponding subdirectories under `out_dir`.

## Details

The recommended workflow before calling this function is:

    gsl          <- list_gmts("path/to/gmt/")
    pa_data      <- merge_PA("path/to/pa_data/")
    ranked       <- deseq2_results$gene_id[order(deseq2_results$stat,
                                                 decreasing = TRUE)]
    gene_lists   <- getgenesPA(pa_data, gsl, ranked, genes = c("all", "le"))
    pa_annot     <- addgenesPA(pa_data, gene_lists)

    heatmap_PA(
      expression_data = vst_counts,
      metadata        = sampledata,
      pa_data_annot   = pa_annot,
      ranked_genes    = ranked,
      plot_genes      = c("all_genes", "le_genes")
    )

## See also

[`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md)
for gene extraction;
[`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md)
to generate `pa_data_annot`;
[`list_gmts()`](https://danielgarbozo.github.io/OmicsKit/reference/list_gmts.md)
to generate the geneset list;
[`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md)
to generate `pa_data`;
[vst_counts](https://danielgarbozo.github.io/OmicsKit/reference/vst_counts.md)
for an example expression matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
data(vst_counts)
data(sampledata)
data(deseq2_results)
data(gsea_results)
data(geneset_list)

ranked    <- deseq2_results$gene_id[order(deseq2_results$stat,
                                          decreasing = TRUE)]

# ── Example 1: GSEA results (all_genes + le_genes) ────
pa_single  <- gsea_results[gsea_results$COMPARISON == "TumorVsNormal", ]
gene_lists <- getgenesPA(pa_single, geneset_list, ranked,
                         genes = c("all", "le"))
pa_annot   <- addgenesPA(pa_single, gene_lists)

heatmap_PA(
  expression_data = vst_counts,
  metadata        = sampledata,
  pa_data_annot   = pa_annot,
  ranked_genes    = ranked,
  plot_genes      = c("all_genes", "le_genes"),
  sample_col      = "patient_id",
  group_col       = "sample_type",
  out_dir         = "heatmaps_gsea",
  pdf             = TRUE,
  jpg             = TRUE
)
# Creates:
#   heatmaps_gsea/pdf/all_genes/<geneset>_heatmap.pdf
#   heatmaps_gsea/pdf/le_genes/<geneset>_heatmap.pdf
#   heatmaps_gsea/jpg/all_genes/<geneset>_heatmap.jpg
#   heatmaps_gsea/jpg/le_genes/<geneset>_heatmap.jpg

# ── Example 2: CAMERA results (all_genes + top_genes)
# camera_results does not contain leading edge information.
# Use genes = "top" with a manually set top fraction instead.
# Note: top_genes are rank-based and do NOT represent true leading edge genes.
data(camera_results)
camera_pa      <- camera_results
colnames(camera_pa)[colnames(camera_pa) == "GeneSet"] <- "NAME"
camera_pa$SIZE <- sapply(camera_pa$NAME,
                         function(x) length(geneset_list[[x]]))
camera_pa$top  <- 0.25

gene_lists_cam <- getgenesPA(camera_pa, geneset_list, ranked,
                             genes = c("all", "top"))
pa_annot_cam   <- addgenesPA(camera_pa, gene_lists_cam)

heatmap_PA(
  expression_data = vst_counts,
  metadata        = sampledata,
  pa_data_annot   = pa_annot_cam,
  ranked_genes    = ranked,
  plot_genes      = c("all_genes", "top_genes"),
  sample_col      = "patient_id",
  group_col       = "sample_type",
  out_dir         = "heatmaps_camera"
)
} # }
```
