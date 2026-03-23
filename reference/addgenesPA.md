# Add gene columns to pathway analysis results

Appends `all_genes`, `le_genes`, and/or `top_genes` columns to a pathway
analysis results data frame based on the output of
[`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md).
Gene symbols within each cell are comma-separated. Automatically detects
which column(s) to add based on the structure of the input.

## Usage

``` r
addgenesPA(pa_data, gene_lists)
```

## Arguments

- pa_data:

  A data frame of pathway analysis results containing a `NAME` column.
  Typically the output of
  [`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md).

- gene_lists:

  Output of
  [`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md).
  Can be:

  - A list with `$all`, `$le`, and/or `$top` slots: when multiple modes
    are requested (e.g.,
    `getgenesPA(..., genes = c("all", "le", "top"))`). Adds the
    corresponding columns.

  - A flat named list with attribute `genes_type`: when a single mode is
    requested. Adds the corresponding column (`all_genes`, `le_genes`,
    or `top_genes`).

## Value

The input `pa_data` data frame with one or more additional columns:

- `all_genes`: comma-separated string of all gene set members ordered by
  rank.

- `le_genes`: comma-separated string of leading edge genes (GSEA only),
  ordered by rank.

- `top_genes`: comma-separated string of top-ranked genes based on the
  user-defined `top` fraction.

Gene sets not found in `gene_lists` receive `NA`.

## See also

[`getgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/getgenesPA.md)
to generate `gene_lists`;
[`heatmap_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_PA.md)
for heatmap visualization;
[`save_results()`](https://danielgarbozo.github.io/OmicsKit/reference/save_results.md)
to export the annotated results.

## Examples

``` r
if (FALSE) { # \dontrun{
data(gsea_results)
data(geneset_list)
data(deseq2_results)

ranked    <- deseq2_results$gene_id[order(deseq2_results$stat,
                                          decreasing = TRUE)]
pa_single <- gsea_results[gsea_results$COMPARISON == "TumorVsNormal", ]
pa_single$top <- 0.30

# Add all three columns
gene_lists <- getgenesPA(pa_single, geneset_list, ranked,
                         genes = c("all", "le", "top"))
pa_annot   <- addgenesPA(pa_single, gene_lists)
head(pa_annot[, c("NAME", "all_genes", "le_genes", "top_genes")])

# Add only leading edge genes
le_only  <- getgenesPA(pa_single, geneset_list, ranked, genes = "le")
pa_annot <- addgenesPA(pa_single, le_only)
head(pa_annot[, c("NAME", "le_genes")])

# CAMERA: add only top and all (no leading edge)
data(camera_results)
camera_pa      <- camera_results
colnames(camera_pa)[colnames(camera_pa) == "GeneSet"] <- "NAME"
camera_pa$SIZE <- sapply(camera_pa$NAME,
                         function(x) length(geneset_list[[x]]))
camera_pa$top  <- 0.25
gene_lists_cam <- getgenesPA(camera_pa, geneset_list, ranked,
                             genes = c("all", "top"))
pa_annot_cam   <- addgenesPA(camera_pa, gene_lists_cam)
head(pa_annot_cam[, c("NAME", "all_genes", "top_genes")])
} # }
```
