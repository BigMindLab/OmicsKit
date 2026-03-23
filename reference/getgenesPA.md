# Extract gene members from pathway analysis results

For each gene set in a pathway analysis results table, retrieves leading
edge genes, a user-defined top fraction of genes, all genes in the gene
set, or any combination. All gene lists are ordered by their rank in the
provided ranked gene list.

## Usage

``` r
getgenesPA(pa_data, geneset_list, ranked_genes, genes = c("all", "le"))
```

## Arguments

- pa_data:

  A data frame of pathway analysis results. Must always contain:

  - `NAME`: gene set name.

  Additionally required depending on `genes`:

  - `SIZE`: number of genes in the gene set. Required for `"le"` and
    `"top"`.

  - `tags`: numeric fraction (0-1) of genes contributing to the
    enrichment score (GSEA leading edge). Produced automatically by
    [`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md).
    Required for `genes = "le"`.

  - `top`: numeric fraction (0-1) defining the proportion of top-ranked
    genes to extract. Set manually by the user (e.g.,
    `pa_data$top <- 0.25` for the top 25%). Required for
    `genes = "top"`.

  Typically the output of
  [`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md).

- geneset_list:

  A named list of gene sets, where each element is a character vector of
  gene symbols. Typically the output of
  [`list_gmts()`](https://danielgarbozo.github.io/OmicsKit/reference/list_gmts.md),
  or use the built-in
  [geneset_list](https://danielgarbozo.github.io/OmicsKit/reference/geneset_list.md)
  for quick testing.

- ranked_genes:

  A character vector of gene symbols ordered by their ranking metric
  (e.g., DESeq2 `stat`, log2FC, or signal-to-noise ratio), from most
  positive to most negative. Non-significant genes fall in the middle of
  the list. Used to order genes within each extracted set.

- genes:

  Character vector specifying which extraction mode(s) to use. Any
  combination of `"all"`, `"le"`, and `"top"`. Default:
  `c("all", "le")`.

## Value

Depends on `genes`:

- Single mode (e.g., `genes = "le"`): a named list where each element is
  a character vector of gene symbols. The list has an attribute
  `genes_type` used by
  [`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md)
  to name the output column.

- Multiple modes (e.g., `genes = c("all", "le", "top")`): a named list
  with one element per requested mode:

  - `$all`: named list of all gene set members.

  - `$le`: named list of leading edge genes (GSEA only).

  - `$top`: named list of top-ranked genes.

## Details

**Three extraction modes:**

- `"le"`: **GSEA output only.** Leading edge genes: the subset of genes
  that drives the enrichment signal. Size is computed as
  `round(SIZE * tags)`, where `tags` is the fraction of gene hits before
  (positive ES) or after (negative ES) the peak in the running
  enrichment score. Definition from the GSEA User Guide: *"The
  percentage of gene hits before (for positive ES) or after (for
  negative ES) the peak in the running enrichment score. This gives an
  indication of the percentage of genes contributing to the enrichment
  score."* (<https://docs.gsea-msigdb.org/#GSEA/GSEA_User_Guide/>).
  Requires columns `SIZE` and `tags` in `pa_data`, produced
  automatically by
  [`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md).

- `"top"`: **Any enrichment result (GSEA, CAMERA, PADOG, etc.).** A
  user-defined top fraction of genes ordered by rank. Size is computed
  as `round(SIZE * top)`, where `top` is a numeric value between 0 and 1
  provided in a `top` column of `pa_data`. This does **not** represent
  true leading edge genes: it is a flexible, rank-based selection
  suitable for exploratory visualization with any pathway analysis
  method. Requires columns `SIZE` and `top` in `pa_data`.

- `"all"`: All genes in the gene set, ordered by rank.

## See also

[`addgenesPA()`](https://danielgarbozo.github.io/OmicsKit/reference/addgenesPA.md)
to append gene columns to pa_data;
[`heatmap_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/heatmap_PA.md)
for heatmap visualization;
[`list_gmts()`](https://danielgarbozo.github.io/OmicsKit/reference/list_gmts.md)
to generate `geneset_list`;
[`merge_PA()`](https://danielgarbozo.github.io/OmicsKit/reference/merge_PA.md)
to generate `pa_data` with the required `tags` column.

## Examples

``` r
if (FALSE) { # \dontrun{
data(gsea_results)
data(geneset_list)
data(deseq2_results)

#or
gsl <- list_gmts("path/to/gmt_folder/")

ranked    <- deseq2_results$gene_id[order(deseq2_results$stat,
                                          decreasing = TRUE)]
pa_single <- gsea_results[gsea_results$COMPARISON == "TumorVsNormal", ]

# ── GSEA results: all three modes available
gene_lists <- getgenesPA(pa_single, geneset_list, ranked,
                         genes = c("all", "le", "top"))

# But first add the top column (e.g. top 30% of genes by rank)
pa_single$top <- 0.30
gene_lists <- getgenesPA(pa_single, geneset_list, ranked,
                         genes = c("all", "le", "top"))

gene_lists$le[["KEGG_APOPTOSIS"]]    # leading edge genes
gene_lists$top[["KEGG_APOPTOSIS"]]   # top 30% by rank
gene_lists$all[["KEGG_APOPTOSIS"]]   # all genes

pa_annot <- addgenesPA(pa_single, gene_lists)
head(pa_annot[, c("NAME", "all_genes", "le_genes", "top_genes")])

# ── CAMERA results: use "top" (no leading edge available) ───
data(camera_results)
camera_pa      <- camera_results
colnames(camera_pa)[colnames(camera_pa) == "GeneSet"] <- "NAME"
camera_pa$SIZE <- sapply(camera_pa$NAME,
                         function(x) length(geneset_list[[x]]))
camera_pa$top  <- 0.25   # top 25% by rank

gene_lists_cam <- getgenesPA(camera_pa, geneset_list, ranked,
                             genes = c("all", "top"))
pa_annot_cam   <- addgenesPA(camera_pa, gene_lists_cam)
head(pa_annot_cam[, c("NAME", "all_genes", "top_genes")])
} # }
```
