## generate_hallmark_heatmap.R
## Static Hallmark heatmap for pathway_analysis.Rmd
## Output:
## vignettes/figures/heatmaps/HALLMARK_ESTROGEN_RESPONSE_LATE_heatmap.jpg

suppressPackageStartupMessages({
  library(OmicsKit)
})

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop(
    "Package 'pheatmap' is required. ",
    "Install it with install.packages('pheatmap') and run this script again."
  )
}

## -----------------------------
## 1. Parameters
## -----------------------------

target_pathway    <- "HALLMARK_ESTROGEN_RESPONSE_LATE"
target_comparison <- "ERpositive_vs_ERnegative"
target_collection <- "Hallmark"

out_dir <- file.path("vignettes", "figures", "heatmaps")
target_file <- file.path(out_dir, paste0(target_pathway, "_heatmap.jpg"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## 2. Load data
## -----------------------------

data("brca_pa_merged")
data("brca_geneset_list")
data("brca_ranked_genes")
data("brca_rna_vst_or_logexpr_small")
data("brca_rna_metadata_tumor_normal")

## -----------------------------
## 3. Helper functions
## -----------------------------

pick_col <- function(df, candidates, label) {
  hits <- candidates[candidates %in% names(df)]

  if (length(hits) == 0) {
    stop(
      "Could not find column for: ", label, "\n",
      "Tried: ", paste(candidates, collapse = ", "), "\n\n",
      "Available columns are:\n",
      paste(names(df), collapse = ", ")
    )
  }

  message(label, " column selected: ", hits[1])
  hits[1]
}

clean_er_status <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x_low <- tolower(x)

  out <- rep(NA_character_, length(x))

  out[grepl("positive|pos|\\ber\\+\\b", x_low)] <- "ER Positive"
  out[grepl("negative|neg|\\ber-\\b", x_low)] <- "ER Negative"

  out
}

## -----------------------------
## 4. Check required columns
## -----------------------------

required_pa_cols <- c("NAME", "COLLECTION", "COMPARISON")
missing_pa_cols <- setdiff(required_pa_cols, names(brca_pa_merged))

if (length(missing_pa_cols) > 0) {
  stop(
    "brca_pa_merged is missing required columns: ",
    paste(missing_pa_cols, collapse = ", ")
  )
}

sample_col <- pick_col(
  brca_rna_metadata_tumor_normal,
  candidates = c(
    "sampleID",
    "sample_id",
    "sample",
    "Sample",
    "barcode",
    "submitter_id",
    "sample_submitter_id",
    "Tumor_Sample_Barcode"
  ),
  label = "Sample ID"
)

er_col <- pick_col(
  brca_rna_metadata_tumor_normal,
  candidates = c(
    "ER_group",
    "ER_status_raw",
    "ER_status",
    "er_status",
    "ER.Status",
    "er.status",
    "estrogen_receptor_status",
    "Estrogen_Receptor_Status",
    "breast_carcinoma_estrogen_receptor_status",
    "breast_carcinoma_estrogen_receptor_status_by_ihc",
    "er_status_by_ihc",
    "ER_Status_By_IHC"
  ),
  label = "ER status"
)

if (is.null(rownames(brca_rna_vst_or_logexpr_small))) {
  stop("Expression matrix must have gene symbols as row names.")
}

if (is.null(colnames(brca_rna_vst_or_logexpr_small))) {
  stop("Expression matrix must have sample IDs as column names.")
}

## -----------------------------
## 5. Select pathway
## -----------------------------

pa_single <- subset(
  brca_pa_merged,
  NAME == target_pathway &
    COLLECTION == target_collection &
    COMPARISON == target_comparison
)

if (nrow(pa_single) == 0) {
  stop(
    "Target pathway was not found in brca_pa_merged: ",
    target_pathway,
    " / ",
    target_comparison
  )
}

## -----------------------------
## 6. Extract leading-edge genes
## -----------------------------
## Use only leading-edge genes.
## This avoids requiring SIZE/top columns.

gene_lists <- getgenesPA(
  pa_data      = pa_single,
  geneset_list = brca_geneset_list,
  ranked_genes = brca_ranked_genes[[target_comparison]],
  genes        = "le"
)

pa_annotated <- addgenesPA(
  pa_data    = pa_single,
  gene_lists = gene_lists
)

if (!"le_genes" %in% names(pa_annotated)) {
  stop("addgenesPA() did not create the 'le_genes' column.")
}

leading_edge_genes <- unique(unlist(strsplit(pa_annotated$le_genes[1], "[,;/|[:space:]]+")))
leading_edge_genes <- trimws(leading_edge_genes)
leading_edge_genes <- leading_edge_genes[nzchar(leading_edge_genes)]

if (length(leading_edge_genes) < 2) {
  stop("Could not extract enough leading-edge genes.")
}

## -----------------------------
## 7. Match genes to expression matrix
## -----------------------------

expr_mat <- as.matrix(brca_rna_vst_or_logexpr_small)

matched_genes <- intersect(leading_edge_genes, rownames(expr_mat))

if (length(matched_genes) < 2) {
  ## Case-insensitive fallback
  rn_upper <- toupper(rownames(expr_mat))
  genes_upper <- toupper(leading_edge_genes)

  idx <- match(genes_upper, rn_upper)
  idx <- idx[!is.na(idx)]

  matched_genes <- unique(rownames(expr_mat)[idx])
}

if (length(matched_genes) < 2) {
  stop(
    "Not enough leading-edge genes matched expression matrix row names. ",
    "Matched genes: ",
    length(matched_genes)
  )
}

message("Matched leading-edge genes: ", length(matched_genes))

expr_sub <- expr_mat[matched_genes, , drop = FALSE]

## -----------------------------
## 8. Match samples and group by ER status
## -----------------------------
## No clustering is used.
## Samples are ordered by ER status, then sample ID.

metadata <- brca_rna_metadata_tumor_normal

metadata[[sample_col]] <- as.character(metadata[[sample_col]])
metadata[[er_col]] <- clean_er_status(metadata[[er_col]])

metadata <- metadata[!is.na(metadata[[er_col]]), , drop = FALSE]

metadata[[er_col]] <- factor(
  metadata[[er_col]],
  levels = c("ER Positive", "ER Negative")
)

common_samples <- intersect(colnames(expr_sub), metadata[[sample_col]])

if (length(common_samples) < 2) {
  stop(
    "Could not match expression columns with metadata sample IDs. ",
    "Check that colnames(expression_data) match metadata$", sample_col, "."
  )
}

metadata_sub <- metadata[match(common_samples, metadata[[sample_col]]), , drop = FALSE]
metadata_sub <- metadata_sub[!is.na(metadata_sub[[er_col]]), , drop = FALSE]

sample_order <- order(metadata_sub[[er_col]], metadata_sub[[sample_col]])

metadata_sub <- metadata_sub[sample_order, , drop = FALSE]
expr_sub <- expr_sub[, metadata_sub[[sample_col]], drop = FALSE]

message("\nSamples per ER group:")
print(table(metadata_sub[[er_col]]))

## -----------------------------
## 9. Scale expression by gene
## -----------------------------
## Z-score by row:
## blue = lower expression
## white = average expression
## red = higher expression

expr_z <- t(scale(t(expr_sub)))
expr_z[is.na(expr_z)] <- 0
expr_z[is.infinite(expr_z)] <- 0

## Clip extreme values to make the color scale readable.
expr_z[expr_z > 2.5] <- 2.5
expr_z[expr_z < -2.5] <- -2.5

## -----------------------------
## 10. Build ER annotation bar
## -----------------------------

annotation_col <- data.frame(
  ER_Status = metadata_sub[[er_col]],
  row.names = metadata_sub[[sample_col]]
)

annotation_colors <- list(
  ER_Status = c(
    "ER Positive" = "#D55E00",
    "ER Negative" = "#0072B2"
  )
)

## -----------------------------
## 11. Blue-white-red heatmap
## -----------------------------

heat_colors <- grDevices::colorRampPalette(
  c("blue", "white", "red")
)(101)

breaks <- seq(-2.5, 2.5, length.out = length(heat_colors) + 1)

jpeg(
  filename = target_file,
  width = 3600,
  height = 2700,
  res = 300,
  quality = 95
)

pheatmap::pheatmap(
  mat = expr_z,

  ## No clustering
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  treeheight_col = 0,

  ## Color scale
  color = heat_colors,
  breaks = breaks,

  ## ER positive / ER negative annotation
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,

  ## Legends
  legend = TRUE,
  annotation_legend = TRUE,

  ## Labels
  show_rownames = TRUE,
  show_colnames = FALSE,

  ## Text sizes
  fontsize = 9,
  fontsize_row = 7,
  fontsize_col = 5,

  ## Layout
  border_color = NA,
  main = paste0(
    target_pathway,
    "\nLeading-edge genes | ER Positive vs ER Negative"
  )
)

dev.off()

## -----------------------------
## 12. Final check
## -----------------------------

if (!file.exists(target_file)) {
  stop("Heatmap file was not created.")
}

message("\nDone. Heatmap saved at:")
message(normalizePath(target_file))
