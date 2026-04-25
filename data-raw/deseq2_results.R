## data-raw/deseq2_results.R
## Generates deseq2_results, , norm_counts, and vst_counts from TCGA-LUAD raw_counts and sampledata.
## Run once to regenerate the three .rda files in data/.
## Requires: DESeq2, dplyr, tibble

# =============================================================================
# DESeq2 pipeline - TCGA LUAD
# Comparison: Tumor vs Normal (sample_type column)
# Source: GDC Data Portal - TCGA-LUAD STAR Counts
# Outputs:
#   data/deseq2_results.rda  : DESeq2 results table (all genes, FDR filtered)
#   data/norm_counts.rda     : Normalized counts matrix (counts(dds, normalized=TRUE))
#   data/vst_counts.rda      : VST-transformed matrix (assay(vst(dds)))
# =============================================================================

library(DESeq2)
library(dplyr)
library(tibble)

data(raw_counts)
data(sampledata)

# Inspect sample_type levels
cat("sample_type levels:\n")
print(table(sampledata$sample_type))

# Reference = "normal" (baseline), comparison = "tumor"

sampledata$sample_type <- factor(
  sampledata$sample_type,
  levels = c("normal", "tumor")
)

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(raw_counts),
  colData   = sampledata,
  design    = ~ sample_type
)

# Filter: keep genes with >= 10 counts in at least 10 samples
keep <- rowSums(counts(dds) >= 10) >= 10
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# Run DESeq2
dds <- DESeq(dds)

# 1. deseq2 results
res <- results(
  dds,
  contrast  = c("sample_type", "tumor", "normal"),
  alpha     = 0.05
)

deseq2_results <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  arrange(padj)

cat("\ndeseq2_results:\n")
cat("Total genes with results:", nrow(deseq2_results), "\n")
cat("Significant (FDR < 0.05):", sum(deseq2_results$padj < 0.05), "\n")
cat("  Columns        :", paste(colnames(deseq2_results), collapse = ", "), "\n")

# 2. norm_counts
# Normalized counts: suitable for nice_VSB, detect_filter, add_annotations
# Counts are divided by DESeq2 size factors to correct for library size.
# Still in counts scale (not log-transformed).
norm_counts <- counts(dds, normalized = TRUE)

cat("\nnorm_counts:\n")
cat("  Dimensions     :", nrow(norm_counts), "genes x", ncol(norm_counts), "samples\n")
cat("  Value range    : [", round(min(norm_counts), 1), ",",
    round(max(norm_counts), 1), "]\n")

# 3. vst_counts
# Variance Stabilizing Transformation: suitable for nice_PCA, nice_UMAP, nice_tSNE.
# VST removes the mean-variance dependence of RNA-seq counts,
# placing all genes on a comparable log2-like scale for dimensionality
# reduction and sample-level clustering.
vst_counts <- assay(vst(dds, blind = TRUE))

cat("\nvst_counts:\n")
cat("  Dimensions     :", nrow(vst_counts), "genes x", ncol(vst_counts), "samples\n")
cat("  Value range    : [", round(min(vst_counts), 2), ",",
    round(max(vst_counts), 2), "]\n")

# Save
usethis::use_data(deseq2_results, compress = "xz", overwrite = TRUE)
usethis::use_data(norm_counts,    compress = "xz", overwrite = TRUE)
usethis::use_data(vst_counts,     compress = "xz", overwrite = TRUE)

message("\nDone. Saved:")
message("  data/deseq2_results.rda")
message("  data/norm_counts.rda")
message("  data/vst_counts.rda")
