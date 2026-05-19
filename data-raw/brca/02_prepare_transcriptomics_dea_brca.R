## 02_prepare_transcriptomics_dea_brca.R
##
## Purpose:
##   Prepare real TCGA-BRCA RNA-seq differential-expression objects from the
##   UCSC Xena HiSeqV2 log2-normalized expression matrix.
##
## Inputs:
##   - inst/extdata/brca/raw_xena/TCGA.BRCA.sampleMap_HiSeqV2/HiSeqV2
##   - data/brca_metadata.rda
##
## Outputs:
##   Package objects:
##   - brca_rna_dea_er_pos_vs_er_neg
##   - brca_rna_dea_tumor_vs_normal
##   - brca_rna_expr_er_shared_filtered
##   - brca_rna_expr_tumor_normal_filtered
##   - brca_rna_metadata_er_shared
##   - brca_rna_metadata_tumor_normal
##   - brca_rna_vst_or_logexpr_small
##
##   External files:
##   - inst/extdata/brca/intermediate/DEA_RNAseq_limma_ERpositive_vs_ERnegative.tsv
##   - inst/extdata/brca/intermediate/DEA_RNAseq_limma_Tumor_vs_Normal.tsv
##   - inst/extdata/brca/intermediate/expr_RNAseq_ERpositive_vs_ERnegative_filtered.rds
##   - inst/extdata/brca/intermediate/expr_RNAseq_Tumor_vs_Normal_filtered.rds
##
## Method:
##   limma on Xena log2-normalized RNA-seq expression with robust = TRUE and
##   trend = TRUE.

options(stringsAsFactors = FALSE)

min_expr <- 1
min_prop_samples <- 0.10
small_matrix_n_genes <- 1000
max_full_matrix_data_mb <- 25

project_file <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    here::here(...)
  } else {
    file.path(...)
  }
}

raw_expr_file <- project_file(
  "inst", "extdata", "brca", "raw_xena",
  "TCGA.BRCA.sampleMap_HiSeqV2", "HiSeqV2"
)
metadata_file <- project_file("data", "brca_metadata.rda")
intermediate_dir <- project_file("inst", "extdata", "brca", "intermediate")

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "Required input file is missing for BRCA RNA-seq DEA: ", label, "\n",
      "Expected location: ", normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }
}

require_file(raw_expr_file, "TCGA.BRCA.sampleMap_HiSeqV2/HiSeqV2")
require_file(metadata_file, "data/brca_metadata.rda")

if (!requireNamespace("limma", quietly = TRUE)) {
  stop(
    "Package \"limma\" is required for TCGA-BRCA RNA-seq DEA.",
    call. = FALSE
  )
}

if (!requireNamespace("statmod", quietly = TRUE)) {
  stop(
    "Package \"statmod\" is required because limma::eBayes() is run with robust = TRUE.",
    call. = FALSE
  )
}

if (!requireNamespace("usethis", quietly = TRUE)) {
  stop(
    "Package \"usethis\" is required to save BRCA RNA-seq package data.",
    call. = FALSE
  )
}

has_value <- function(x) {
  y <- trimws(as.character(x))
  !is.na(y) &
    nzchar(y) &
    !(tolower(y) %in% c(
      "na", "n/a", "nan", "null", "none", "unknown", "not available",
      "not reported", "not applicable", "[not available]",
      "[not applicable]", "[unknown]", "--"
    ))
}

clean_tcga_barcode <- function(x, level = c("sample", "patient")) {
  level <- match.arg(level)
  y <- toupper(trimws(as.character(x)))
  y <- gsub("\\.", "-", y)
  y[!has_value(y)] <- NA_character_

  if (identical(level, "patient")) {
    return(ifelse(!is.na(y) & nchar(y) >= 12L, substr(y, 1L, 12L), NA_character_))
  }

  ifelse(
    !is.na(y) & nchar(y) >= 16L,
    substr(y, 1L, 16L),
    ifelse(!is.na(y) & nchar(y) >= 15L, substr(y, 1L, 15L), NA_character_)
  )
}

tcga_sample_code <- function(x) {
  y <- clean_tcga_barcode(x, level = "sample")
  ifelse(!is.na(y) & nchar(y) >= 15L, substr(y, 14L, 15L), NA_character_)
}

normalize_sample_code <- function(x) {
  y <- trimws(as.character(x))
  y[!has_value(y)] <- NA_character_
  ifelse(!is.na(y) & nchar(y) == 1L, paste0("0", y), y)
}

row_variance <- function(x) {
  apply(x, 1L, stats::var, na.rm = TRUE)
}

collapse_duplicate_rows_by_mean <- function(mat, row_ids) {
  keep <- has_value(row_ids)
  mat <- mat[keep, , drop = FALSE]
  row_ids <- trimws(as.character(row_ids[keep]))

  unique_ids <- unique(row_ids)
  collapsed <- vapply(
    unique_ids,
    function(id) {
      colMeans(mat[row_ids == id, , drop = FALSE], na.rm = TRUE)
    },
    numeric(ncol(mat))
  )

  collapsed <- t(collapsed)
  rownames(collapsed) <- unique_ids
  colnames(collapsed) <- colnames(mat)
  collapsed
}

collapse_duplicate_columns_by_mean <- function(mat, col_ids) {
  keep <- has_value(col_ids)
  mat <- mat[, keep, drop = FALSE]
  col_ids <- as.character(col_ids[keep])

  unique_ids <- unique(col_ids)
  collapsed <- vapply(
    unique_ids,
    function(id) {
      rowMeans(mat[, col_ids == id, drop = FALSE], na.rm = TRUE)
    },
    numeric(nrow(mat))
  )

  rownames(collapsed) <- rownames(mat)
  colnames(collapsed) <- unique_ids
  collapsed
}

filter_expression_matrix <- function(expr, min_expr, min_prop_samples) {
  expressed <- rowMeans(expr >= min_expr, na.rm = TRUE) >= min_prop_samples
  vars <- row_variance(expr)
  variable <- is.finite(vars) & vars > 0
  expr[expressed & variable, , drop = FALSE]
}

format_limma_results <- function(top_table, comparison, reference, contrast) {
  out <- data.frame(
    gene_symbol = rownames(top_table),
    top_table,
    row.names = NULL,
    check.names = FALSE
  )

  out$gene_id <- out$gene_symbol
  out$padj <- out$adj.P.Val
  out$stat <- out$t
  out$comparison <- comparison
  out$reference <- reference
  out$contrast <- contrast
  out$method <- "limma_trend_robust_on_xena_log2_expression"

  out
}

fit_limma_contrast <- function(expr, group, levels, contrast_expr, comparison, reference) {
  group <- factor(group, levels = levels)
  group_counts <- table(group)

  if (any(group_counts < 2L)) {
    stop(
      "Cannot fit ", comparison, ": each group must have at least 2 samples. ",
      "Observed counts: ", paste(names(group_counts), group_counts, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }

  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit <- limma::lmFit(expr, design)
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_expr, levels = design)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  top_table <- limma::topTable(fit, coef = 1L, number = Inf, sort.by = "P")
  format_limma_results(
    top_table = top_table,
    comparison = comparison,
    reference = reference,
    contrast = contrast_expr
  )
}

compressed_rds_size_mb <- function(object) {
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)
  saveRDS(object, temp_file, compress = "xz")
  file.info(temp_file)$size / 1024^2
}

load(metadata_file)
if (!exists("brca_metadata")) {
  stop(
    "data/brca_metadata.rda must contain an object named brca_metadata.",
    call. = FALSE
  )
}

required_metadata_cols <- c(
  "sample16", "patient", "sample_code", "tumor_normal", "ER_group"
)
missing_metadata_cols <- setdiff(required_metadata_cols, names(brca_metadata))
if (length(missing_metadata_cols) > 0L) {
  stop(
    "brca_metadata is missing required columns for RNA-seq DEA: ",
    paste(missing_metadata_cols, collapse = ", "),
    call. = FALSE
  )
}

message("Reading TCGA-BRCA Xena HiSeqV2 expression matrix...")
rna_raw <- utils::read.delim(
  raw_expr_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE
)

if (ncol(rna_raw) < 3L) {
  stop(
    "RNA-seq expression file must contain a gene column and at least two sample columns.",
    call. = FALSE
  )
}

gene_symbol <- trimws(as.character(rna_raw[[1L]]))
gene_symbol <- sub("\\|.*$", "", gene_symbol)

rna_expr <- as.matrix(rna_raw[, -1L, drop = FALSE])
mode(rna_expr) <- "numeric"

if (all(is.na(rna_expr))) {
  stop(
    "RNA-seq expression values could not be converted to numeric values.",
    call. = FALSE
  )
}

rna_expr <- collapse_duplicate_rows_by_mean(rna_expr, gene_symbol)

sample16 <- clean_tcga_barcode(colnames(rna_expr), level = "sample")
rna_expr <- collapse_duplicate_columns_by_mean(rna_expr, sample16)

message(
  "Expression matrix after duplicate collapse: ",
  nrow(rna_expr), " genes x ", ncol(rna_expr), " samples"
)

brca_metadata$sample_code <- normalize_sample_code(brca_metadata$sample_code)

if ("RNA_genomic_id" %in% names(brca_metadata)) {
  brca_metadata$rna_sample16 <- clean_tcga_barcode(
    brca_metadata$RNA_genomic_id,
    level = "sample"
  )
  missing_rna_id <- !has_value(brca_metadata$RNA_genomic_id)
  brca_metadata$rna_sample16[missing_rna_id] <- clean_tcga_barcode(
    brca_metadata$sample16[missing_rna_id],
    level = "sample"
  )
} else {
  warning(
    "RNA_genomic_id is missing from brca_metadata; matching RNA samples by sample16.",
    call. = FALSE
  )
  brca_metadata$rna_sample16 <- clean_tcga_barcode(brca_metadata$sample16, level = "sample")
}

brca_metadata <- brca_metadata[has_value(brca_metadata$rna_sample16), , drop = FALSE]
brca_metadata <- brca_metadata[!duplicated(brca_metadata$rna_sample16), , drop = FALSE]

make_matched_metadata <- function(metadata, expr, sample_ids) {
  metadata <- metadata[metadata$rna_sample16 %in% sample_ids, , drop = FALSE]
  metadata <- metadata[match(sample_ids, metadata$rna_sample16), , drop = FALSE]
  rownames(metadata) <- metadata$rna_sample16
  metadata
}

## Comparison A: ER_positive vs ER_negative, primary tumor only.
er_keep <- brca_metadata$sample_code == "01" &
  brca_metadata$ER_group %in% c("ER_negative", "ER_positive") &
  brca_metadata$rna_sample16 %in% colnames(rna_expr)

if ("RPPA_genomic_id" %in% names(brca_metadata)) {
  er_keep <- er_keep & has_value(brca_metadata$RPPA_genomic_id)
} else {
  warning(
    "RPPA_genomic_id is missing from brca_metadata; ER RNA subset cannot be restricted to RNA/RPPA-shared samples.",
    call. = FALSE
  )
}

er_sample_ids <- brca_metadata$rna_sample16[er_keep]
if (length(er_sample_ids) == 0L) {
  stop(
    "No RNA-seq samples are available for ER_positive vs ER_negative after applying primary tumor, ER_group, RNA, and RPPA-sharing filters.",
    call. = FALSE
  )
}

brca_rna_metadata_er_shared <- make_matched_metadata(
  metadata = brca_metadata[er_keep, , drop = FALSE],
  expr = rna_expr,
  sample_ids = er_sample_ids
)

brca_rna_expr_er_shared_filtered <- rna_expr[, rownames(brca_rna_metadata_er_shared), drop = FALSE]
brca_rna_expr_er_shared_filtered <- filter_expression_matrix(
  brca_rna_expr_er_shared_filtered,
  min_expr = min_expr,
  min_prop_samples = min_prop_samples
)

if (nrow(brca_rna_expr_er_shared_filtered) == 0L) {
  stop(
    "No genes remain for ER_positive vs ER_negative after expression and variance filtering.",
    call. = FALSE
  )
}

brca_rna_dea_er_pos_vs_er_neg <- fit_limma_contrast(
  expr = brca_rna_expr_er_shared_filtered,
  group = brca_rna_metadata_er_shared$ER_group,
  levels = c("ER_negative", "ER_positive"),
  contrast_expr = "ER_positive - ER_negative",
  comparison = "ER_positive_vs_ER_negative",
  reference = "ER_negative"
)

## Comparison B: Tumor vs Normal, RNA-seq only.
tumor_normal_group <- brca_metadata$tumor_normal
tumor_normal_group[is.na(tumor_normal_group) & brca_metadata$sample_code %in% sprintf("%02d", 1:9)] <- "Tumor"
tumor_normal_group[is.na(tumor_normal_group) & brca_metadata$sample_code %in% sprintf("%02d", 10:19)] <- "Normal"

tn_keep <- tumor_normal_group %in% c("Normal", "Tumor") &
  brca_metadata$rna_sample16 %in% colnames(rna_expr)

tn_sample_ids <- brca_metadata$rna_sample16[tn_keep]
if (length(tn_sample_ids) == 0L) {
  stop(
    "No RNA-seq samples are available for Tumor vs Normal after matching brca_metadata to the expression matrix.",
    call. = FALSE
  )
}

brca_rna_metadata_tumor_normal <- make_matched_metadata(
  metadata = brca_metadata[tn_keep, , drop = FALSE],
  expr = rna_expr,
  sample_ids = tn_sample_ids
)
brca_rna_metadata_tumor_normal$rna_tumor_normal <- tumor_normal_group[tn_keep][
  match(rownames(brca_rna_metadata_tumor_normal), brca_metadata$rna_sample16[tn_keep])
]

brca_rna_expr_tumor_normal_filtered <- rna_expr[, rownames(brca_rna_metadata_tumor_normal), drop = FALSE]
brca_rna_expr_tumor_normal_filtered <- filter_expression_matrix(
  brca_rna_expr_tumor_normal_filtered,
  min_expr = min_expr,
  min_prop_samples = min_prop_samples
)

if (nrow(brca_rna_expr_tumor_normal_filtered) == 0L) {
  stop(
    "No genes remain for Tumor vs Normal after expression and variance filtering.",
    call. = FALSE
  )
}

brca_rna_dea_tumor_vs_normal <- fit_limma_contrast(
  expr = brca_rna_expr_tumor_normal_filtered,
  group = brca_rna_metadata_tumor_normal$rna_tumor_normal,
  levels = c("Normal", "Tumor"),
  contrast_expr = "Tumor - Normal",
  comparison = "Tumor_vs_Normal",
  reference = "Normal"
)

er_vars <- row_variance(brca_rna_expr_er_shared_filtered)
top_genes <- names(sort(er_vars, decreasing = TRUE))[seq_len(min(small_matrix_n_genes, length(er_vars)))]
brca_rna_vst_or_logexpr_small <- brca_rna_expr_er_shared_filtered[top_genes, , drop = FALSE]

dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.table(
  brca_rna_dea_er_pos_vs_er_neg,
  file = file.path(intermediate_dir, "DEA_RNAseq_limma_ERpositive_vs_ERnegative.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

utils::write.table(
  brca_rna_dea_tumor_vs_normal,
  file = file.path(intermediate_dir, "DEA_RNAseq_limma_Tumor_vs_Normal.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

saveRDS(
  brca_rna_expr_er_shared_filtered,
  file = file.path(intermediate_dir, "expr_RNAseq_ERpositive_vs_ERnegative_filtered.rds"),
  compress = "xz"
)

saveRDS(
  brca_rna_expr_tumor_normal_filtered,
  file = file.path(intermediate_dir, "expr_RNAseq_Tumor_vs_Normal_filtered.rds"),
  compress = "xz"
)

er_matrix_mb <- compressed_rds_size_mb(brca_rna_expr_er_shared_filtered)
tn_matrix_mb <- compressed_rds_size_mb(brca_rna_expr_tumor_normal_filtered)
save_full_matrices_to_data <- er_matrix_mb <= max_full_matrix_data_mb &&
  tn_matrix_mb <= max_full_matrix_data_mb

if (save_full_matrices_to_data) {
  usethis::use_data(
    brca_rna_dea_er_pos_vs_er_neg,
    brca_rna_dea_tumor_vs_normal,
    brca_rna_expr_er_shared_filtered,
    brca_rna_expr_tumor_normal_filtered,
    brca_rna_metadata_er_shared,
    brca_rna_metadata_tumor_normal,
    brca_rna_vst_or_logexpr_small,
    compress = "xz",
    overwrite = TRUE
  )
} else {
  warning(
    "Full filtered RNA matrices are larger than ",
    max_full_matrix_data_mb,
    " MB compressed and were not saved into data/. ",
    "They were written as RDS files under inst/extdata/brca/intermediate/.",
    call. = FALSE
  )

  usethis::use_data(
    brca_rna_dea_er_pos_vs_er_neg,
    brca_rna_dea_tumor_vs_normal,
    brca_rna_metadata_er_shared,
    brca_rna_metadata_tumor_normal,
    brca_rna_vst_or_logexpr_small,
    compress = "xz",
    overwrite = TRUE
  )
}

message("Saved RNA-seq DEA TSV files to inst/extdata/brca/intermediate/.")
message("Saved filtered RNA-seq expression RDS files to inst/extdata/brca/intermediate/.")
message("Saved package RNA-seq objects to data/.")
message(
  "Compressed full matrix sizes: ER shared = ", round(er_matrix_mb, 2),
  " MB; Tumor/Normal = ", round(tn_matrix_mb, 2), " MB."
)
