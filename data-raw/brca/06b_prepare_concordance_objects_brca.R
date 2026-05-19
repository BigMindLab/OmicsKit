## 06b_prepare_concordance_objects_brca.R
##
## Purpose:
##   Prepare real RNA-RPPA concordance objects for Omics_layers_BRCA.Rmd.
##
## Inputs:
##   - data/brca_rna_dea_er_pos_vs_er_neg.rda
##   - data/brca_rppa_dea_gene_er_pos_vs_er_neg.rda
##   - data/brca_rna_expr_er_shared_filtered.rda or
##     inst/extdata/brca/intermediate/expr_RNAseq_ERpositive_vs_ERnegative_filtered.rds
##   - data/brca_rppa_expr_gene_er_shared.rda or
##     inst/extdata/brca/intermediate/expr_RPPA_gene_ERpositive_vs_ERnegative.rds
##
## Outputs:
##   Package objects:
##   - brca_concordance_rna_rppa_er
##   - brca_crosslayercorr_rna_rppa_er
##
##   External files:
##   - inst/extdata/brca/intermediate/Concordance_RNAseq_RPPA_ERpositive_vs_ERnegative.tsv
##   - inst/extdata/brca/intermediate/CrossLayerCorr_RNAseq_RPPA_samples.tsv
##
## Notes:
##   This script does not rerun RNA or RPPA DEA. It only loads already prepared
##   DEA and expression objects.

options(stringsAsFactors = FALSE)

project_file <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    here::here(...)
  } else {
    file.path(...)
  }
}

data_dir <- project_file("data")
intermediate_dir <- project_file("inst", "extdata", "brca", "intermediate")

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "Required input is missing for BRCA RNA-RPPA concordance: ", label, "\n",
      "Expected location: ", normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }
}

load_rda_object <- function(object_name, path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)

  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop(
      "File does not contain expected object `", object_name, "`: ",
      normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }

  get(object_name, envir = env, inherits = FALSE)
}

find_rda_object <- function(object_name) {
  preferred_path <- file.path(data_dir, paste0(object_name, ".rda"))

  if (file.exists(preferred_path)) {
    return(load_rda_object(object_name, preferred_path))
  }

  all_rda <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  for (path in all_rda) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (exists(object_name, envir = env, inherits = FALSE)) {
      return(get(object_name, envir = env, inherits = FALSE))
    }
  }

  stop(
    "Could not find required object `", object_name, "` in data/*.rda.",
    call. = FALSE
  )
}

load_matrix_object <- function(object_name, rds_path) {
  preferred_path <- file.path(data_dir, paste0(object_name, ".rda"))

  if (file.exists(preferred_path)) {
    return(load_rda_object(object_name, preferred_path))
  }

  all_rda <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  for (path in all_rda) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (exists(object_name, envir = env, inherits = FALSE)) {
      return(get(object_name, envir = env, inherits = FALSE))
    }
  }

  require_file(rds_path, paste0(object_name, " RDS fallback"))
  readRDS(rds_path)
}

first_existing_column <- function(data, candidates, label) {
  found <- candidates[candidates %in% names(data)]
  if (length(found) == 0L) {
    stop(
      "Could not find ", label, " column. Expected one of: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }

  found[1L]
}

as_concordance_de <- function(data, layer) {
  data <- as.data.frame(data)

  gene_col <- first_existing_column(
    data,
    candidates = c("gene", "gene_symbol", "gene_id"),
    label = paste(layer, "gene")
  )
  logfc_col <- first_existing_column(
    data,
    candidates = c("logFC", "log2FoldChange"),
    label = paste(layer, "log fold-change")
  )
  padj_col <- first_existing_column(
    data,
    candidates = c("padj", "adj.P.Val", "FDR", "q_value"),
    label = paste(layer, "adjusted p-value")
  )

  out <- data.frame(
    gene = as.character(data[[gene_col]]),
    logFC = suppressWarnings(as.numeric(data[[logfc_col]])),
    padj = suppressWarnings(as.numeric(data[[padj_col]])),
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$gene) & nzchar(trimws(out$gene)), , drop = FALSE]
  out <- out[!is.na(out$logFC) & !is.na(out$padj), , drop = FALSE]

  if (anyDuplicated(out$gene) > 0L) {
    stop(
      layer,
      " DEA table must contain one row per gene before concordance analysis.",
      call. = FALSE
    )
  }

  out
}

standardize_matrix <- function(mat, label) {
  mat <- as.matrix(mat)

  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop(label, " matrix must have row names and column names.", call. = FALSE)
  }

  storage.mode(mat) <- "numeric"

  if (all(is.na(mat))) {
    stop(label, " matrix contains no numeric values.", call. = FALSE)
  }

  mat
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop(
    "Package \"dplyr\" is required because concordanceDE() uses dplyr::case_when().",
    call. = FALSE
  )
}

source(project_file("R", "concordanceDE.R"))
source(project_file("R", "crossLayerCorr.R"))

brca_rna_dea_er_pos_vs_er_neg <- find_rda_object("brca_rna_dea_er_pos_vs_er_neg")
brca_rppa_dea_gene_er_pos_vs_er_neg <- find_rda_object("brca_rppa_dea_gene_er_pos_vs_er_neg")

brca_rna_expr_er_shared_filtered <- load_matrix_object(
  "brca_rna_expr_er_shared_filtered",
  file.path(intermediate_dir, "expr_RNAseq_ERpositive_vs_ERnegative_filtered.rds")
)

brca_rppa_expr_gene_er_shared <- load_matrix_object(
  "brca_rppa_expr_gene_er_shared",
  file.path(intermediate_dir, "expr_RPPA_gene_ERpositive_vs_ERnegative.rds")
)

rna_de <- as_concordance_de(brca_rna_dea_er_pos_vs_er_neg, "RNA")
rppa_de <- as_concordance_de(brca_rppa_dea_gene_er_pos_vs_er_neg, "RPPA")

brca_concordance_rna_rppa_er <- concordanceDE(
  de_x = rna_de,
  de_y = rppa_de,
  gene_col = "gene",
  logfc_col = "logFC",
  padj_col = "padj",
  padj_threshold = 0.05,
  logfc_threshold = c(1, 0.20)
)

rna_expr <- standardize_matrix(brca_rna_expr_er_shared_filtered, "RNA")
rppa_expr <- standardize_matrix(brca_rppa_expr_gene_er_shared, "RPPA")

shared_genes <- intersect(rownames(rna_expr), rownames(rppa_expr))
if (length(shared_genes) < 2L) {
  stop(
    "At least two shared RNA/RPPA genes are required for crossLayerCorr().",
    call. = FALSE
  )
}

top_n <- min(length(shared_genes), 146L)

brca_crosslayercorr_rna_rppa_er <- crossLayerCorr(
  mat_x = rna_expr,
  mat_y = rppa_expr,
  method = "spearman",
  top_n = top_n,
  plot = FALSE
)

dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.table(
  brca_concordance_rna_rppa_er$table,
  file = file.path(intermediate_dir, "Concordance_RNAseq_RPPA_ERpositive_vs_ERnegative.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

crosslayer_table <- brca_crosslayercorr_rna_rppa_er$correlations
crosslayer_table$median_r <- brca_crosslayercorr_rna_rppa_er$median_r
crosslayer_table$n_shared_samples <- brca_crosslayercorr_rna_rppa_er$n_shared_samples
crosslayer_table$n_shared_features <- brca_crosslayercorr_rna_rppa_er$n_shared_features
crosslayer_table$n_features_used <- brca_crosslayercorr_rna_rppa_er$n_features_used

utils::write.table(
  crosslayer_table,
  file = file.path(intermediate_dir, "CrossLayerCorr_RNAseq_RPPA_samples.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

if (!requireNamespace("usethis", quietly = TRUE)) {
  stop(
    "Package \"usethis\" is required to save BRCA concordance package data.",
    call. = FALSE
  )
}

usethis::use_data(
  brca_concordance_rna_rppa_er,
  brca_crosslayercorr_rna_rppa_er,
  compress = "xz",
  overwrite = TRUE
)

message("Saved BRCA RNA-RPPA concordance objects to data/.")
message("Wrote concordance and cross-layer correlation TSV files to inst/extdata/brca/intermediate/.")
