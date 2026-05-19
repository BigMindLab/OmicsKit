## 03_prepare_proteomics_dea_brca.R
##
## Purpose:
##   Prepare real TCGA-BRCA RPPA/proteomics differential-expression objects
##   for ER_positive vs ER_negative primary tumors shared across RNA and RPPA.
##
## Inputs:
##   - inst/extdata/brca/raw_xena/[RPPA matrix file]
##   - data/brca_metadata.rda
##   - data/brca_rna_expr_er_shared_filtered.rda or
##     inst/extdata/brca/intermediate/expr_RNAseq_ERpositive_vs_ERnegative_filtered.rds
##
## Outputs:
##   Package objects:
##   - brca_rppa_dea_feature_er_pos_vs_er_neg
##   - brca_rppa_dea_gene_er_pos_vs_er_neg
##   - brca_rppa_expr_feature_er_shared
##   - brca_rppa_expr_gene_er_shared
##   - brca_rppa_feature_map
##
##   External files:
##   - inst/extdata/brca/intermediate/DEA_RPPA_feature_limma_ERpositive_vs_ERnegative.tsv
##   - inst/extdata/brca/intermediate/DEA_RPPA_gene_limma_ERpositive_vs_ERnegative.tsv
##   - inst/extdata/brca/intermediate/expr_RPPA_gene_ERpositive_vs_ERnegative.rds
##
## Method:
##   limma on the RPPA matrix with robust = TRUE.
##   Reference: ER_negative.
##   Contrast: ER_positive - ER_negative.
##   Interpretation: logFC > 0 means higher RPPA signal in ER_positive tumors.

options(stringsAsFactors = FALSE)

## If discovery is ambiguous, set this to a project-relative path under
## inst/extdata/brca/raw_xena/, for example:
## manual_rppa_matrix_file <- "inst/extdata/brca/raw_xena/TCGA.BRCA.sampleMap_RPPA/RPPA"
manual_rppa_matrix_file <- NULL

max_rppa_matrix_data_mb <- 25

project_file <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    here::here(...)
  } else {
    file.path(...)
  }
}

raw_dir <- project_file("inst", "extdata", "brca", "raw_xena")
metadata_file <- project_file("data", "brca_metadata.rda")
intermediate_dir <- project_file("inst", "extdata", "brca", "intermediate")

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "Required input file is missing for BRCA RPPA DEA: ", label, "\n",
      "Expected location: ", normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }
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

  sample_match <- regexpr(
    "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}",
    y,
    perl = TRUE
  )
  out <- rep(NA_character_, length(y))
  matched <- !is.na(y) & sample_match > 0L
  out[matched] <- substr(
    y[matched],
    sample_match[matched],
    sample_match[matched] + attr(sample_match, "match.length")[matched] - 1L
  )
  out
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

is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\\\]", path) || grepl("^[/\\\\]", path)
}

find_rppa_matrix_file <- function(raw_dir, manual_path = NULL) {
  if (!is.null(manual_path) && has_value(manual_path)) {
    if (is_absolute_path(manual_path)) {
      stop(
        "manual_rppa_matrix_file must be project-relative, not absolute: ",
        manual_path,
        call. = FALSE
      )
    }

    path <- project_file(manual_path)
    require_file(path, "manual_rppa_matrix_file")
    return(path)
  }

  known_candidates <- c(
    file.path(raw_dir, "TCGA.BRCA.sampleMap_RPPA", "RPPA"),
    file.path(raw_dir, "TCGA.BRCA.sampleMap_RPPA", "RPPA_RBN"),
    file.path(raw_dir, "TCGA.BRCA.sampleMap_RPPA_RBN", "RPPA_RBN"),
    file.path(raw_dir, "TCGA.BRCA.sampleMap_RPPA"),
    file.path(raw_dir, "TCGA.BRCA.sampleMap_RPPA_RBN")
  )
  known_candidates <- known_candidates[file.exists(known_candidates)]

  discovered <- character()
  if (dir.exists(raw_dir)) {
    all_files <- list.files(raw_dir, recursive = TRUE, full.names = TRUE, no.. = TRUE)
    all_files <- all_files[file.exists(all_files) & !dir.exists(all_files)]
    discovered <- all_files[grepl("rppa|protein|proteomic", all_files, ignore.case = TRUE)]
  }

  candidates <- unique(c(known_candidates, discovered))

  if (length(candidates) == 0L) {
    stop(
      "No RPPA matrix file was found for BRCA proteomics DEA.\n",
      "Place the real RPPA matrix under inst/extdata/brca/raw_xena/ or set ",
      "manual_rppa_matrix_file to a project-relative path in this script.",
      call. = FALSE
    )
  }

  if (length(candidates) > 1L) {
    stop(
      "Multiple possible RPPA matrix files were found. Set manual_rppa_matrix_file ",
      "to the intended project-relative file.\nCandidates:\n- ",
      paste(normalizePath(candidates, mustWork = FALSE), collapse = "\n- "),
      call. = FALSE
    )
  }

  candidates
}

read_xena_matrix <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
}

matrix_from_xena_table <- function(tbl) {
  if (ncol(tbl) < 3L) {
    stop(
      "RPPA matrix file must contain identifiers and at least two data columns.",
      call. = FALSE
    )
  }

  header_sample_ids <- clean_tcga_barcode(names(tbl)[-1L], level = "sample")
  first_col_sample_ids <- clean_tcga_barcode(tbl[[1L]], level = "sample")

  n_header_samples <- sum(has_value(header_sample_ids))
  n_first_col_samples <- sum(has_value(first_col_sample_ids))

  if (n_header_samples >= 2L && n_header_samples >= n_first_col_samples) {
    sample_col_idx <- which(c(FALSE, has_value(header_sample_ids)))
    annotation_col_idx <- setdiff(seq_along(tbl), sample_col_idx)
    feature_id_col <- annotation_col_idx[1L]

    feature_id <- trimws(as.character(tbl[[feature_id_col]]))
    mat <- as.matrix(tbl[, sample_col_idx, drop = FALSE])
    mode(mat) <- "numeric"
    rownames(mat) <- feature_id
    colnames(mat) <- clean_tcga_barcode(names(tbl)[sample_col_idx], level = "sample")
    attr(mat, "feature_annotation") <- tbl[, annotation_col_idx, drop = FALSE]
    return(mat)
  }

  if (n_first_col_samples >= 2L) {
    sample_ids <- first_col_sample_ids
    mat <- as.matrix(tbl[, -1L, drop = FALSE])
    mode(mat) <- "numeric"
    mat <- t(mat)
    rownames(mat) <- names(tbl)[-1L]
    colnames(mat) <- sample_ids
    return(mat)
  }

  stop(
    "Could not detect RPPA matrix orientation. Expected TCGA barcodes either ",
    "in column names or in the first column.",
    call. = FALSE
  )
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

row_variance <- function(x) {
  apply(x, 1L, stats::var, na.rm = TRUE)
}

remove_zero_variance_rows <- function(mat) {
  vars <- row_variance(mat)
  mat[is.finite(vars) & vars > 0, , drop = FALSE]
}

extract_single_explicit_gene <- function(feature_id) {
  tokens <- unlist(strsplit(feature_id, "[|;]"))
  tokens <- trimws(tokens)

  paren_tokens <- unlist(regmatches(
    feature_id,
    gregexpr("\\(([A-Z0-9-]{2,20})\\)", feature_id, perl = TRUE)
  ))
  paren_tokens <- gsub("^\\(|\\)$", "", paren_tokens)

  candidates <- unique(c(tokens, paren_tokens))
  candidates <- candidates[grepl("^[A-Z][A-Z0-9-]{1,19}$", candidates)]

  if (length(candidates) == 1L) {
    candidates
  } else {
    NA_character_
  }
}

normalize_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

find_column <- function(data, candidates) {
  data_names <- names(data)
  data_norm <- normalize_name(data_names)
  candidate_norm <- normalize_name(candidates)

  matched <- match(candidate_norm, data_norm)
  if (any(!is.na(matched))) {
    return(data_names[matched[which(!is.na(matched))[1L]]])
  }

  NULL
}

clean_single_gene <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[!has_value(y)] <- NA_character_
  y[!grepl("^[A-Z][A-Z0-9-]{1,19}$", y)] <- NA_character_
  y
}

make_feature_map <- function(feature_ids, feature_annotation = NULL) {
  feature_ids <- trimws(as.character(feature_ids))
  gene <- rep(NA_character_, length(feature_ids))
  mapping_source <- rep(NA_character_, length(feature_ids))

  protein_label <- feature_ids

  if (!is.null(feature_annotation) && nrow(feature_annotation) == length(feature_ids)) {
    gene_col <- find_column(
      feature_annotation,
      c("gene", "Gene", "genes", "Genes", "gene_symbol", "Gene Symbol",
        "Hugo_Symbol", "HUGO_SYMBOL", "HGNC", "hgnc_symbol")
    )
    protein_col <- find_column(
      feature_annotation,
      c("protein_label", "Protein", "protein", "antibody", "Antibody",
        "Composite.Element.REF", "feature", "Feature")
    )

    if (!is.null(gene_col)) {
      gene <- clean_single_gene(feature_annotation[[gene_col]])
      mapping_source[has_value(gene)] <- paste0("annotation_column:", gene_col)
    }

    if (!is.null(protein_col)) {
      protein_label <- trimws(as.character(feature_annotation[[protein_col]]))
      protein_label[!has_value(protein_label)] <- feature_ids[!has_value(protein_label)]
    }
  }

  parsed_gene <- vapply(feature_ids, extract_single_explicit_gene, character(1L))
  fill_from_feature_id <- !has_value(gene) & has_value(parsed_gene)
  gene[fill_from_feature_id] <- parsed_gene[fill_from_feature_id]
  mapping_source[fill_from_feature_id] <- "feature_id_explicit_gene_token"

  has_explicit_delimiter <- grepl("[|;]", feature_ids) |
    grepl("\\([A-Z0-9-]{2,20}\\)", feature_ids)
  has_ambiguous_mapping <- has_explicit_delimiter & !has_value(gene)

  feature_type <- ifelse(
    grepl("phospho|phosphoryl|_p[A-Z0-9]|-p[A-Z0-9]", feature_ids, ignore.case = TRUE),
    "phospho",
    "total_or_other"
  )

  mapping_class <- ifelse(
    has_value(gene),
    "single_explicit_gene",
    ifelse(has_ambiguous_mapping, "ambiguous_explicit_mapping", "unmapped")
  )

  data.frame(
    feature_id = feature_ids,
    protein_label = protein_label,
    protein_base = sub("[|;].*$", "", protein_label),
    gene = gene,
    feature_type = feature_type,
    mapping_class = mapping_class,
    mapping_source = mapping_source,
    eligible_for_rna_protein_concordance = has_value(gene),
    stringsAsFactors = FALSE
  )
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
  fit <- limma::eBayes(fit, robust = TRUE)

  top_table <- limma::topTable(fit, coef = 1L, number = Inf, sort.by = "P")
  out <- data.frame(
    id = rownames(top_table),
    top_table,
    row.names = NULL,
    check.names = FALSE
  )
  out$padj <- out$adj.P.Val
  out$stat <- out$t
  out$comparison <- comparison
  out$reference <- reference
  out$contrast <- contrast_expr
  out$method <- "limma_robust_on_rppa_matrix"
  out
}

compressed_rds_size_mb <- function(object) {
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file), add = TRUE)
  saveRDS(object, temp_file, compress = "xz")
  file.info(temp_file)$size / 1024^2
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

load_rna_er_shared_sample_ids <- function() {
  object_name <- "brca_rna_expr_er_shared_filtered"
  preferred_path <- project_file("data", paste0(object_name, ".rda"))

  if (file.exists(preferred_path)) {
    return(colnames(load_rda_object(object_name, preferred_path)))
  }

  data_dir <- project_file("data")
  all_rda <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  for (path in all_rda) {
    env <- new.env(parent = emptyenv())
    load(path, envir = env)
    if (exists(object_name, envir = env, inherits = FALSE)) {
      return(colnames(get(object_name, envir = env, inherits = FALSE)))
    }
  }

  rds_path <- file.path(
    intermediate_dir,
    "expr_RNAseq_ERpositive_vs_ERnegative_filtered.rds"
  )
  require_file(rds_path, "expr_RNAseq_ERpositive_vs_ERnegative_filtered.rds")
  colnames(readRDS(rds_path))
}

require_file(metadata_file, "data/brca_metadata.rda")

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Package \"limma\" is required for TCGA-BRCA RPPA DEA.", call. = FALSE)
}

if (!requireNamespace("statmod", quietly = TRUE)) {
  stop(
    "Package \"statmod\" is required because limma::eBayes() is run with robust = TRUE.",
    call. = FALSE
  )
}

if (!requireNamespace("usethis", quietly = TRUE)) {
  stop(
    "Package \"usethis\" is required to save BRCA RPPA package data.",
    call. = FALSE
  )
}

rppa_matrix_file <- find_rppa_matrix_file(raw_dir, manual_rppa_matrix_file)
message("Reading TCGA-BRCA RPPA matrix: ", rppa_matrix_file)

load(metadata_file)
if (!exists("brca_metadata")) {
  stop(
    "data/brca_metadata.rda must contain an object named brca_metadata.",
    call. = FALSE
  )
}

required_metadata_cols <- c(
  "sample16", "sample_code", "sample_type", "ER_group"
)
missing_metadata_cols <- setdiff(required_metadata_cols, names(brca_metadata))
if (length(missing_metadata_cols) > 0L) {
  stop(
    "brca_metadata is missing required columns for RPPA DEA: ",
    paste(missing_metadata_cols, collapse = ", "),
    call. = FALSE
  )
}

rppa_tbl <- read_xena_matrix(rppa_matrix_file)
rppa_expr <- matrix_from_xena_table(rppa_tbl)
raw_feature_annotation <- attr(rppa_expr, "feature_annotation")
raw_feature_map <- make_feature_map(rownames(rppa_expr), raw_feature_annotation)

if (all(is.na(rppa_expr))) {
  stop("RPPA matrix values could not be converted to numeric values.", call. = FALSE)
}

rppa_expr <- collapse_duplicate_rows_by_mean(rppa_expr, rownames(rppa_expr))
rppa_expr <- collapse_duplicate_columns_by_mean(
  rppa_expr,
  clean_tcga_barcode(colnames(rppa_expr), level = "sample")
)

message("RPPA matrix columns after TCGA barcode cleaning: ", ncol(rppa_expr))

brca_rppa_feature_map <- raw_feature_map[
  match(rownames(rppa_expr), raw_feature_map$feature_id),
  ,
  drop = FALSE
]
rownames(brca_rppa_feature_map) <- NULL

brca_metadata$sample16 <- clean_tcga_barcode(brca_metadata$sample16, level = "sample")
brca_metadata$sample_code <- tcga_sample_code(brca_metadata$sample16)
brca_metadata$rppa_sample16 <- brca_metadata$sample16
brca_metadata <- brca_metadata[has_value(brca_metadata$sample16), , drop = FALSE]
brca_metadata <- brca_metadata[!duplicated(brca_metadata$sample16), , drop = FALSE]

rna_er_sample_ids <- clean_tcga_barcode(load_rna_er_shared_sample_ids(), level = "sample")
rna_er_sample_ids <- rna_er_sample_ids[has_value(rna_er_sample_ids)]

message(
  "RPPA columns matching brca_metadata$sample16: ",
  sum(colnames(rppa_expr) %in% brca_metadata$sample16)
)
message(
  "RPPA columns shared with RNA ER expression samples: ",
  sum(colnames(rppa_expr) %in% rna_er_sample_ids)
)

er_keep <- brca_metadata$sample_type == "Primary Tumor" &
  brca_metadata$sample_code == "01" &
  brca_metadata$ER_group %in% c("ER_negative", "ER_positive") &
  has_value(brca_metadata$rppa_sample16) &
  brca_metadata$rppa_sample16 %in% rna_er_sample_ids &
  brca_metadata$rppa_sample16 %in% colnames(rppa_expr)

er_metadata <- brca_metadata[er_keep, , drop = FALSE]
er_metadata <- er_metadata[!duplicated(er_metadata$rppa_sample16), , drop = FALSE]

if (nrow(er_metadata) == 0L) {
  stop(
    "No RPPA samples are available for ER_positive vs ER_negative after applying Primary Tumor, sample_code == \"01\", ER_group, RPPA sample16 matching, and RNA+RPPA sharing filters.",
    call. = FALSE
  )
}

er_metadata <- er_metadata[order(er_metadata$ER_group, er_metadata$rppa_sample16), , drop = FALSE]
rownames(er_metadata) <- er_metadata$rppa_sample16

brca_rppa_expr_feature_er_shared <- rppa_expr[, rownames(er_metadata), drop = FALSE]
brca_rppa_expr_feature_er_shared <- remove_zero_variance_rows(brca_rppa_expr_feature_er_shared)

if (nrow(brca_rppa_expr_feature_er_shared) == 0L) {
  stop(
    "No RPPA features remain for ER_positive vs ER_negative after zero-variance filtering.",
    call. = FALSE
  )
}

feature_dea <- fit_limma_contrast(
  expr = brca_rppa_expr_feature_er_shared,
  group = er_metadata$ER_group,
  levels = c("ER_negative", "ER_positive"),
  contrast_expr = "ER_positive - ER_negative",
  comparison = "ER_positive_vs_ER_negative",
  reference = "ER_negative"
)

brca_rppa_dea_feature_er_pos_vs_er_neg <- merge(
  feature_dea,
  brca_rppa_feature_map,
  by.x = "id",
  by.y = "feature_id",
  all.x = TRUE,
  sort = FALSE
)
names(brca_rppa_dea_feature_er_pos_vs_er_neg)[
  names(brca_rppa_dea_feature_er_pos_vs_er_neg) == "id"
] <- "feature_id"

gene_features <- brca_rppa_feature_map$eligible_for_rna_protein_concordance &
  brca_rppa_feature_map$feature_id %in% rownames(brca_rppa_expr_feature_er_shared)

if (!any(gene_features)) {
  stop(
    "No RPPA features have explicit single-gene mappings. Gene-level RPPA DEA ",
    "cannot be created without inventing gene mappings.",
    call. = FALSE
  )
}

mapped_feature_ids <- brca_rppa_feature_map$feature_id[gene_features]
mapped_genes <- brca_rppa_feature_map$gene[gene_features]
mapped_expr <- brca_rppa_expr_feature_er_shared[mapped_feature_ids, , drop = FALSE]

brca_rppa_expr_gene_er_shared <- collapse_duplicate_rows_by_mean(
  mapped_expr,
  mapped_genes
)
brca_rppa_expr_gene_er_shared <- remove_zero_variance_rows(brca_rppa_expr_gene_er_shared)

if (nrow(brca_rppa_expr_gene_er_shared) == 0L) {
  stop(
    "No gene-level RPPA rows remain after collapsing mapped features and filtering zero-variance rows.",
    call. = FALSE
  )
}

gene_dea <- fit_limma_contrast(
  expr = brca_rppa_expr_gene_er_shared,
  group = er_metadata$ER_group,
  levels = c("ER_negative", "ER_positive"),
  contrast_expr = "ER_positive - ER_negative",
  comparison = "ER_positive_vs_ER_negative",
  reference = "ER_negative"
)

names(gene_dea)[names(gene_dea) == "id"] <- "gene"
brca_rppa_dea_gene_er_pos_vs_er_neg <- gene_dea

gene_matrix_mb <- compressed_rds_size_mb(brca_rppa_expr_gene_er_shared)
if (gene_matrix_mb > max_rppa_matrix_data_mb) {
  stop(
    "The gene-level RPPA matrix is unexpectedly large (",
    round(gene_matrix_mb, 2),
    " MB compressed). Review before saving it into data/.",
    call. = FALSE
  )
}

dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.table(
  brca_rppa_dea_feature_er_pos_vs_er_neg,
  file = file.path(intermediate_dir, "DEA_RPPA_feature_limma_ERpositive_vs_ERnegative.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

utils::write.table(
  brca_rppa_dea_gene_er_pos_vs_er_neg,
  file = file.path(intermediate_dir, "DEA_RPPA_gene_limma_ERpositive_vs_ERnegative.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

saveRDS(
  brca_rppa_expr_gene_er_shared,
  file = file.path(intermediate_dir, "expr_RPPA_gene_ERpositive_vs_ERnegative.rds"),
  compress = "xz"
)

usethis::use_data(
  brca_rppa_dea_feature_er_pos_vs_er_neg,
  brca_rppa_dea_gene_er_pos_vs_er_neg,
  brca_rppa_expr_feature_er_shared,
  brca_rppa_expr_gene_er_shared,
  brca_rppa_feature_map,
  compress = "xz",
  overwrite = TRUE
)

message("Saved RPPA DEA TSV files to inst/extdata/brca/intermediate/.")
message("Saved gene-level RPPA expression RDS to inst/extdata/brca/intermediate/.")
message("Saved RPPA package objects to data/.")
message("Compressed gene-level RPPA matrix size: ", round(gene_matrix_mb, 2), " MB.")
