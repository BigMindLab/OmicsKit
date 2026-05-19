## 05_import_gsea_pipeline_outputs_brca.R
##
## Purpose:
##   Import real TCGA-BRCA GSEA outputs generated externally by
##   BigMindLab/GSEA-Pipeline.
##
## Inputs:
##   Expected files in inst/extdata/brca/gsea_outputs/:
##   - GSEA_Tumor_vs_Normal_Hallmark.tsv
##   - GSEA_Tumor_vs_Normal_GO_BP.tsv
##   - GSEA_ERpositive_vs_ERnegative_Hallmark.tsv
##   - GSEA_ERpositive_vs_ERnegative_GO_BP.tsv
##
## Outputs:
##   Cleaned imported TSV files in inst/extdata/brca/intermediate/:
##   - brca_gsea_tumor_vs_normal_hallmark.tsv
##   - brca_gsea_tumor_vs_normal_go_bp.tsv
##   - brca_gsea_er_pos_vs_er_neg_hallmark.tsv
##   - brca_gsea_er_pos_vs_er_neg_go_bp.tsv
##   - brca_gsea_all_imported_cleaned.tsv
##
## Expected file locations:
##   - inst/extdata/brca/gsea_outputs/
##   - inst/extdata/brca/intermediate/
##
## Notes:
##   - This script imports already-generated real GSEA outputs only.
##   - It does not run GSEA-Pipeline.
##   - It does not download MSigDB or gene-set files.

options(stringsAsFactors = FALSE)

project_file <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    here::here(...)
  } else {
    file.path(getwd(), ...)
  }
}

gsea_dir <- project_file("inst", "extdata", "brca", "gsea_outputs")
intermediate_dir <- project_file("inst", "extdata", "brca", "intermediate")

expected_files <- data.frame(
  file = c(
    "GSEA_Tumor_vs_Normal_Hallmark.tsv",
    "GSEA_Tumor_vs_Normal_GO_BP.tsv",
    "GSEA_ERpositive_vs_ERnegative_Hallmark.tsv",
    "GSEA_ERpositive_vs_ERnegative_GO_BP.tsv"
  ),
  comparison = c(
    "Tumor_vs_Normal",
    "Tumor_vs_Normal",
    "ER_positive_vs_ER_negative",
    "ER_positive_vs_ER_negative"
  ),
  collection = c(
    "HALLMARK",
    "GO_BP",
    "HALLMARK",
    "GO_BP"
  ),
  output_file = c(
    "brca_gsea_tumor_vs_normal_hallmark.tsv",
    "brca_gsea_tumor_vs_normal_go_bp.tsv",
    "brca_gsea_er_pos_vs_er_neg_hallmark.tsv",
    "brca_gsea_er_pos_vs_er_neg_go_bp.tsv"
  ),
  stringsAsFactors = FALSE
)

required_gsea_columns <- c(
  "NAME",
  "SIZE",
  "ES",
  "NES",
  "NOM p-val",
  "FDR q-val",
  "FWER p-val",
  "RANK AT MAX",
  "LEADING EDGE"
)

known_optional_columns <- c(
  "GS<br> follow link to MSigDB",
  "GS DETAILS",
  "DESCRIPTION",
  "URL",
  "Comparison",
  "COMPARISON",
  "CORE ENRICHMENT",
  "CORE_ENRICHMENT",
  "core_enrichment",
  "LEADING_EDGE_GENES",
  "leading_edge_genes",
  "leadingEdgeGenes",
  "GENES",
  "genes",
  "gene_symbols",
  "MEMBERS",
  "members",
  "...12"
)

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "Required GSEA-Pipeline output is missing for BRCA pathway import: ",
      label, "\nExpected location: ", path,
      call. = FALSE
    )
  }
  invisible(path)
}

read_tsv <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
}

write_tsv <- function(data, path) {
  utils::write.table(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = ""
  )
  invisible(path)
}

validate_columns <- function(data, file_label) {
  observed <- names(data)
  allowed <- c(required_gsea_columns, known_optional_columns)
  missing <- setdiff(required_gsea_columns, observed)
  unexpected <- setdiff(observed, allowed)

  if (length(missing) > 0L || length(unexpected) > 0L) {
    stop(
      "Unexpected columns in ", file_label, ".\n",
      "Expected required columns:\n- ",
      paste(required_gsea_columns, collapse = "\n- "),
      "\nAllowed optional columns:\n- ",
      paste(known_optional_columns, collapse = "\n- "),
      "\nObserved columns:\n- ",
      paste(observed, collapse = "\n- "),
      if (length(missing) > 0L) {
        paste0("\nMissing required columns:\n- ", paste(missing, collapse = "\n- "))
      } else {
        ""
      },
      if (length(unexpected) > 0L) {
        paste0("\nUnexpected columns:\n- ", paste(unexpected, collapse = "\n- "))
      } else {
        ""
      },
      call. = FALSE
    )
  }

  invisible(TRUE)
}

as_numeric_checked <- function(x, col, file_label) {
  y <- suppressWarnings(as.numeric(x))
  bad <- !is.na(x) & nzchar(trimws(as.character(x))) & is.na(y)
  if (any(bad)) {
    stop(
      "Column `", col, "` in ", file_label,
      " contains non-numeric values. First offending row: ",
      which(bad)[1L],
      call. = FALSE
    )
  }
  y
}

parse_leading_edge_fraction <- function(x, key) {
  pattern <- paste0(key, "=([0-9.]+)%")
  out <- rep(NA_real_, length(x))
  hit <- regexec(pattern, x)
  pieces <- regmatches(x, hit)
  ok <- lengths(pieces) >= 2L
  out[ok] <- as.numeric(vapply(pieces[ok], `[`, character(1), 2L)) / 100
  out
}

standardize_gsea <- function(data, comparison, collection, file_label) {
  validate_columns(data, file_label)

  if ("...12" %in% names(data)) {
    data[["...12"]] <- NULL
  }

  numeric_cols <- c(
    "SIZE", "ES", "NES", "NOM p-val",
    "FDR q-val", "FWER p-val", "RANK AT MAX"
  )

  for (col in numeric_cols) {
    data[[col]] <- as_numeric_checked(data[[col]], col, file_label)
  }

  if (any(data[["FDR q-val"]] == 0, na.rm = TRUE)) {
    data[["FDR q-val"]][data[["FDR q-val"]] == 0] <- 0.001
  }

  data$COLLECTION <- collection
  data$COMPARISON <- comparison
  data$FDR <- data[["FDR q-val"]]
  data$Log10FDR <- -log10(data$FDR)
  data$tags <- parse_leading_edge_fraction(data[["LEADING EDGE"]], "tags")
  data$list <- parse_leading_edge_fraction(data[["LEADING EDGE"]], "list")
  data$signal <- parse_leading_edge_fraction(data[["LEADING EDGE"]], "signal")
  data$source_file <- basename(file_label)

  remove_cols <- intersect(
    c("GS<br> follow link to MSigDB", "GS DETAILS", "Comparison"),
    names(data)
  )
  data[remove_cols] <- NULL

  preferred <- c(
    "NAME", "COLLECTION", "COMPARISON", "SIZE", "ES", "NES", "NOM p-val",
    "FDR", "FWER p-val", "Log10FDR", "RANK AT MAX", "LEADING EDGE",
    "tags", "list", "signal", "source_file"
  )
  data[, c(intersect(preferred, names(data)), setdiff(names(data), preferred)), drop = FALSE]
}

require_file(gsea_dir, "inst/extdata/brca/gsea_outputs/")
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

imported_tables <- vector("list", nrow(expected_files))

for (i in seq_len(nrow(expected_files))) {
  spec <- expected_files[i, , drop = FALSE]
  input_path <- file.path(gsea_dir, spec$file)
  output_path <- file.path(intermediate_dir, spec$output_file)

  require_file(input_path, spec$file)
  raw <- read_tsv(input_path)

  cleaned <- standardize_gsea(
    data = raw,
    comparison = spec$comparison,
    collection = spec$collection,
    file_label = spec$file
  )

  write_tsv(cleaned, output_path)
  imported_tables[[i]] <- cleaned
}

brca_gsea_all_imported_cleaned <- do.call(rbind, imported_tables)
rownames(brca_gsea_all_imported_cleaned) <- NULL

write_tsv(
  brca_gsea_all_imported_cleaned,
  file.path(intermediate_dir, "brca_gsea_all_imported_cleaned.tsv")
)

message(
  "Imported and cleaned ", nrow(brca_gsea_all_imported_cleaned),
  " BRCA GSEA rows from ", nrow(expected_files), " files."
)
