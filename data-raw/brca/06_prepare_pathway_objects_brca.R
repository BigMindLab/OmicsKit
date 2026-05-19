## 06_prepare_pathway_objects_brca.R
##
## Purpose:
##   Convert imported real TCGA-BRCA GSEA-Pipeline outputs into compact
##   OmicsKit pathway-analysis example objects.
##
## Inputs:
##   Cleaned imported TSV files in inst/extdata/brca/intermediate/:
##   - brca_gsea_tumor_vs_normal_hallmark.tsv
##   - brca_gsea_tumor_vs_normal_go_bp.tsv
##   - brca_gsea_er_pos_vs_er_neg_hallmark.tsv
##   - brca_gsea_er_pos_vs_er_neg_go_bp.tsv
##
## Outputs:
##   Package objects in data/:
##   - brca_gsea_tumor_vs_normal_hallmark_go
##   - brca_gsea_er_pos_vs_er_neg_hallmark_go
##   - brca_pa_merged
##   - brca_pa_clustering
##
## Expected file locations:
##   - inst/extdata/brca/intermediate/
##   - data/
##
## Notes:
##   - This script does not run GSEA-Pipeline.
##   - This script does not download MSigDB.
##   - Gene-set membership for clustering is extracted only from real gene
##     membership/core-enrichment columns present in the imported TSVs.

options(stringsAsFactors = FALSE)

project_file <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    here::here(...)
  } else {
    file.path(getwd(), ...)
  }
}

intermediate_dir <- project_file("inst", "extdata", "brca", "intermediate")

expected_imported_files <- data.frame(
  file = c(
    "brca_gsea_tumor_vs_normal_hallmark.tsv",
    "brca_gsea_tumor_vs_normal_go_bp.tsv",
    "brca_gsea_er_pos_vs_er_neg_hallmark.tsv",
    "brca_gsea_er_pos_vs_er_neg_go_bp.tsv"
  ),
  object_group = c(
    "tumor_vs_normal",
    "tumor_vs_normal",
    "er_pos_vs_er_neg",
    "er_pos_vs_er_neg"
  ),
  stringsAsFactors = FALSE
)

required_clean_columns <- c(
  "NAME",
  "COLLECTION",
  "COMPARISON",
  "SIZE",
  "ES",
  "NES",
  "NOM p-val",
  "FDR",
  "FWER p-val",
  "Log10FDR",
  "RANK AT MAX",
  "LEADING EDGE"
)

gene_membership_columns <- c(
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
  "members"
)

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "Required imported BRCA pathway file is missing: ", label, "\n",
      "Expected location: ", path, "\n",
      "Run data-raw/brca/05_import_gsea_pipeline_outputs_brca.R first.",
      call. = FALSE
    )
  }
  invisible(path)
}

require_package <- function(package, why) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package \"", package, "\" is required ", why, ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
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

validate_clean_columns <- function(data, file_label) {
  observed <- names(data)
  missing <- setdiff(required_clean_columns, observed)

  if (length(missing) > 0L) {
    stop(
      "Imported BRCA pathway file has unexpected columns: ", file_label, ".\n",
      "Expected columns include:\n- ",
      paste(required_clean_columns, collapse = "\n- "),
      "\nObserved columns:\n- ",
      paste(observed, collapse = "\n- "),
      "\nMissing required columns:\n- ",
      paste(missing, collapse = "\n- "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

compact_pa_table <- function(data) {
  keep <- unique(c(
    required_clean_columns,
    "tags",
    "list",
    "signal",
    "source_file",
    intersect(gene_membership_columns, names(data))
  ))

  data[, intersect(keep, names(data)), drop = FALSE]
}

extract_gene_vector <- function(x) {
  if (is.na(x) || !nzchar(trimws(as.character(x)))) {
    return(character())
  }

  genes <- unlist(strsplit(as.character(x), "[,;/|[:space:]]+"))
  genes <- trimws(genes)
  genes <- genes[nzchar(genes)]
  unique(genes)
}

make_geneset_list <- function(data) {
  gene_col <- intersect(gene_membership_columns, names(data))[1L]

  if (is.na(gene_col)) {
    stop(
      "Cannot create brca_pa_clustering because none of the imported GSEA ",
      "tables contains a real gene membership/core-enrichment column.\n",
      "Expected one of:\n- ",
      paste(gene_membership_columns, collapse = "\n- "),
      "\nDo not download MSigDB or invent memberships inside OmicsKit. ",
      "Export a gene membership column from BigMindLab/GSEA-Pipeline and rerun ",
      "the import script.",
      call. = FALSE
    )
  }

  split_names <- paste(data$COMPARISON, data$COLLECTION, data$NAME, sep = "::")
  gene_sets <- stats::setNames(
    lapply(data[[gene_col]], extract_gene_vector),
    split_names
  )

  gene_sets <- gene_sets[lengths(gene_sets) > 0L]

  if (length(gene_sets) < 3L) {
    stop(
      "At least three real gene sets with gene memberships are required for ",
      "brca_pa_clustering. Observed: ", length(gene_sets),
      call. = FALSE
    )
  }

  gene_sets
}

source(project_file("R", "doclust_PA.R"))

require_file(intermediate_dir, "inst/extdata/brca/intermediate/")

imported <- vector("list", nrow(expected_imported_files))
for (i in seq_len(nrow(expected_imported_files))) {
  file_name <- expected_imported_files$file[[i]]
  path <- file.path(intermediate_dir, file_name)
  require_file(path, file_name)

  data <- read_tsv(path)
  validate_clean_columns(data, file_name)

  imported[[i]] <- compact_pa_table(data)
}

names(imported) <- expected_imported_files$file

brca_gsea_tumor_vs_normal_hallmark_go <- do.call(
  rbind,
  imported[expected_imported_files$object_group == "tumor_vs_normal"]
)
rownames(brca_gsea_tumor_vs_normal_hallmark_go) <- NULL

brca_gsea_er_pos_vs_er_neg_hallmark_go <- do.call(
  rbind,
  imported[expected_imported_files$object_group == "er_pos_vs_er_neg"]
)
rownames(brca_gsea_er_pos_vs_er_neg_hallmark_go) <- NULL

brca_pa_merged <- rbind(
  brca_gsea_tumor_vs_normal_hallmark_go,
  brca_gsea_er_pos_vs_er_neg_hallmark_go
)
rownames(brca_pa_merged) <- NULL

brca_pa_results <- brca_pa_merged
brca_pathway_results_merged <- brca_pa_merged

write_tsv(
  brca_gsea_tumor_vs_normal_hallmark_go,
  file.path(intermediate_dir, "brca_gsea_tumor_vs_normal_hallmark_go.tsv")
)
write_tsv(
  brca_gsea_er_pos_vs_er_neg_hallmark_go,
  file.path(intermediate_dir, "brca_gsea_er_pos_vs_er_neg_hallmark_go.tsv")
)
write_tsv(
  brca_pa_merged,
  file.path(intermediate_dir, "brca_pa_merged.tsv")
)

## Clustering object -----------------------------------------------------------
## OmicsKit's pathway clustering functions require real gene-set memberships.
## Use memberships exported by the external GSEA pipeline when present; otherwise
## stop instead of fabricating a gene-set list from pathway names.
clustering_input <- brca_pa_merged[brca_pa_merged$FDR < 0.05, , drop = FALSE]
if (nrow(clustering_input) < 3L) {
  stop(
    "At least three significant pathways with FDR < 0.05 are required for ",
    "brca_pa_clustering. Observed: ", nrow(clustering_input),
    call. = FALSE
  )
}

brca_geneset_list <- make_geneset_list(clustering_input)

clustering_results <- data.frame(
  GeneSet = names(brca_geneset_list),
  FDR = clustering_input$FDR[
    match(
      names(brca_geneset_list),
      paste(
        clustering_input$COMPARISON,
        clustering_input$COLLECTION,
        clustering_input$NAME,
        sep = "::"
      )
    )
  ],
  stringsAsFactors = FALSE
)

brca_pa_similarity <- geneset_similarity(
  geneset_list = brca_geneset_list,
  results = clustering_results,
  fdr_th = 0.05
)

brca_pa_clustering <- do_clust(brca_pa_similarity)

## Keep the saved object compact. The network plotting functions need hclust and
## cluster_assignments; the heatmap can be regenerated from brca_pa_similarity.
if ("heatmap" %in% names(brca_pa_clustering)) {
  brca_pa_clustering$heatmap <- NULL
}

require_package("usethis", "to save BRCA pathway objects with usethis::use_data()")

usethis::use_data(
  brca_gsea_tumor_vs_normal_hallmark_go,
  brca_gsea_er_pos_vs_er_neg_hallmark_go,
  brca_pa_merged,
  brca_pa_results,
  brca_pathway_results_merged,
  brca_pa_clustering,
  brca_pa_similarity,
  brca_geneset_list,
  compress = "xz",
  overwrite = TRUE
)

message("Saved BRCA pathway example objects to data/.")
message("Wrote compact BRCA pathway TSV files to inst/extdata/brca/intermediate/.")
