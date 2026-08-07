## 09_render_brca_figures.R
##
## Purpose:
##   Render real prerendered TCGA-BRCA vignette figures.
##
## Inputs:
##   - Package-ready BRCA objects in data/*.rda.
##   - Real omics-layer track files in inst/extdata/brca/omics_layers/.
##
## Outputs:
##   - DEA, concordance, omics-layer, and modeling figures in
##     vignettes/figures/.
##   - Pathway-analysis and pathway-clustering figures in
##     vignettes/figures_PA/.
##
## Expected file locations:
##   - data/
##   - inst/extdata/brca/omics_layers/
##   - vignettes/figures/
##   - vignettes/figures_PA/
##
## Notes:
##   - This script renders figures only.
##   - It does not download inputs, rerun DEA, edit vignettes, or create
##     synthetic biological results.

project_file <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    here::here(...)
  } else {
    file.path(getwd(), ...)
  }
}

require_file <- function(path, label = path) {
  if (!file.exists(path)) {
    stop(
      "Required input is missing for BRCA figure rendering: ", label, "\n",
      "Expected location: ", path,
      call. = FALSE
    )
  }
  invisible(path)
}

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package \"", package, "\" is required to render BRCA vignette figures.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

source_function_file <- function(path) {
  require_file(path, paste0("R function file ", path))
  source(path, local = .GlobalEnv)
}

has_object <- function(name) {
  exists(name, envir = .GlobalEnv, inherits = FALSE)
}

get_object <- function(name) {
  get(name, envir = .GlobalEnv, inherits = FALSE)
}

require_object <- function(name, label = name) {
  if (!has_object(name)) {
    stop(
      "Required object is missing for BRCA figure rendering: ", label, "\n",
      "Expected object name: ", name, "\n",
      "Create it from the appropriate data-raw/brca preparation script first.",
      call. = FALSE
    )
  }
  get_object(name)
}

get_first_object <- function(candidates, label) {
  found <- candidates[vapply(candidates, has_object, logical(1))]
  if (length(found) == 0L) {
    stop(
      "Required object is missing for BRCA figure rendering: ", label, "\n",
      "Expected one of: ", paste(candidates, collapse = ", "), "\n",
      "Only real BRCA objects are accepted; simulated package examples are not used.",
      call. = FALSE
    )
  }
  get_object(found[[1L]])
}

assert_columns <- function(data, cols, label) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0L) {
    stop(
      label, " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

resolve_path <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("Track path is missing or empty.", call. = FALSE)
  }
  if (file.exists(path)) {
    return(path)
  }
  project_file(path)
}

load_all_package_data <- function(data_dir) {
  require_file(data_dir, "data/ directory")
  data_files <- list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)
  if (length(data_files) == 0L) {
    stop("No .rda files found in data/.", call. = FALSE)
  }

  for (file in data_files) {
    load(file, envir = .GlobalEnv)
  }

  invisible(data_files)
}

save_gg_png <- function(plot, path, width = 7, height = 5) {
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  invisible(path)
}

save_gg_pdf <- function(plot, path, width = 7, height = 5) {
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
  invisible(path)
}

open_png_device <- function(path, width = 7, height = 5) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = width, height = height, units = "in", res = 300)
  } else {
    grDevices::png(path, width = width, height = height, units = "in", res = 300)
  }
  invisible(TRUE)
}

save_device_png <- function(path, width = 7, height = 5, expr) {
  open_png_device(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  invisible(path)
}

save_device_pdf <- function(path, width = 7, height = 5, expr) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
  invisible(path)
}

prepare_dea_table <- function(data, label) {
  data <- as.data.frame(data)

  if (!"logFC" %in% names(data)) {
    stop(label, " must contain a logFC column.", call. = FALSE)
  }

  if (!"padj" %in% names(data)) {
    if ("adj.P.Val" %in% names(data)) {
      data$padj <- data$adj.P.Val
    } else {
      stop(label, " must contain padj or adj.P.Val.", call. = FALSE)
    }
  }

  if (!"gene" %in% names(data)) {
    gene_col <- intersect(c("gene_symbol", "gene_id", "feature_id", "id"), names(data))[1L]
    if (is.na(gene_col)) {
      stop(label, " must contain a gene, gene_symbol, gene_id, feature_id, or id column.", call. = FALSE)
    }
    data$gene <- data[[gene_col]]
  }

  data
}

matrix_with_required_gene <- function(candidates, gene, label) {
  for (name in candidates) {
    if (!has_object(name)) {
      next
    }
    mat <- as.matrix(get_object(name))
    if (gene %in% rownames(mat)) {
      return(mat)
    }
  }

  stop(
    label, " requires gene/feature `", gene, "` in one of: ",
    paste(candidates, collapse = ", "), ".",
    call. = FALSE
  )
}

make_er_annotations <- function(metadata, matrix_colnames, label) {
  metadata <- as.data.frame(metadata)
  assert_columns(metadata, "ER_group", label)

  clean_sample16 <- function(x) {
    x <- gsub("\\.", "-", as.character(x))
    x[is.na(x) | !nzchar(x)] <- NA_character_
    substring(x, 1L, 16L)
  }

  id_candidates <- list(rownames = rownames(metadata))
  for (col in intersect(
    c(
      "rna_sample16", "rppa_sample16", "sample16", "sampleID", "sample",
      "RNA_genomic_id", "RPPA_genomic_id"
    ),
    names(metadata)
  )) {
    id_candidates[[col]] <- clean_sample16(metadata[[col]])
  }

  match_counts <- vapply(
    id_candidates,
    function(ids) sum(as.character(ids) %in% matrix_colnames, na.rm = TRUE),
    numeric(1)
  )

  best <- names(match_counts)[which.max(match_counts)]
  if (length(best) == 0L || match_counts[[best]] == 0L) {
    stop(label, " could not be matched to matrix column names.", call. = FALSE)
  }

  metadata$id <- as.character(id_candidates[[best]])
  metadata <- metadata[metadata$id %in% matrix_colnames, , drop = FALSE]
  metadata <- metadata[match(matrix_colnames, metadata$id), , drop = FALSE]

  if (any(is.na(metadata$id))) {
    stop(label, " could not be matched to all matrix columns.", call. = FALSE)
  }

  metadata$sample_type <- factor(
    metadata$ER_group,
    levels = c("ER_negative", "ER_positive")
  )

  metadata
}

positive_for_vsb <- function(mat, gene, label) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  values <- mat[gene, ]
  if (any(!is.finite(values))) {
    stop(label, " contains non-finite values for ", gene, ".", call. = FALSE)
  }

  min_value <- min(mat, na.rm = TRUE)
  if (is.finite(min_value) && min_value <= 0) {
    warning(
      label,
      " contains values <= 0; applying a display-only positive shift before ",
      "calling nice_VSB(), which log2-transforms its input.",
      call. = FALSE
    )
    mat <- mat - min_value + 1e-3
  }

  mat
}

palette_for_levels <- function(x) {
  n <- length(levels(factor(x)))
  grDevices::hcl.colors(max(n, 2L), palette = "Dark 3")
}

standardize_pa_table <- function(data, label) {
  data <- as.data.frame(data)

  rename_if_present <- function(target, aliases) {
    if (target %in% names(data)) {
      return(invisible(TRUE))
    }
    hit <- aliases[aliases %in% names(data)][1L]
    if (!is.na(hit)) {
      names(data)[names(data) == hit] <<- target
    }
    invisible(TRUE)
  }

  rename_if_present("NAME", c("GeneSet", "geneset", "pathway", "Pathway"))
  rename_if_present("COLLECTION", c("collection", "MSigDB", "msigdb"))
  rename_if_present("NES", c("nes", "normalized_enrichment_score"))
  rename_if_present("FDR", c("fdr", "padj", "adj.P.Val", "FDR q-val", "qval"))
  rename_if_present("COMPARISON", c("comparison", "Comparison", "contrast"))

  assert_columns(data, c("NAME", "COLLECTION", "NES", "FDR"), label)
  data$FDR <- suppressWarnings(as.numeric(data$FDR))
  data$NES <- suppressWarnings(as.numeric(data$NES))
  data <- data[is.finite(data$FDR) & is.finite(data$NES), , drop = FALSE]

  if (nrow(data) == 0L) {
    stop(label, " has no finite NES/FDR rows.", call. = FALSE)
  }

  data
}

select_pa_comparison <- function(data, patterns, label) {
  data <- standardize_pa_table(data, label)

  if (!"COMPARISON" %in% names(data)) {
    return(data)
  }

  comparisons <- as.character(data$COMPARISON)
  keep <- Reduce(
    `|`,
    lapply(patterns, function(pattern) grepl(pattern, comparisons, ignore.case = TRUE))
  )

  if (!any(keep)) {
    stop(
      label,
      " could not find a matching comparison in COMPARISON. Available values: ",
      paste(unique(comparisons), collapse = ", "),
      call. = FALSE
    )
  }

  data[keep, , drop = FALSE]
}

top_pa_rows <- function(data, n = 20L) {
  data <- data[order(data$FDR, -abs(data$NES)), , drop = FALSE]
  data[seq_len(min(n, nrow(data))), , drop = FALSE]
}

make_brca1_annotation <- function() {
  data.frame(
    geneID = "ENSG00000012048",
    symbol = "BRCA1",
    chromosome = "17",
    gene_start = 43044295,
    gene_end = 43125483,
    strand = "-",
    stringsAsFactors = FALSE
  )
}

make_circos_tracks <- function(track_files) {
  list(
    `RNA-seq` = list(
      path = resolve_path(track_files$rna_1Mb_bed %||% track_files$rna_1Mb_bw),
      type = "histogram",
      color = "#457B9D",
      height = 0.10
    ),
    CNV = list(
      path = resolve_path(track_files$cnv_bed %||% track_files$cnv_bw),
      type = "line",
      color = "#E63946",
      height = 0.10
    ),
    `Methylation 450k` = list(
      path = resolve_path(track_files$methyl_1Mb_bed %||% track_files$methyl_1Mb_bw),
      type = "points",
      color = "#2A9D8F",
      height = 0.10
    ),
    Mutations = list(
      path = resolve_path(track_files$mutation_bed),
      type = "mutation",
      color = "#E9C46A",
      height = 0.08
    )
  )
}

`%||%` <- function(x, y) {
  if (!is.null(x) && length(x) > 0L && !is.na(x) && nzchar(x)) {
    x
  } else {
    y
  }
}

## Setup -----------------------------------------------------------------------
required_packages <- c(
  "ggplot2", "magrittr", "dplyr", "tibble", "matrixStats", "scales",
  "ggrepel", "ggpubr", "umap", "tsne", "survminer", "Gviz",
  "GenomicRanges", "GenomeInfoDb", "rtracklayer", "circlize",
  "ComplexHeatmap", "cowplot", "patchwork", "cluster", "igraph",
  "ggraph", "tidygraph"
)
invisible(lapply(required_packages, require_package))

suppressPackageStartupMessages({
  library(ggplot2)
  library(magrittr)
})

function_files <- c(
  "nice_PCA.R",
  "nice_UMAP.R",
  "nice_tSNE.R",
  "nice_Volcano.R",
  "nice_VSB.R",
  "nice_ConcordanceScatter.R",
  "crossLayerCorr.R",
  "nice_GenomeTrack.R",
  "nice_circos.R",
  "nice_KM.R",
  "nice_forest.R",
  "nice_ROC.R",
  "plot_PA.R",
  "doclust_PA.R",
  "plotclust_PA.R"
)

invisible(lapply(project_file("R", function_files), source_function_file))

figures_dir <- project_file("vignettes", "figures")
figures_pa_dir <- project_file("vignettes", "figures_PA")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_pa_dir, recursive = TRUE, showWarnings = FALSE)

load_all_package_data(project_file("data"))

## DEA figures -----------------------------------------------------------------
brca_rna_dea_er_pos_vs_er_neg <- require_object("brca_rna_dea_er_pos_vs_er_neg")
brca_rna_dea_tumor_vs_normal <- require_object("brca_rna_dea_tumor_vs_normal")
brca_rppa_dea_gene_er_pos_vs_er_neg <- get_first_object(
  c("brca_rppa_dea_gene_er_pos_vs_er_neg", "brca_rppa_dea_feature_er_pos_vs_er_neg"),
  "RPPA ER_positive vs ER_negative DEA results"
)
brca_rna_vst_or_logexpr_small <- require_object("brca_rna_vst_or_logexpr_small")
brca_rna_metadata_er_shared <- require_object("brca_rna_metadata_er_shared")

rna_small <- as.matrix(brca_rna_vst_or_logexpr_small)
rna_meta_er <- make_er_annotations(
  brca_rna_metadata_er_shared,
  colnames(rna_small),
  "brca_rna_metadata_er_shared"
)

er_colors <- c("ER_negative" = "#457B9D", "ER_positive" = "#E76F51")

save_gg_png(
  nice_PCA(
    object = rna_small,
    annotations = rna_meta_er,
    variables = c(fill = "sample_type"),
    legend_names = c(fill = "ER status"),
    colors = er_colors,
    size = 3.5,
    title = "TCGA-BRCA RNA-seq PCA by ER status"
  ),
  file.path(figures_dir, "brca_pca_rna_er_status.png"),
  width = 6.5,
  height = 5.5
)

save_gg_png(
  nice_UMAP(
    object = rna_small,
    annotations = rna_meta_er,
    neighbors = 10,
    epochs = 1000,
    seed = 2025,
    variables = c(fill = "sample_type"),
    legend_names = c(fill = "ER status"),
    colors = er_colors,
    size = 3.5,
    title = "TCGA-BRCA RNA-seq UMAP by ER status",
    transform = FALSE
  ),
  file.path(figures_dir, "brca_umap_rna_er_status.png"),
  width = 6.5,
  height = 5.5
)

save_gg_png(
  nice_tSNE(
    object = rna_small,
    annotations = rna_meta_er,
    perplexity = max(3, min(30, floor((ncol(rna_small) - 1L) / 3L))),
    max_iterations = 2000,
    seed = 2025,
    variables = c(fill = "sample_type"),
    legend_names = c(fill = "ER status"),
    colors = er_colors,
    size = 3.5,
    title = "TCGA-BRCA RNA-seq t-SNE by ER status",
    transform = FALSE
  ),
  file.path(figures_dir, "brca_tsne_rna_er_status.png"),
  width = 6.5,
  height = 5.5
)

rna_de_er <- prepare_dea_table(
  brca_rna_dea_er_pos_vs_er_neg,
  "brca_rna_dea_er_pos_vs_er_neg"
)
rna_de_tn <- prepare_dea_table(
  brca_rna_dea_tumor_vs_normal,
  "brca_rna_dea_tumor_vs_normal"
)
rppa_de_er <- prepare_dea_table(
  brca_rppa_dea_gene_er_pos_vs_er_neg,
  "brca_rppa_dea_gene_er_pos_vs_er_neg"
)

save_gg_png(
  nice_Volcano(
    results = rna_de_er,
    x_var = "logFC",
    y_var = "padj",
    label_var = "gene",
    title = "RNA-seq: ER-positive vs ER-negative",
    cutoff_y = 0.05,
    cutoff_x = 1,
    genes = intersect(c("ESR1", "PGR", "GATA3", "ERBB2", "MKI67"), rna_de_er$gene)
  ),
  file.path(figures_dir, "brca_volcano_rna_er_pos_vs_er_neg.png"),
  width = 8,
  height = 6
)

save_gg_png(
  nice_Volcano(
    results = rna_de_tn,
    x_var = "logFC",
    y_var = "padj",
    label_var = "gene",
    title = "RNA-seq: Tumor vs Normal",
    cutoff_y = 0.05,
    cutoff_x = 1,
    genes = intersect(c("ESR1", "ERBB2", "BRCA1", "MKI67", "EPCAM"), rna_de_tn$gene)
  ),
  file.path(figures_dir, "brca_volcano_rna_tumor_vs_normal.png"),
  width = 8,
  height = 6
)

save_gg_png(
  nice_Volcano(
    results = rppa_de_er,
    x_var = "logFC",
    y_var = "padj",
    label_var = "gene",
    title = "RPPA: ER-positive vs ER-negative",
    cutoff_y = 0.05,
    cutoff_x = 0.20,
    genes = intersect(c("ERBB2", "ESR1", "PGR", "GATA3"), rppa_de_er$gene)
  ),
  file.path(figures_dir, "brca_volcano_rppa_er_pos_vs_er_neg.png"),
  width = 8,
  height = 6
)

rna_vsb_mat <- matrix_with_required_gene(
  c("brca_rna_vst_or_logexpr_small", "brca_rna_expr_er_shared_filtered"),
  "ESR1",
  "RNA ESR1 VSB figure"
)
rna_vsb_meta <- make_er_annotations(
  brca_rna_metadata_er_shared,
  colnames(rna_vsb_mat),
  "brca_rna_metadata_er_shared"
)

save_gg_png(
  nice_VSB(
    object = positive_for_vsb(rna_vsb_mat, "ESR1", "RNA VSB matrix"),
    annotations = rna_vsb_meta,
    variables = c(fill = "sample_type"),
    genename = "ESR1",
    symbol = "RNA-seq",
    labels = c("ER-", "ER+"),
    categories = c("ER_negative", "ER_positive"),
    colors = er_colors,
    shapes = 21,
    markersize = 2.5
  ),
  file.path(figures_dir, "brca_vsb_rna_ESR1.png"),
  width = 5.5,
  height = 5
)

brca_rppa_expr_gene_er_shared <- require_object("brca_rppa_expr_gene_er_shared")
rppa_vsb_mat <- matrix_with_required_gene(
  c("brca_rppa_expr_gene_er_shared"),
  "ERBB2",
  "RPPA ERBB2 VSB figure"
)
rppa_vsb_meta <- make_er_annotations(
  brca_rna_metadata_er_shared,
  colnames(rppa_vsb_mat),
  "brca_rna_metadata_er_shared for RPPA VSB"
)

save_gg_png(
  nice_VSB(
    object = positive_for_vsb(rppa_vsb_mat, "ERBB2", "RPPA VSB matrix"),
    annotations = rppa_vsb_meta,
    variables = c(fill = "sample_type"),
    genename = "ERBB2",
    symbol = "RPPA",
    labels = c("ER-", "ER+"),
    categories = c("ER_negative", "ER_positive"),
    colors = er_colors,
    shapes = 21,
    markersize = 2.5
  ),
  file.path(figures_dir, "brca_vsb_rppa_ERBB2.png"),
  width = 5.5,
  height = 5
)

## Concordance figures ---------------------------------------------------------
brca_concordance_rna_rppa_er <- require_object("brca_concordance_rna_rppa_er")
brca_crosslayercorr_rna_rppa_er <- require_object("brca_crosslayercorr_rna_rppa_er")
brca_rna_expr_er_shared_filtered <- require_object("brca_rna_expr_er_shared_filtered")

save_gg_png(
  nice_ConcordanceScatter(
    brca_concordance_rna_rppa_er,
    x_label = "RNA-seq log2FC",
    y_label = "RPPA log2FC",
    genes_label = intersect(c("ESR1", "ERBB2", "PGR", "GATA3"), rppa_de_er$gene),
    method = "spearman",
    logfc_col = "logFC"
  ),
  file.path(figures_dir, "brca_concordance_scatter_rna_rppa_er.png"),
  width = 6.5,
  height = 5.5
)

shared_rna_rppa_genes <- intersect(
  rownames(brca_rna_expr_er_shared_filtered),
  rownames(brca_rppa_expr_gene_er_shared)
)

if (length(shared_rna_rppa_genes) < 2L) {
  stop(
    "At least two shared RNA/RPPA genes are required to render brca_crosslayercorr_rna_rppa_er.",
    call. = FALSE
  )
}

crosslayer_for_plot <- crossLayerCorr(
  mat_x = brca_rna_expr_er_shared_filtered,
  mat_y = brca_rppa_expr_gene_er_shared,
  method = "spearman",
  top_n = min(length(shared_rna_rppa_genes), 146L),
  plot = TRUE
)

save_gg_png(
  crosslayer_for_plot$plot,
  file.path(figures_dir, "brca_crosslayercorr_rna_rppa_er.png"),
  width = 8,
  height = 4.8
)

## Omics-layer figures ---------------------------------------------------------
brca_omics_layer_tracks <- require_object("brca_omics_layer_tracks")
track_files <- brca_omics_layer_tracks$track_files

required_track_keys <- c(
  "rna_bw", "rna_1Mb_bw", "rna_1Mb_bed", "cnv_bw", "cnv_bed", "methyl_brca1_bw",
  "methyl_1Mb_bw", "methyl_1Mb_bed", "mutation_bed"
)
missing_track_keys <- setdiff(required_track_keys, names(track_files))
if (length(missing_track_keys) > 0L) {
  stop(
    "brca_omics_layer_tracks$track_files is missing required key(s): ",
    paste(missing_track_keys, collapse = ", "),
    call. = FALSE
  )
}

resolved_track_files <- lapply(track_files, resolve_path)
invisible(lapply(names(resolved_track_files), function(key) {
  require_file(
    resolved_track_files[[key]],
    paste0("brca_omics_layer_tracks$track_files$", key)
  )
}))

brca1_region <- brca_omics_layer_tracks$brca1_region_hg38
if (is.null(brca1_region)) {
  brca1_region <- c(chr = "chr17", start = 42944295, end = 43610338)
}

genome_track_args <- list(
  region = brca1_region,
  genome_label = "hg38",
  organism = "hsapiens_gene_ensembl",
  ensembl_version = "current",
  annotations = make_brca1_annotation(),
  tracks = list(
    `RNA-seq` = resolved_track_files$rna_bw,
    CNV = resolved_track_files$cnv_bw,
    `Methylation 450k` = resolved_track_files$methyl_brca1_bw,
    Mutations = resolved_track_files$mutation_bed
  ),
  highlight_genes = "BRCA1",
  track_sizes = c(1, 3, 2, 2, 2, 1.5)
)

save_device_png(
  file.path(figures_dir, "brca_genome_track_brca1_TCGA-A1-A0SH-01.png"),
  width = 13,
  height = 9,
  do.call(nice_GenomeTrack, genome_track_args)
)

do.call(
  nice_GenomeTrack,
  c(
    genome_track_args,
    list(export_pdf = file.path(figures_dir, "brca_genome_track_brca1_TCGA-A1-A0SH-01.pdf"))
  )
)

circos_tracks <- make_circos_tracks(resolved_track_files)
invisible(lapply(circos_tracks, function(track) {
  require_file(track$path, paste0("circos track ", track$path))
}))

circos_args <- list(
  genome_build = "hg38",
  data_tracks = circos_tracks,
  ideogram = FALSE,
  chromosome_index = paste0("chr", c(1:22, "X", "Y")),
  show_legend = TRUE,
  show_labels = TRUE,
  track_height = 0.09
)

save_device_png(
  file.path(figures_dir, "brca_circos_TCGA-A1-A0SH-01.png"),
  width = 8,
  height = 8,
  do.call(nice_circos, circos_args)
)

do.call(
  nice_circos,
  c(
    circos_args,
    list(export_pdf = file.path(figures_dir, "brca_circos_TCGA-A1-A0SH-01.pdf"))
  )
)

## Modeling figures ------------------------------------------------------------
brca_clinical_modeling_data <- require_object("brca_clinical_modeling_data")
brca_cox_adjusted_clinical <- require_object("brca_cox_adjusted_clinical")
brca_cox_univariable_clinical <- require_object("brca_cox_univariable_clinical")
brca_roc_stage_clinical <- require_object("brca_roc_stage_clinical")

assert_columns(
  brca_clinical_modeling_data,
  c("age_group", "OS.time", "OS"),
  "brca_clinical_modeling_data"
)

km_data <- brca_clinical_modeling_data[
  !is.na(brca_clinical_modeling_data$age_group) &
    is.finite(brca_clinical_modeling_data$OS.time) &
    !is.na(brca_clinical_modeling_data$OS),
  ,
  drop = FALSE
]

if (length(unique(km_data$age_group)) < 2L) {
  stop("age_group must have at least two observed levels for brca_km_age_OS.", call. = FALSE)
}

save_gg_png(
  nice_KM(
    data = km_data,
    gene = "age_group",
    time_var = "OS.time",
    event_var = "OS",
    title_prefix = "",
    colors = palette_for_levels(km_data$age_group),
    conf_int = TRUE,
    risk_table = FALSE,
    legend_pos = c(0.75, 0.85)
  ),
  file.path(figures_dir, "brca_km_age_OS.png"),
  width = 6.8,
  height = 5.4
)

forest_input <- brca_cox_adjusted_clinical
if (!is.data.frame(forest_input) ||
    nrow(forest_input) == 0L ||
    "skip_reason" %in% names(forest_input)) {
  forest_input <- brca_cox_univariable_clinical
}

forest_plot <- nice_forest(
  forest_input,
  title = "TCGA-BRCA clinical Cox model for PFI",
  sort_by = "p.value",
  base_size = 11
)

save_gg_png(
  forest_plot,
  file.path(figures_dir, "brca_forest_cox_clinical_pfi.png"),
  width = 8,
  height = 6
)

save_gg_pdf(
  forest_plot,
  file.path(figures_dir, "brca_forest_cox_clinical_pfi.pdf"),
  width = 8,
  height = 6
)

roc_plot <- if (is.list(brca_roc_stage_clinical) &&
                !is.null(brca_roc_stage_clinical$plot)) {
  brca_roc_stage_clinical$plot
} else if (inherits(brca_roc_stage_clinical, "ggplot")) {
  brca_roc_stage_clinical
} else {
  stop(
    "brca_roc_stage_clinical must be a ggplot object or a list with a $plot element.",
    call. = FALSE
  )
}

save_gg_png(
  roc_plot,
  file.path(figures_dir, "brca_roc_stage_clinical.png"),
  width = 6.5,
  height = 5.5
)

## Pathway figures -------------------------------------------------------------
brca_pa_results <- get_first_object(
  c(
    "brca_pa_results",
    "brca_pathway_results",
    "brca_gsea_results",
    "brca_pathway_results_merged"
  ),
  "BRCA pathway analysis result table"
)

pa_single_tn <- if (has_object("brca_pa_tumor_vs_normal")) {
  get_object("brca_pa_tumor_vs_normal")
} else {
  select_pa_comparison(
    brca_pa_results,
    c("tumor.*normal", "tumour.*normal", "tumor_vs_normal", "tumorvsnormal"),
    "BRCA pathway Tumor vs Normal results"
  )
}

pa_single_tn <- top_pa_rows(standardize_pa_table(pa_single_tn, "BRCA Tumor vs Normal pathway table"), 18L)

save_gg_png(
  splot_PA(
    data = pa_single_tn,
    geneset_col = "NAME",
    collection_col = "COLLECTION",
    nes_col = "NES",
    fdr_col = "FDR",
    fill_limits = c(0, min(5, max(-log10(pa_single_tn$FDR), na.rm = TRUE))),
    theme_params = list(
      side_label_size = 9,
      geneset_text_size = 2.6,
      collection_text_size = 2.4,
      axis_title_size = 11,
      axis_text_size_x = 9,
      panel_widths = c(2, 18, 12, 3, 5, 3)
    )
  ),
  file.path(figures_pa_dir, "brca_pa_single_tumor_vs_normal.png"),
  width = 8,
  height = 6
)

pa_multi_er <- if (has_object("brca_pa_er_status")) {
  get_object("brca_pa_er_status")
} else {
  select_pa_comparison(
    brca_pa_results,
    c("er.*positive", "er_positive", "erpositive", "er.*negative", "er_negative", "ernegative"),
    "BRCA pathway ER status results"
  )
}

pa_multi_er <- standardize_pa_table(pa_multi_er, "BRCA ER status pathway table")
top_multi_names <- unique(top_pa_rows(pa_multi_er, 12L)$NAME)
pa_multi_er_plot <- pa_multi_er[pa_multi_er$NAME %in% top_multi_names, , drop = FALSE]

save_gg_png(
  multiplot_PA(
    data = pa_multi_er_plot,
    comparison_col = "COMPARISON",
    facet_col = "NAME",
    axis_y = "NES",
    fdr_col = "FDR",
    ncol_wrap = 3,
    fill_limits = c(0, min(5, max(-log10(pa_multi_er_plot$FDR), na.rm = TRUE))),
    theme_params = list(
      axis_title_size = 12,
      axis_text_size_x = 8,
      axis_text_size_y = 8,
      strip_text_size = 8,
      hline_size = 0.4
    )
  ),
  file.path(figures_pa_dir, "brca_pa_multi_er_status.png"),
  width = 9,
  height = 6.5
)

brca_geneset_list <- get_first_object(
  c("brca_geneset_list", "brca_pa_geneset_list", "brca_msigdb_geneset_list"),
  "BRCA pathway gene-set list"
)

pa_for_clustering <- standardize_pa_table(pa_multi_er, "BRCA pathway table for clustering")
pa_for_similarity <- data.frame(
  GeneSet = pa_for_clustering$NAME,
  FDR = pa_for_clustering$FDR,
  stringsAsFactors = FALSE
)

brca_pa_similarity <- if (has_object("brca_pa_similarity")) {
  get_object("brca_pa_similarity")
} else if (has_object("brca_pa_jaccard")) {
  get_object("brca_pa_jaccard")
} else {
  geneset_similarity(
    geneset_list = brca_geneset_list,
    results = pa_for_similarity,
    fdr_th = 0.05
  )
}

brca_pa_clustering <- if (has_object("brca_pa_clustering")) {
  get_object("brca_pa_clustering")
} else if (has_object("brca_pa_clustering_result")) {
  get_object("brca_pa_clustering_result")
} else {
  do_clust(brca_pa_similarity)
}

clustering_plot_list <- network_clust_gg(
  brca_pa_similarity,
  clust_result = brca_pa_clustering,
  jaccard_threshold = 0.25,
  min_degree = 1,
  type = "combined",
  seed = 2025
)

save_gg_png(
  clustering_plot_list$combined,
  file.path(figures_pa_dir, "brca_pa_clustering.png"),
  width = 8,
  height = 7
)

message("Rendered all requested TCGA-BRCA vignette figures.")
