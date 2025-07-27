#########################
# Function heatmap_GSEA #
#########################

#' Plot leading edge heatmaps from GSEA results.
#'
#' Generates heatmaps of leading edge genes for each gene set from GSEA output.
#'
#' @param main_dir Optional base directory. If supplied, it will be prepended to all relative file paths.
#' @param expression_file Path to the expression data file (tab-delimited) or relative to main_dir.
#' @param metadata_file Path to the metadata file (Excel) or relative to main_dir.
#' @param gmt_file Path to the GMT file defining gene sets or relative to main_dir.
#' @param ranked_genes_file Path to the ranked genes list file or relative to main_dir.
#' @param gsea_file Path to the GSEA results file with leading edge genes or relative to main_dir.
#' @param output_dir Directory to save heatmaps and optional TSV; default "leading_edge_heatmaps".
#' @param sample_col Name of the sample ID column in metadata; default "Sample".
#' @param group_col Name of the group column in metadata; default "group".
#' @param save_dataframe Logical; if TRUE, saves the merged data frame as TSV before plotting.
#' @return Saves one PDF and one JPG heatmap per gene set under output_dir; optionally saves intermediate TSV.
#' @export

heatmap_GSEA <- function(main_dir = NULL, expression_file, metadata_file, gmt_file,
                         ranked_genes_file, gsea_file, output_dir = "leading_edge_heatmaps",
                         sample_col = "Sample", group_col = "group", save_dataframe = FALSE)
{
  # Ensure required packages are installed
  if (!requireNamespace("readr", quietly = TRUE)) stop("Package \"readr\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("grDevices", quietly = TRUE)) stop("Package \"grDevices\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("tidyselect", quietly = TRUE)) stop("Package \"tidyselect\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Package \"openxlsx\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Package \"pheatmap\" must be installed to use this function.", call. = FALSE)

  # Prepend base directory if provided
  if (!is.null(main_dir)) {
    expression_file <- file.path(main_dir, expression_file)
    metadata_file <- file.path(main_dir, metadata_file)
    gmt_file <- file.path(main_dir, gmt_file)
    ranked_genes_file <- file.path(main_dir, ranked_genes_file)
    gsea_file <- file.path(main_dir, gsea_file)
    output_dir <- file.path(main_dir, output_dir)
  }

  # 1) Read and process GMT
  gmt_data <- readLines(gmt_file) %>%
    strsplit("\t") %>%
    lapply(function(x) data.frame(NAME = x[1], DESCRIPTION = x[2], GENES = paste(x[-c(1,2)], collapse = ","), stringsAsFactors = FALSE)) %>%
    dplyr::bind_rows()

  # 2) Read GSEA results and join genes
  gsea_df <- readr::read_tsv(gsea_file, show_col_types = FALSE) %>%
    dplyr::left_join(gmt_data %>% dplyr::select(NAME, GENES), by = "NAME")

  # 3) Read ranked genes list
  ranked_df <- readr::read_tsv(ranked_genes_file, show_col_types = FALSE)
  ranked_vector <- ranked_df[[1]]

  # 4) Internal helper: extract top-n genes from leading edge
  extract_top_n <- function(genes_str, n) {
    if (is.na(genes_str) || n <= 0) return(NA_character_)
    glist <- unlist(strsplit(genes_str, ","))
    glist <- glist[order(match(glist, ranked_vector), na.last = TRUE)]
    paste(utils::head(glist, n), collapse = ",")
  }

  # 5) Compute leading edge size and genes
  gsea_df <- gsea_df %>%
    dplyr::mutate(L.EDGE_size = ifelse(is.na(SIZE * tags), NA, ifelse((SIZE * tags) %% 1 <= 0.5, floor(SIZE * tags), ceiling(SIZE * tags)))) %>%
    dplyr::rowwise() %>% dplyr::mutate(LEADING_EDGE_GENES = extract_top_n(GENES, L.EDGE_size)) %>%
    dplyr::ungroup()

  # Save intermediate dataframe if requested
  if (save_dataframe) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    intermediate_file <- file.path(output_dir, "leading_edge_genes_df.tsv")
    readr::write_tsv(gsea_df, intermediate_file)
    message("Saved data frame to: ", intermediate_file)
  }

  # 6) Read metadata and prepare annotation
  meta <- openxlsx::read.xlsx(metadata_file) %>%
    dplyr::select(tidyselect::all_of(c(sample_col, group_col))) %>%
    dplyr::rename(Sample = tidyselect::all_of(sample_col), Group = tidyselect::all_of(group_col)) %>%
    as.data.frame()
  rownames(meta) <- meta$Sample

  # 7) Read expression data
  expr_raw <- utils::read.table(expression_file, header = TRUE, sep = "\t",
                                stringsAsFactors = FALSE, check.names = FALSE)
  # Determine gene-name column
  if ("NAME" %in% colnames(expr_raw)) {
    rownames(expr_raw) <- expr_raw$NAME
    expr_mat <- expr_raw[, setdiff(colnames(expr_raw), "NAME"), drop = FALSE]
  } else {
    gene_col <- colnames(expr_raw)[1]
    rownames(expr_raw) <- expr_raw[[gene_col]]
    expr_mat <- expr_raw[, -1, drop = FALSE]
  }
  # Clean sample names
  colnames(expr_mat) <- sub("^X", "", colnames(expr_mat))

  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 8) Loop through each gene set and plot heatmap
  for (i in seq_len(nrow(gsea_df))) {
    geneset_name <- gsea_df$NAME[i]
    leading_genes <- unlist(strsplit(gsea_df$LEADING_EDGE_GENES[i], ","))
    genes_present <- leading_genes[leading_genes %in% rownames(expr_mat)]
    if (length(genes_present) == 0) next

    heatmap_mat <- expr_mat[genes_present, , drop = FALSE]
    common_samps <- intersect(colnames(heatmap_mat), rownames(meta))
    if (length(common_samps) == 0) next

    heatmap_mat <- heatmap_mat[, common_samps, drop = FALSE]
    annot_col <- data.frame(Group = meta[common_samps, "Group"])
    rownames(annot_col) <- common_samps

    # Dynamic sizing
    w <- 10
    h <- max(5, nrow(heatmap_mat) * 0.1 + 2)

    # PDF output
    grDevices::pdf(file.path(output_dir, paste0(geneset_name, "_heatmap.pdf")), width = w, height = h)
    pheatmap::pheatmap(
      heatmap_mat,
      main = geneset_name,
      color = grDevices::colorRampPalette(c("blue","white","red"))(30),
      scale = "row",
      clustering_distance_rows = "euclidean",
      cluster_cols = FALSE,
      clustering_method = "complete",
      fontsize_row = 6,
      fontsize_col = 7,
      annotation_col = annot_col,
      border_color = NA,
      cellheight = 5,
      cellwidth = 8
    )
    grDevices::dev.off()

    # JPG output
    grDevices::jpeg(file.path(output_dir, paste0(geneset_name, "_heatmap.jpg")),
                    width = w * 100, height = h * 100, res = 150)
    pheatmap::pheatmap(
      heatmap_mat,
      main = geneset_name,
      color = grDevices::colorRampPalette(c("blue","white","red"))(30),
      scale = "row",
      clustering_distance_rows = "euclidean",
      cluster_cols = FALSE,
      clustering_method  = "complete",
      fontsize_row = 6,
      fontsize_col = 7,
      annotation_col = annot_col,
      border_color = NA,
      cellheight = 5,
      cellwidth = 8
    )
    grDevices::dev.off()
  }

  message("Heatmaps saved in: ", normalizePath(output_dir))
  return(TRUE)
}
