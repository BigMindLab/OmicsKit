####################
# Function splot_PA #
####################

utils::globalVariables(c(
  "NAME", "GENES", "SIZE", "tags", "L.EDGE_size",
  "numeric_cols", "LEADING EDGE", "signal",
  "FDR q-val", "Log10FDR", "FWER p-val", "Comparison"
))

#' Pathway analysis visualization for a single comparison
#'
#' Generates a publication-quality multi-panel pathway enrichment plot for a
#' single comparison using patchwork. Gene sets appear on the y-axis grouped
#' by MSigDB collection, NES on the x-axis, and -log10(FDR) as fill color.
#' Six panels are assembled side by side: a "Pathways" label, gene set names,
#' the NES bar chart, collection labels, a "MSigDB" label, and the color legend.
#'
#' For visualizing enrichment across multiple comparisons, use
#' [multiplot_PA()] instead.
#'
#' @param data A data frame of pathway analysis results for a single
#'   comparison. Typically the output of [merge_PA()] filtered to one value
#'   of the `COMPARISON` column, or results from a single CAMERA/GSEA run.
#'   Must contain the columns specified by `geneset_col`, `collection_col`,
#'   `nes_col`, and `fdr_col`.
#' @param geneset_col Name of the column containing gene set labels shown on
#'   the y-axis. Default: `"NAME"`.
#' @param collection_col Name of the column containing MSigDB collection
#'   labels used to group gene sets (e.g., `"KEGG"`, `"HALLMARK"`, `"GO"`).
#'   Default: `"COLLECTION"`.
#' @param nes_col Name of the column containing NES values (x-axis).
#'   Default: `"NES"`.
#' @param fdr_col Name of the column containing FDR values. `-log10(FDR)` is
#'   computed internally and used as the fill color. Default: `"FDR"`.
#' @param order One of `"desc"` or `"asc"`. Sort order for NES values on the
#'   y-axis. Default: `"desc"`.
#' @param fill_limits Numeric vector of length 2 setting the color scale range
#'   for `-log10(FDR)`. Values outside this range are clamped to the nearest
#'   limit. For example, `fill_limits = c(0, 5)` maps all gene sets with
#'   `-log10(FDR) >= 5` (i.e., FDR <= 0.00001) to the maximum color (red),
#'   and any value below 0 to the minimum color (white). Useful when a few
#'   gene sets have extreme significance that washes out color variation in the
#'   rest. Default: `NULL` (auto  uses the actual data range).
#' @param fill_palette Character vector of two colors for the fill gradient
#'   (low to high -log10(FDR)). Default: `c("white", "red")`.
#' @param theme_params Named list to override default theme parameters.
#'   See Details.
#'
#' @details
#'   `theme_params` accepts any of the following named elements:
#'   \describe{
#'     \item{`side_label_size`}{Size for "Pathways" and "MSigDB" labels.
#'       Default: `35`.}
#'     \item{`geneset_text_size`}{Text size for gene set labels. Default: `5`.}
#'     \item{`collection_text_size`}{Text size for collection labels.
#'       Default: `5`.}
#'     \item{`panel_widths`}{Patchwork relative widths for the 6 panels.
#'       Default: `c(4, 25, 15, 3, 10, 3)`.}
#'     \item{`col_size`}{Border linewidth for `geom_col`. Default: `1`.}
#'     \item{`axis_title_size`}{Font size for axis titles. Default: `45`.}
#'     \item{`axis_text_size_x`}{Font size for x-axis labels. Default: `30`.}
#'     \item{`tick_size`}{Linewidth for axis ticks. Default: `1.5`.}
#'     \item{`tick_length`}{Length of axis ticks in cm. Default: `0.3`.}
#'     \item{`panel_spacing_single`}{Spacing between facets. Default: `4`.}
#'   }
#'
#' @return A `patchwork` object combining six ggplot2 panels.
#'
#' @examples
#' \dontrun{
#' gsea_results <- merge_PA("path/to/gsea_results/")
#'
#' # Filter to one comparison
#' single <- gsea_results[gsea_results$COMPARISON == "TumorVsNormal", ]
#'
#' splot_PA(
#'   data           = single,
#'   geneset_col    = "NAME",
#'   collection_col = "COLLECTION",
#'   nes_col        = "NES",
#'   fdr_col        = "FDR"
#' )
#'
#' # Cap color scale at -log10(FDR) = 5 so subtle differences are visible
#' # (gene sets with FDR <= 0.00001 all get the same max red color)
#' splot_PA(
#'   data        = single,
#'   geneset_col = "NAME", collection_col = "COLLECTION",
#'   nes_col     = "NES",  fdr_col        = "FDR",
#'   fill_limits = c(0, 5)
#' )
#' }
#'
#' @seealso [multiplot_PA()] for multi-comparison faceted barplots;
#'   [merge_PA()] to generate the input data frame;
#'   [camera_results] for a minimal example dataset.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom patchwork plot_layout
#' @importFrom utils modifyList
#' @export

splot_PA <- function(data,
                     geneset_col    = "NAME",
                     collection_col = "COLLECTION",
                     nes_col        = "NES",
                     fdr_col        = "FDR",
                     order          = "desc",
                     fill_limits    = NULL,
                     fill_palette   = c("white", "red"),
                     theme_params   = list()) {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package \"patchwork\" must be installed to use this function.", call. = FALSE)
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    stop("Package \"cowplot\" must be installed to use this function.", call. = FALSE)
  }

  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)

  for (col in c(geneset_col, collection_col, nes_col, fdr_col)) {
    if (!col %in% colnames(data)) {
      stop("Column '", col, "' not found in `data`.", call. = FALSE)
    }
  }

  order <- match.arg(order, c("desc", "asc"))

  defaults <- list(
    side_label_size      = 35,
    geneset_text_size    = 5,
    collection_text_size = 5,
    panel_widths         = c(4, 25, 15, 3, 10, 3),
    col_size             = 1,
    axis_title_size      = 45,
    axis_text_size_x     = 30,
    tick_size            = 1.5,
    tick_length          = 0.3,
    panel_spacing_single = 4
  )
  params <- utils::modifyList(defaults, theme_params)

  # Always compute -log10(FDR) internally
  data$tmp_log10FDR <- -log10(data[[fdr_col]])

  data <- data[order(data[[nes_col]], decreasing = (order == "desc")), ]

  df <- data[, c(geneset_col, collection_col, nes_col, "tmp_log10FDR")]
  colnames(df) <- c("Geneset", "Collection", "NES", "tmp_log10FDR")
  df$Geneset    <- factor(df$Geneset,    levels = rev(unique(df$Geneset)))
  df$Collection <- factor(df$Collection, levels = unique(df$Collection))

  plot_text_pathways <- ggplot() +
    annotate("text", label = "Pathways", fontface = "bold.italic", angle = 90,
             size = params$side_label_size, x = 0, y = 0.5) +
    theme_void()

  plot_left <- ggplot(df, aes(y = .data$Geneset, x = 0)) +
    geom_text(aes(label = .data$Geneset), hjust = 1,
              size = params$geneset_text_size) +
    theme_void() +
    theme(axis.text.y = element_blank(), plot.margin = margin(0, 0, 0, -50))

  plot_center <- ggplot(df, aes(x = .data$NES, y = .data$Geneset,
                                fill = .data$tmp_log10FDR)) +
    geom_col(color = "black", linewidth = params$col_size) +
    scale_fill_gradient(low = fill_palette[1], high = fill_palette[2],
                        limits = fill_limits,
                        breaks = scales::pretty_breaks()) +
    scale_y_discrete(position = "right") +
    facet_grid(Collection ~ ., scales = "free_y", space = "free_y") +
    theme_bw() + labs(x = "NES", y = "") +
    theme(
      axis.text.y       = element_blank(),
      strip.background  = element_rect(fill = "white", color = "black", linewidth = 1),
      axis.ticks.y      = element_line(linewidth = params$tick_size),
      axis.ticks.length = grid::unit(params$tick_length, "cm"),
      strip.text.y      = element_text(size = 1),
      legend.position   = "none",
      axis.title.x      = element_text(size = params$axis_title_size),
      axis.text.x       = element_text(size = params$axis_text_size_x),
      panel.spacing     = grid::unit(params$panel_spacing_single, "lines")
    )

  plot_text_msigdb <- ggplot() +
    annotate("text", label = "MSigDB", fontface = "bold.italic", angle = 90,
             size = params$side_label_size, x = 0, y = 0.5) +
    theme_void()

  plot_right <- ggplot(df, aes(y = .data$Geneset, x = 1.5)) +
    geom_text(
      aes(label = ifelse(duplicated(.data$Collection), "",
                         as.character(.data$Collection))),
      hjust = 0.5, size = params$collection_text_size, fontface = "bold"
    ) +
    facet_grid(Collection ~ ., scales = "free_y", space = "free", switch = "y") +
    theme_void() +
    theme(strip.text.y  = element_text(size = params$collection_text_size),
          panel.spacing = grid::unit(1, "lines"))

  plot_legend <- ggplot(df, aes(x = .data$NES, y = .data$Geneset,
                                fill = .data$tmp_log10FDR)) +
    geom_tile() +
    scale_fill_gradient(
      low   = fill_palette[1], high = fill_palette[2],
      name  = expression(-log[10] ~ FDR), limits = fill_limits,
      guide = guide_colorbar(ticks.colour = "black", ticks.linewidth = 1.5,
                             draw.ulim = TRUE, draw.llim = TRUE)
    ) +
    theme_bw() +
    theme(
      legend.title      = element_text(size = 44, face = "bold"),
      legend.text       = element_text(size = 30),
      legend.key.size   = grid::unit(1.5, "cm"),
      legend.key.height = grid::unit(2, "cm"),
      legend.spacing    = grid::unit(3.5, "cm"),
      legend.box.margin = margin(10, 20, 10, 10)
    )

  plot_right_legend <- cowplot::get_legend(plot_legend)

  final_plot <- plot_text_pathways + plot_left + plot_center + plot_right +
    plot_text_msigdb + plot_right_legend +
    patchwork::plot_layout(ncol = 6, widths = params$panel_widths)

  return(final_plot)
}


########################
# Function multiplot_PA #
########################

#' Pathway analysis visualization across multiple comparisons
#'
#' Generates a faceted barplot showing NES values across multiple comparisons
#' for a set of gene sets. Each facet represents one gene set and bars
#' represent the NES per comparison, colored by -log10(FDR). This layout makes
#' it easy to compare how enrichment of gene sets changes across conditions
#' (e.g., TumorVsNormal, MetastasisVsNormal).
#'
#' All comparisons must be combined in a single data frame with a column
#' identifying each comparison  as produced by [merge_PA()].
#'
#' For visualizing a single comparison with full collection grouping, use
#' [splot_PA()] instead.
#'
#' @param data A data frame of pathway analysis results containing two or more
#'   comparisons. Typically the output of [merge_PA()].
#' @param comparison_col Name of the column identifying each comparison.
#'   Appears on the x-axis of each facet. Default: `"COMPARISON"`.
#' @param facet_col Name of the column used to define facets  one facet per
#'   unique value. Can be the original gene set name column (e.g., `"NAME"`)
#'   or a manually curated column with cleaner or shorter labels
#'   (e.g., `"clean_name"`). Default: `"NAME"`.
#' @param axis_y Name of the column to use for the y-axis. Default: `"NES"`.
#' @param fdr_col Name of the column containing FDR values. `-log10(FDR)` is
#'   computed internally and used as the fill color. Default: `"FDR"`.
#' @param comparison_order Character vector specifying the left-to-right order
#'   of comparisons on the x-axis of each facet. For example,
#'   `comparison_order = c("BvsA", "CvsA")` places `BvsA` on the left and
#'   `CvsA` on the right. If `NULL` (default), the order follows the factor
#'   levels of `comparison_col` as they appear in `data`.
#' @param custom_labels Named character vector of x-axis tick labels. Useful
#'   for shortening comparison names on the axis. For example,
#'   `custom_labels = c(TumorVsNormal = "Tumor", MetastasisVsNormal = "Mets")`.
#'   Default: `NULL`.
#' @param ncol_wrap Integer. Number of columns in `facet_wrap`. Default: `2`.
#' @param free_y Logical. If `TRUE`, each facet uses its own y-axis scale.
#'   Default: `TRUE`.
#' @param fill_limits Numeric vector of length 2 setting the color scale range
#'   for `-log10(FDR)`. Values outside this range are clamped to the nearest
#'   limit. For example, `fill_limits = c(0, 5)` maps all gene sets with
#'   `-log10(FDR) >= 5` (FDR <= 0.00001) to maximum red, and any value below
#'   0 to white. Useful when one gene set has extreme significance that makes
#'   the rest appear uniform. Default: `NULL` (auto).
#' @param fill_palette Character vector of two colors for the fill gradient
#'   (low to high -log10(FDR)). Default: `c("white", "red")`.
#' @param theme_params Named list to override default theme parameters.
#'   See Details.
#'
#' @details
#'   `theme_params` accepts any of the following named elements:
#'   \describe{
#'     \item{`bar_col`}{Bar border color. Default: `"black"`.}
#'     \item{`bar_size`}{Bar border linewidth. Default: `0.5`.}
#'     \item{`bar_width`}{Bar width. Default: `0.6`.}
#'     \item{`hline_size`}{Linewidth for horizontal line at y = 0. Default: `2`.}
#'     \item{`axis_title_size`}{Font size for axis titles. Default: `45`.}
#'     \item{`axis_text_size_x`}{Font size for x-axis labels. Default: `30`.}
#'     \item{`axis_text_size_y`}{Font size for y-axis labels. Default: `50`.}
#'     \item{`tick_size`}{Linewidth for axis ticks. Default: `1.5`.}
#'     \item{`tick_length`}{Length of axis ticks in cm. Default: `0.3`.}
#'     \item{`strip_text_size`}{Font size for facet strip labels. Default: `50`.}
#'     \item{`panel_spacing_multi`}{Spacing between facets. Default: `0.6`.}
#'   }
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' gsea_results <- merge_PA("path/to/gsea_results/")
#'
#' # Basic multi-comparison plot
#' multiplot_PA(
#'   data           = gsea_results,
#'   comparison_col = "COMPARISON",
#'   facet_col      = "NAME",
#'   fdr_col        = "FDR",
#'   ncol_wrap      = 3
#' )
#'
#' # Control left-to-right order of comparisons on the x-axis
#' multiplot_PA(
#'   data             = gsea_results,
#'   comparison_col   = "COMPARISON",
#'   facet_col        = "NAME",
#'   fdr_col          = "FDR",
#'   comparison_order = c("BvsA", "CvsA")   # BvsA on the left, CvsA on the right
#' )
#'
#' # Use cleaner facet labels and shorten x-axis tick names
#' gsea_results$clean_name <- gsub("_", " ", gsea_results$NAME)
#'
#' multiplot_PA(
#'   data             = gsea_results,
#'   comparison_col   = "COMPARISON",
#'   facet_col        = "clean_name",
#'   fdr_col          = "FDR",
#'   comparison_order = c("BvsA", "CvsA"),
#'   custom_labels    = c(BvsA = "Tumor", CvsA = "Metastasis")
#' )
#' }
#'
#' @seealso [splot_PA()] for single-comparison patchwork plots;
#'   [merge_PA()] to generate the input data frame;
#'   [camera_results] for a minimal example dataset.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom utils modifyList
#' @export

multiplot_PA <- function(data,
                         comparison_col   = "COMPARISON",
                         facet_col        = "NAME",
                         axis_y           = "NES",
                         fdr_col          = "FDR",
                         comparison_order = NULL,
                         custom_labels    = NULL,
                         ncol_wrap        = 2,
                         free_y           = TRUE,
                         fill_limits      = NULL,
                         fill_palette     = c("white", "red"),
                         theme_params     = list()) {

  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)

  for (col in c(comparison_col, axis_y, fdr_col)) {
    if (!col %in% colnames(data)) {
      stop("Column '", col, "' not found in `data`.", call. = FALSE)
    }
  }

  if (!facet_col %in% colnames(data)) {
    stop(
      "Column '", facet_col, "' not found in `data`. ",
      "`facet_col` can be the original gene set name column (e.g., 'NAME') ",
      "or a manually curated column with cleaner labels.",
      call. = FALSE
    )
  }

  defaults <- list(
    bar_col             = "black",
    bar_size            = 0.5,
    bar_width           = 0.6,
    hline_size          = 2,
    axis_title_size     = 45,
    axis_text_size_x    = 30,
    axis_text_size_y    = 50,
    tick_size           = 1.5,
    tick_length         = 0.3,
    strip_text_size     = 50,
    panel_spacing_multi = 0.6
  )
  params <- utils::modifyList(defaults, theme_params)

  # Always compute -log10(FDR) internally
  data$tmp_log10FDR <- -log10(data[[fdr_col]])

  # Apply comparison order if specified
  if (!is.null(comparison_order)) {
    missing_comps <- setdiff(comparison_order, unique(data[[comparison_col]]))
    if (length(missing_comps) > 0) {
      warning(
        "The following values in `comparison_order` were not found in `",
        comparison_col, "`: ",
        paste(missing_comps, collapse = ", "),
        call. = FALSE
      )
    }
    data[[comparison_col]] <- factor(data[[comparison_col]],
                                     levels = comparison_order)
  }

  p <- ggplot(data, aes(x    = .data[[comparison_col]],
                        y    = .data[[axis_y]],
                        fill = .data$tmp_log10FDR)) +
    geom_bar(stat = "identity",
             color     = params$bar_col,
             linewidth = params$bar_size,
             width     = params$bar_width) +
    scale_fill_gradient(
      low    = fill_palette[1], high = fill_palette[2],
      limits = fill_limits, oob = scales::squish,
      name   = expression(-log[10] ~ FDR),
      guide  = guide_colorbar(barwidth = 3, barheight = 18)
    ) +
    labs(x = "Comparisons", y = axis_y) +
    theme_bw() +
    theme(
      axis.line.x       = element_blank(),
      axis.line         = element_line(linewidth = 0.5),
      axis.title.x      = element_text(size = params$axis_title_size),
      axis.title.y      = element_text(size = params$axis_title_size),
      axis.text.x       = element_text(size = params$axis_text_size_x),
      axis.text.y       = element_text(size = params$axis_text_size_y),
      axis.ticks        = element_line(linewidth = params$tick_size),
      axis.ticks.length = grid::unit(params$tick_length, "cm"),
      strip.text        = element_text(size = params$strip_text_size),
      panel.spacing     = grid::unit(params$panel_spacing_multi, "lines")
    ) +
    geom_hline(yintercept = 0, linewidth = params$hline_size) +
    expand_limits(y = 0) +
    facet_wrap(~ .data[[facet_col]], ncol = ncol_wrap,
               scales = if (free_y) "free_y" else "fixed")

  if (!is.null(custom_labels)) {
    p <- p + scale_x_discrete(labels = custom_labels)
  }

  return(p)
}


########################
# Function heatmap_PA  #
########################

#' Plot leading edge heatmaps from pathway analysis results
#'
#' Generates heatmaps of gene expression for each gene set in `pa_data_annot`,
#' using the `all_genes`, `le_genes` (GSEA output only), and/or `top_genes`
#' columns produced by [addgenesPA()]. Genes within each heatmap are ordered
#' by their position in `ranked_genes`.
#'
#' The recommended workflow before calling this function is:
#' ```r
#' gsl          <- list_gmts("path/to/gmt/")
#' pa_data      <- merge_PA("path/to/pa_data/")
#' ranked       <- deseq2_results$gene_id[order(deseq2_results$stat,
#'                                              decreasing = TRUE)]
#' gene_lists   <- getgenesPA(pa_data, gsl, ranked, genes = c("all", "le"))
#' pa_annot     <- addgenesPA(pa_data, gene_lists)
#'
#' heatmap_PA(
#'   expression_data = vst_counts,
#'   metadata        = sampledata,
#'   pa_data_annot   = pa_annot,
#'   ranked_genes    = ranked,
#'   plot_genes      = c("all_genes", "le_genes")
#' )
#' ```
#'
#' @param expression_data A numeric matrix or data frame of expression values
#'   with gene symbols or Ensembl IDs as row names and sample IDs as column
#'   names. Recommended input: VST-transformed counts from [vst_counts] or
#'   normalized coutns [norm_counts].
#' @param metadata A data frame of sample annotations. Must contain a column
#'   matching `sample_col` (sample IDs) and a column matching `group_col`
#'   (condition labels, e.g., `"Control"`, `"Treatment"`).
#' @param pa_data_annot  A data frame of pathway analysis results enriched with
#'   gene columns. Must contain the column `NAME` and at least one of
#'   `all_genes`, `le_genes`, or `top_genes` (comma-separated gene symbols per
#'   gene set). Typically the output of [addgenesPA()].
#' @param ranked_genes A character vector of gene symbols ordered by their
#'   ranking metric (e.g., stat, log2FC or signal-to-noise ratio), used to sort
#'   genes within each heatmap row.
#' @param plot_genes Character vector specifying which gene columns to plot.
#'   One or both of `"all_genes"` and `"le_genes"`, and `"top_genes"`. Each
#'   selection produces its own set of output files in a dedicated subfolder.
#'   Default: `c("all_genes", "le_genes")`.
#' @param sample_col Name of the sample ID column in `metadata`.
#'   Default: `"Sample"`.
#' @param group_col Name of the condition/group column in `metadata`
#'   (e.g., `"Control"` vs `"Treatment"`). Used for heatmap column
#'   annotations. Default: `"group"`.
#' @param out_dir Character. Path to the output directory. Subdirectories are
#'   created automatically based on `pdf`, `jpg`, and `plot_genes`:
#'   * `<out_dir>/pdf/all_genes/`
#'   * `<out_dir>/pdf/le_genes/`
#'   * `<out_dir>/pdf/top_genes/`
#'   * `<out_dir>/jpg/all_genes/`
#'   * `<out_dir>/jpg/le_genes/`
#'   * `<out_dir>/jpg/top_genes/`
#'   Default: `"heatmaps_PA"`.
#' @param pdf Logical. If `TRUE`, saves PDF heatmaps. Default: `TRUE`.
#' @param jpg Logical. If `TRUE`, saves JPG heatmaps. Default: `TRUE`.
#'
#' @return Invisibly returns `TRUE` upon completion. Saves heatmap files to
#'   the corresponding subdirectories under `out_dir`.
#'
#' @examples
#' \dontrun{
#' data(vst_counts)
#' data(sampledata)
#' data(deseq2_results)
#' data(gsea_results)
#' data(geneset_list)
#'
#' ranked    <- deseq2_results$gene_id[order(deseq2_results$stat,
#'                                           decreasing = TRUE)]
#'
#' # ── Example 1: GSEA results (all_genes + le_genes) ────
#' pa_single  <- gsea_results[gsea_results$COMPARISON == "TumorVsNormal", ]
#' gene_lists <- getgenesPA(pa_single, geneset_list, ranked,
#'                          genes = c("all", "le"))
#' pa_annot   <- addgenesPA(pa_single, gene_lists)
#'
#' heatmap_PA(
#'   expression_data = vst_counts,
#'   metadata        = sampledata,
#'   pa_data_annot   = pa_annot,
#'   ranked_genes    = ranked,
#'   plot_genes      = c("all_genes", "le_genes"),
#'   sample_col      = "patient_id",
#'   group_col       = "sample_type",
#'   out_dir         = "heatmaps_gsea",
#'   pdf             = TRUE,
#'   jpg             = TRUE
#' )
#' # Creates:
#' #   heatmaps_gsea/pdf/all_genes/<geneset>_heatmap.pdf
#' #   heatmaps_gsea/pdf/le_genes/<geneset>_heatmap.pdf
#' #   heatmaps_gsea/jpg/all_genes/<geneset>_heatmap.jpg
#' #   heatmaps_gsea/jpg/le_genes/<geneset>_heatmap.jpg
#'
#' # ── Example 2: CAMERA results (all_genes + top_genes)
#' # camera_results does not contain leading edge information.
#' # Use genes = "top" with a manually set top fraction instead.
#' # Note: top_genes are rank-based and do NOT represent true leading edge genes.
#' data(camera_results)
#' camera_pa      <- camera_results
#' colnames(camera_pa)[colnames(camera_pa) == "GeneSet"] <- "NAME"
#' camera_pa$SIZE <- sapply(camera_pa$NAME,
#'                          function(x) length(geneset_list[[x]]))
#' camera_pa$top  <- 0.25
#'
#' gene_lists_cam <- getgenesPA(camera_pa, geneset_list, ranked,
#'                              genes = c("all", "top"))
#' pa_annot_cam   <- addgenesPA(camera_pa, gene_lists_cam)
#'
#' heatmap_PA(
#'   expression_data = vst_counts,
#'   metadata        = sampledata,
#'   pa_data_annot   = pa_annot_cam,
#'   ranked_genes    = ranked,
#'   plot_genes      = c("all_genes", "top_genes"),
#'   sample_col      = "patient_id",
#'   group_col       = "sample_type",
#'   out_dir         = "heatmaps_camera"
#' )
#' }
#'
#' @seealso [getgenesPA()] for gene extraction;
#'   [addgenesPA()] to generate `pa_data_annot`;
#'   [list_gmts()] to generate the geneset list;
#'   [merge_PA()] to generate `pa_data`;
#'   [vst_counts] for an example expression matrix.
#'
#' @export

heatmap_PA <- function(expression_data,
                       metadata,
                       pa_data_annot,
                       ranked_genes,
                       plot_genes  = c("all_genes", "le_genes"),
                       sample_col  = "Sample",
                       group_col   = "group",
                       out_dir     = "heatmaps_PA",
                       pdf         = TRUE,
                       jpg         = TRUE) {

  if (!requireNamespace("pheatmap",  quietly = TRUE)) {
    stop("Package \"pheatmap\" must be installed to use this function.",  call. = FALSE)
  }
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("Package \"grDevices\" must be installed to use this function.", call. = FALSE)
  }

  # --- Input validation ----
  if (!is.matrix(expression_data) && !is.data.frame(expression_data)) {
    stop("`expression_data` must be a matrix or data frame.", call. = FALSE)
  }
  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data frame.", call. = FALSE)
  }
  if (!sample_col %in% colnames(metadata)) {
    stop("`sample_col` ('", sample_col, "') not found in `metadata`.", call. = FALSE)
  }
  if (!group_col %in% colnames(metadata)) {
    stop("`group_col` ('", group_col, "') not found in `metadata`.", call. = FALSE)
  }
  if (!is.data.frame(pa_data_annot)) {
    stop("`pa_data_annot` must be a data frame.", call. = FALSE)
  }
  if (!"NAME" %in% colnames(pa_data_annot)) {
    stop("`pa_data_annot` must contain a column named 'NAME'.", call. = FALSE)
  }

  plot_genes <- match.arg(plot_genes,
                          choices      = c("all_genes", "le_genes", "top_genes"),
                          several.ok   = TRUE)

  missing_cols <- setdiff(plot_genes, colnames(pa_data_annot))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) not found in `pa_data_annot`: ",
      paste(missing_cols, collapse = ", "),
      ". Run addgenesPA() first.",
      call. = FALSE
    )
  }

  if (!is.character(ranked_genes) || length(ranked_genes) == 0) {
    stop("`ranked_genes` must be a non-empty character vector.", call. = FALSE)
  }

  # --- Prepare metadata ----
  meta <- metadata[, c(sample_col, group_col), drop = FALSE]
  colnames(meta) <- c("Sample", "Group")
  rownames(meta) <- meta$Sample

  expr_mat     <- as.matrix(expression_data)
  common_samps <- intersect(colnames(expr_mat), rownames(meta))

  if (length(common_samps) == 0) {
    stop(
      "No common samples between `expression_data` columns and `metadata` rows. ",
      "Check that `sample_col` matches the column names of `expression_data`.",
      call. = FALSE
    )
  }

  # --- Internal heatmap drawer --
  .draw_heatmap <- function(mat, title) {
    annot_col        <- data.frame(Group = meta[colnames(mat), "Group"])
    rownames(annot_col) <- colnames(mat)
    pheatmap::pheatmap(
      mat,
      main                     = title,
      color                    = grDevices::colorRampPalette(
        c("blue", "white", "red"))(30),
      scale                    = "row",
      clustering_distance_rows = "euclidean",
      cluster_cols             = FALSE,
      clustering_method        = "complete",
      fontsize_row             = 6,
      fontsize_col             = 7,
      annotation_col           = annot_col,
      border_color             = NA,
      cellheight               = 5,
      cellwidth                = 8
    )
  }

  # --- Loop over requested gene column types ---
  n_plotted <- 0L

  for (gene_col in plot_genes) {

    if (pdf) {
      pdf_sub <- file.path(out_dir, "pdf", gene_col)
      if (!dir.exists(pdf_sub)) dir.create(pdf_sub, recursive = TRUE)
    }
    if (jpg) {
      jpg_sub <- file.path(out_dir, "jpg", gene_col)
      if (!dir.exists(jpg_sub)) dir.create(jpg_sub, recursive = TRUE)
    }

    for (i in seq_len(nrow(pa_data_annot))) {

      gs        <- pa_data_annot$NAME[i]
      genes_str <- pa_data_annot[[gene_col]][i]

      if (is.na(genes_str) || nchar(genes_str) == 0) next

      genes_vec     <- strsplit(genes_str, ",")[[1]]
      genes_in_mat  <- genes_vec[genes_vec %in% rownames(expr_mat)]

      if (length(genes_in_mat) == 0) next

      gene_ranks    <- match(genes_in_mat, ranked_genes)
      genes_ordered <- genes_in_mat[order(gene_ranks, na.last = TRUE)]
      heatmap_mat   <- expr_mat[genes_ordered, common_samps, drop = FALSE]

      w <- 10
      h <- max(5, nrow(heatmap_mat) * 0.1 + 2)

      if (pdf) {
        grDevices::pdf(file.path(pdf_sub, paste0(gs, "_heatmap.pdf")),
                       width = w, height = h)
        .draw_heatmap(heatmap_mat, gs)
        grDevices::dev.off()
      }

      if (jpg) {
        grDevices::jpeg(file.path(jpg_sub, paste0(gs, "_heatmap.jpg")),
                        width = w * 100, height = h * 100, res = 150)
        .draw_heatmap(heatmap_mat, gs)
        grDevices::dev.off()
      }

      n_plotted <- n_plotted + 1L
    }
  }

  message("Done. ", n_plotted, " heatmap(s) saved in: ",
          normalizePath(out_dir))
  invisible(TRUE)
}

############################
# Function heatmap_path_PA #
############################

utils::globalVariables(c(
  "NAME", "GENES", "SIZE", "tags", "L.EDGE_size"
))

#' Plot leading edge heatmaps from GSEA analysis results using file paths
#'
#' Generates one heatmap per gene set from GSEA/CAMERA/PADOG output by reading
#' all required inputs from file paths. For each gene set, the leading edge
#' genes are extracted, ordered by their rank in the ranked gene list, and
#' plotted as a scaled row heatmap against the expression matrix.
#'
#' This function is the file-path-based alternative to [heatmap_PA()], which
#' accepts R objects directly. Use this version when working from raw output
#' files on disk (e.g., directly after running `GSEA_merge.sh`).
#'
#' @param main_dir Character or `NULL`. Optional base directory prepended to
#'   all relative file paths. If `NULL` (default), all paths are used as-is.
#' @param expression_file Character. Path to a tab-delimited expression data
#'   file. Rows are genes (first column or a column named `NAME` used as row
#'   names), columns are sample IDs. Recommended input: VST-transformed counts.
#' @param metadata_file Character. Path to an Excel (`.xlsx`) metadata file.
#'   Must contain a column matching `sample_col` (sample IDs) and a column
#'   matching `group_col` (condition labels, e.g., `"Control"`,
#'   `"Treatment"`).
#' @param gmt_file Character. Path to a `.gmt` file defining gene sets. Each
#'   row contains: gene set name (column 1), description (column 2, ignored),
#'   and gene symbols (columns 3+).
#' @param ranked_genes_file Character. Path to a tab-delimited file where the
#'   first column contains gene symbols ordered by their ranking metric (e.g.,
#'   log2FC or signal-to-noise ratio), from most positive to most negative.
#'   Used to order leading edge genes within each heatmap.
#' @param gsea_file Character. Path to a GSEA results `.tsv` file containing
#'   at least the columns `NAME`, `SIZE`, and `tags` (from the `LEADING EDGE`
#'   column parsed by [merge_PA()]).
#' @param output_dir Character. Directory where heatmap files are saved.
#'   Created automatically if it does not exist. Default:
#'   `"leading_edge_heatmaps"`.
#' @param sample_col Name of the sample ID column in the metadata file.
#'   Default: `"Sample"`.
#' @param group_col Name of the condition/group column in the metadata file
#'   (e.g., `"Control"` vs `"Treatment"`). Used for heatmap column
#'   annotations. Default: `"group"`.
#' @param save_dataframe Logical. If `TRUE`, saves the intermediate data frame
#'   (gene sets with computed leading edge genes) as a `.tsv` file in
#'   `output_dir` before plotting. Useful for inspection or reuse.
#'   Default: `FALSE`.
#'
#' @return Invisibly returns `TRUE` upon completion. Saves two files per gene
#'   set in `output_dir`:
#'   * `<geneset_name>_heatmap.pdf`
#'   * `<geneset_name>_heatmap.jpg`
#'
#'   If `save_dataframe = TRUE`, also saves
#'   `<output_dir>/leading_edge_genes_df.tsv`.
#'
#' @note For a more flexible workflow that accepts R objects directly (avoiding
#'   repeated file reads), use [heatmap_PA()] instead, which takes
#'   `expression_data`, `metadata`, and `pa_data_annot` as R objects and
#'   integrates with [getgenesPA()] and [addgenesPA()].
#'
#' @examples
#' \dontrun{
#' # Run with all files in a single base directory
#' heatmap_path_PA(
#'   main_dir          = "path/to/analysis/",
#'   expression_file   = "vst_expression.tsv",
#'   metadata_file     = "metadata.xlsx",
#'   gmt_file          = "genesets.gmt",
#'   ranked_genes_file = "ranked_genes.tsv",
#'   gsea_file         = "gsea_results.tsv",
#'   output_dir        = "leading_edge_heatmaps",
#'   sample_col        = "Sample",
#'   group_col         = "group",
#'   save_dataframe    = TRUE
#' )
#' # Saves:
#' #   leading_edge_heatmaps/<geneset>_heatmap.pdf
#' #   leading_edge_heatmaps/<geneset>_heatmap.jpg
#' #   leading_edge_heatmaps/leading_edge_genes_df.tsv  (if save_dataframe = TRUE)
#'
#' # Run with absolute paths (no main_dir)
#' heatmap_path_PA(
#'   expression_file   = "/data/vst_counts.tsv",
#'   metadata_file     = "/data/metadata.xlsx",
#'   gmt_file          = "/data/h.all.v2023.gmt",
#'   ranked_genes_file = "/data/ranked_genes.tsv",
#'   gsea_file         = "/data/gsea_results.tsv"
#' )
#' }
#'
#' @seealso [heatmap_PA()] for the R-object-based alternative;
#'   [getgenesPA()] and [addgenesPA()] for extracting leading edge genes
#'   from R objects; [merge_PA()] to generate the GSEA results input;
#'   [list_gmts()] to load GMT files as R objects.
#'
#' @importFrom magrittr %>%
#' @export

heatmap_path_PA <- function(main_dir = NULL,
                            expression_file,
                            metadata_file,
                            gmt_file,
                            ranked_genes_file,
                            gsea_file,
                            output_dir     = "leading_edge_heatmaps",
                            sample_col     = "Sample",
                            group_col      = "group",
                            save_dataframe = FALSE) {

  if (!requireNamespace("readr",      quietly = TRUE)) stop("Package \"readr\" must be installed to use this function.",      call. = FALSE)
  if (!requireNamespace("grDevices",  quietly = TRUE)) stop("Package \"grDevices\" must be installed to use this function.",  call. = FALSE)
  if (!requireNamespace("tidyselect", quietly = TRUE)) stop("Package \"tidyselect\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("openxlsx",   quietly = TRUE)) stop("Package \"openxlsx\" must be installed to use this function.",   call. = FALSE)
  if (!requireNamespace("pheatmap",   quietly = TRUE)) stop("Package \"pheatmap\" must be installed to use this function.",   call. = FALSE)

  # Prepend base directory if provided
  if (!is.null(main_dir)) {
    expression_file   <- file.path(main_dir, expression_file)
    metadata_file     <- file.path(main_dir, metadata_file)
    gmt_file          <- file.path(main_dir, gmt_file)
    ranked_genes_file <- file.path(main_dir, ranked_genes_file)
    gsea_file         <- file.path(main_dir, gsea_file)
    output_dir        <- file.path(main_dir, output_dir)
  }

  # 1) Read and process GMT
  gmt_data <- readLines(gmt_file) %>%
    strsplit("\t") %>%
    lapply(function(x) data.frame(
      NAME  = x[1],
      GENES = paste(x[-c(1, 2)], collapse = ","),
      stringsAsFactors = FALSE
    )) %>%
    dplyr::bind_rows()

  # 2) Read GSEA results and join genes
  gsea_df <- readr::read_tsv(gsea_file, show_col_types = FALSE) %>%
    dplyr::left_join(
      dplyr::select(gmt_data, NAME, GENES),
      by = "NAME"
    )

  # 3) Read ranked genes list
  ranked_df     <- readr::read_tsv(ranked_genes_file, show_col_types = FALSE)
  ranked_vector <- ranked_df[[1]]

  # 4) Internal helper: extract top-n genes ordered by rank
  extract_top_n <- function(genes_str, n) {
    if (is.na(genes_str) || n <= 0) return(NA_character_)
    glist <- unlist(strsplit(genes_str, ","))
    glist <- glist[order(match(glist, ranked_vector), na.last = TRUE)]
    paste(utils::head(glist, n), collapse = ",")
  }

  # 5) Compute leading edge size and extract genes
  gsea_df <- gsea_df %>%
    dplyr::mutate(
      L.EDGE_size = ifelse(
        is.na(SIZE * tags), NA,
        ifelse((SIZE * tags) %% 1 <= 0.5,
               floor(SIZE * tags),
               ceiling(SIZE * tags))
      )
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(LEADING_EDGE_GENES = extract_top_n(GENES, L.EDGE_size)) %>%
    dplyr::ungroup()

  # Optionally save intermediate data frame
  if (save_dataframe) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    intermediate_file <- file.path(output_dir, "leading_edge_genes_df.tsv")
    readr::write_tsv(gsea_df, intermediate_file)
    message("Saved intermediate data frame to: ", intermediate_file)
  }

  # 6) Read metadata and prepare annotation
  meta <- openxlsx::read.xlsx(metadata_file) %>%
    dplyr::select(tidyselect::all_of(c(sample_col, group_col))) %>%
    dplyr::rename(
      Sample = tidyselect::all_of(sample_col),
      Group  = tidyselect::all_of(group_col)
    ) %>%
    as.data.frame()
  rownames(meta) <- meta$Sample

  # 7) Read expression matrix
  expr_raw <- utils::read.table(expression_file, header = TRUE, sep = "\t",
                                stringsAsFactors = FALSE, check.names = FALSE)

  if ("NAME" %in% colnames(expr_raw)) {
    rownames(expr_raw) <- expr_raw$NAME
    expr_mat <- expr_raw[, setdiff(colnames(expr_raw), "NAME"), drop = FALSE]
  } else {
    gene_col <- colnames(expr_raw)[1]
    rownames(expr_raw) <- expr_raw[[gene_col]]
    expr_mat <- expr_raw[, -1, drop = FALSE]
  }

  # Clean sample names (remove leading "X" added by R)
  colnames(expr_mat) <- sub("^X", "", colnames(expr_mat))

  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 8) Loop through each gene set and plot heatmap
  for (i in seq_len(nrow(gsea_df))) {

    geneset_name  <- gsea_df$NAME[i]
    leading_genes <- unlist(strsplit(gsea_df$LEADING_EDGE_GENES[i], ","))
    genes_present <- leading_genes[leading_genes %in% rownames(expr_mat)]

    if (length(genes_present) == 0) next

    common_samps <- intersect(colnames(expr_mat), rownames(meta))
    if (length(common_samps) == 0) next

    heatmap_mat <- expr_mat[genes_present, common_samps, drop = FALSE]
    annot_col   <- data.frame(Group = meta[common_samps, "Group"])
    rownames(annot_col) <- common_samps

    w <- 10
    h <- max(5, nrow(heatmap_mat) * 0.1 + 2)

    .draw <- function() {
      pheatmap::pheatmap(
        heatmap_mat,
        main                     = geneset_name,
        color                    = grDevices::colorRampPalette(c("blue", "white", "red"))(30),
        scale                    = "row",
        clustering_distance_rows = "euclidean",
        cluster_cols             = FALSE,
        clustering_method        = "complete",
        fontsize_row             = 6,
        fontsize_col             = 7,
        annotation_col           = annot_col,
        border_color             = NA,
        cellheight               = 5,
        cellwidth                = 8
      )
    }

    grDevices::pdf(file.path(output_dir, paste0(geneset_name, "_heatmap.pdf")),
                   width = w, height = h)
    .draw()
    grDevices::dev.off()

    grDevices::jpeg(file.path(output_dir, paste0(geneset_name, "_heatmap.jpg")),
                    width = w * 100, height = h * 100, res = 150)
    .draw()
    grDevices::dev.off()
  }

  message("Heatmaps saved in: ", normalizePath(output_dir))
  invisible(TRUE)
}

