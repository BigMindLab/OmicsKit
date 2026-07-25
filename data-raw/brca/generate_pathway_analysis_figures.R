## data-raw/brca/generate_pathway_analysis_figures.R
## Pre-render figures for pathway_analysis.Rmd

suppressPackageStartupMessages({
  library(OmicsKit)
  library(ggplot2)
})

# suppressPackageStartupMessages({
#   library(ggplot2)
# })
#
# devtools::load_all(".", quiet = TRUE)

out_dir <- file.path("vignettes", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## 1. Load data
## -----------------------------

data("brca_pa_merged")
data("brca_geneset_list")
data("brca_ranked_genes")
data("brca_pa_similarity")
data("brca_pa_clustering")

## -----------------------------
## 2. Helpers
## -----------------------------

save_plot_safely <- function(plot_object,
                             filename,
                             width = 8,
                             height = 6,
                             dpi = 300) {
  out_file <- file.path(out_dir, filename)

  if (inherits(plot_object, "ggplot")) {
    ggplot2::ggsave(
      filename = out_file,
      plot     = plot_object,
      width    = width,
      height   = height,
      dpi      = dpi
    )
  } else {
    pdf(out_file, width = width, height = height, units = "in", res = dpi)
    print(plot_object)
    dev.off()
  }

  message("Saved: ", normalizePath(out_file))
  invisible(out_file)
}

extract_plot <- function(x,
                         preferred_names = c(
                           "combined",
                           "superterms",
                           "labels",
                           "clean",
                           "plot",
                           "network"
                         )) {
  if (inherits(x, "ggplot")) {
    return(x)
  }

  if (is.list(x)) {
    for (nm in preferred_names) {
      if (!is.null(x[[nm]]) && inherits(x[[nm]], "ggplot")) {
        return(x[[nm]])
      }
    }

    for (i in seq_along(x)) {
      if (inherits(x[[i]], "ggplot")) {
        return(x[[i]])
      }
    }
  }

  stop("Could not extract a ggplot object from this result.")
}

safe_network_plot <- function(type_value,
                              superterms_value,
                              filename,
                              preferred_names,
                              width = 9,
                              height = 7) {
  message("\nTrying network plot: ", filename)

  plot_result <- tryCatch(
    {
      network_clust_gg(
        x                 = brca_pa_similarity,
        clust_result      = brca_pa_clustering,
        jaccard_threshold = 0.5,
        min_degree        = 5,
        superterms        = superterms_value,
        superterm_data    = if (superterms_value) net_results$superterms else NULL,
        type              = type_value,
        seed              = 174
      )
    },
    error = function(e) {
      message("Skipped ", filename, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(plot_result)) {
    return(invisible(NULL))
  }

  p <- tryCatch(
    extract_plot(plot_result, preferred_names = preferred_names),
    error = function(e) {
      message("Could not extract plot for ", filename, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(p)) {
    return(invisible(NULL))
  }

  save_plot_safely(
    p,
    filename = filename,
    width    = width,
    height   = height
  )

  invisible(p)
}

## -----------------------------
## 3. Single-comparison pathway plot
## -----------------------------

pa_single <- subset(
  brca_pa_merged,
  COMPARISON == "ERpositive_vs_ERnegative" &
    COLLECTION == "Hallmark"
)

if (nrow(pa_single) == 0) {
  stop("No Hallmark pathways found for ERpositive_vs_ERnegative.")
}

p_splot <- splot_PA(
  data           = pa_single,
  geneset_col    = "NAME",
  collection_col = "COLLECTION",
  nes_col        = "NES",
  fdr_col        = "FDR",
  order          = "desc",
  fill_limits    = c(0, 5),
  top_n_per_direction = 15,
  clean_labels   = TRUE,
  label_wrap     = 32,
  theme_params   = list(
    axis_title_size   = 12,
    axis_text_size_x  = 10,
    axis_text_size_y  = 7,
    strip_text_size   = 10,
    legend_title_size = 10,
    legend_text_size  = 9,
    bar_size          = 0.3,
    bar_width         = 0.72,
    hline_size        = 0.45
  )
)

save_plot_safely(
  p_splot,
  filename = "pathway_splot_er_status.pdf",
  width    = 8,
  height   = 8.5
)

## -----------------------------
## 4. Multi-comparison pathway plot
## -----------------------------

pa_multi_subset <- subset(
  brca_pa_merged,
  COMPARISON %in% c("Tumor_vs_Normal", "ERpositive_vs_ERnegative") &
    NAME %in% c(
      "HALLMARK_ESTROGEN_RESPONSE_EARLY",
      "HALLMARK_ESTROGEN_RESPONSE_LATE",
      "HALLMARK_P53_PATHWAY",
      "HALLMARK_IL2_STAT5_SIGNALING"
    )
)

if (nrow(pa_multi_subset) == 0) {
  stop("No data found for selected Hallmark pathways in pa_multi_subset.")
}

## Clean labels for the x-axis
comparison_labels <- c(
  "Tumor_vs_Normal" = "Tumor vs\nNormal",
  "ERpositive_vs_ERnegative" = "ER+ vs\nER-"
)

## Clean labels for facet strips
pa_multi_subset$PATHWAY_LABEL <- dplyr::recode(
  pa_multi_subset$NAME,
  "HALLMARK_ESTROGEN_RESPONSE_EARLY" = "Estrogen response\nearly",
  "HALLMARK_ESTROGEN_RESPONSE_LATE"  = "Estrogen response\nlate",
  "HALLMARK_P53_PATHWAY" = "p53 pathway",
  "HALLMARK_IL2_STAT5_SIGNALING" = "IL2-STAT5\nPathway"
)

## Keep a biologically meaningful order
pa_multi_subset$COMPARISON <- factor(
  pa_multi_subset$COMPARISON,
  levels = c("Tumor_vs_Normal", "ERpositive_vs_ERnegative")
)

pa_multi_subset$PATHWAY_LABEL <- factor(
  pa_multi_subset$PATHWAY_LABEL,
  levels = c(
    "Estrogen response\nearly",
    "Estrogen response\nlate",
    "p53 pathway","IL2-STAT5\nPathway"
  )
)

p_multi <- multiplot_PA(
  data           = pa_multi_subset,
  comparison_col = "COMPARISON",
  facet_col      = "PATHWAY_LABEL",
  axis_y         = "NES",
  fdr_col        = "FDR",
  comparison_order = c("Tumor_vs_Normal", "ERpositive_vs_ERnegative"),
  custom_labels = comparison_labels,
  ncol_wrap      = 2,
  free_y         = FALSE,
  fill_limits    = c(0, 5),
  theme_params   = list(
    bar_col             = "black",
    bar_size            = 0.35,
    bar_width           = 0.55,
    hline_size          = 0.6,
    axis_title_size     = 13,
    axis_text_size_x    = 10,
    axis_text_size_y    = 10,
    tick_size           = 0.4,
    tick_length         = 0.12,
    strip_text_size     = 11,
    panel_spacing_multi = 0.8
  )
) +
  ggplot2::labs(
    x = NULL,
    y = "Normalized enrichment score"
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      title.position = "top",
      barwidth = 0.5,
      barheight = 3.2
    )
  ) +
  ggplot2::theme(
    plot.margin = ggplot2::margin(10, 12, 10, 10),
    axis.text.x = ggplot2::element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    ),
    legend.title = ggplot2::element_text(size = 10),
    legend.text  = ggplot2::element_text(size = 9),
    legend.key.height = grid::unit(0.35, "cm"),
    strip.background = ggplot2::element_rect(fill = "grey90", color = "grey35"),
    panel.grid.minor = ggplot2::element_blank()
  )

save_plot_safely(
  p_multi,
  filename = "pathway_multiplot_estrogen_response.pdf",
  width    = 7.5,
  height   = 7.5
)

## -----------------------------
## 5. Silhouette plot
## -----------------------------

if (!is.null(brca_pa_clustering$silhouette_plot)) {
  save_plot_safely(
    brca_pa_clustering$silhouette_plot,
    filename = "pathway_silhouette_plot.pdf",
    width    = 8,
    height   = 5
  )
} else {
  warning("brca_pa_clustering$silhouette_plot is NULL. Skipping silhouette plot.")
}

## -----------------------------
## 6. ER-only network input
## -----------------------------
## The pre-computed brca_pa_similarity object contains multiple comparisons.
## Here we keep only the ERpositive_vs_ERnegative nodes.

target_network_comparison <- "ERpositive_vs_ERnegative"

subset_jaccard_result_by_comparison <- function(x,
                                                comparison,
                                                remove_comparison_prefix = TRUE,
                                                keep_collection_prefix = TRUE) {
  if (is.null(x$jaccard_sim)) {
    stop("The input object must contain a $jaccard_sim matrix.")
  }

  sim <- x$jaccard_sim

  if (is.null(rownames(sim)) || is.null(colnames(sim))) {
    stop("$jaccard_sim must have row and column names.")
  }

  keep <- grepl(paste0("^", comparison, "::"), rownames(sim))

  if (!any(keep)) {
    stop(
      "No nodes found for comparison: ", comparison, "\n",
      "First available node names:\n",
      paste(head(rownames(sim), 10), collapse = "\n")
    )
  }

  sim_sub <- sim[keep, keep, drop = FALSE]

  if (remove_comparison_prefix) {
    new_names <- sub(paste0("^", comparison, "::"), "", rownames(sim_sub))

    if (!keep_collection_prefix) {
      ## Keep only the last element after "::"
      new_names <- sub("^.*::", "", new_names)
    }

    new_names <- make.unique(new_names)

    rownames(sim_sub) <- new_names
    colnames(sim_sub) <- new_names
  }

  x_sub <- x
  x_sub$jaccard_sim <- sim_sub

  ## Rebuild distance matrix from the filtered similarity matrix.
  ## This avoids carrying distances from the two-comparison object.
  x_sub$dist_mat <- as.dist(1 - sim_sub)

  ## If the object has extra result tables, keep them only if possible.
  optional_table_slots <- c(
    "results",
    "filtered_results",
    "filtered",
    "results_filtered"
  )

  for (slot in optional_table_slots) {
    if (!is.null(x_sub[[slot]]) && is.data.frame(x_sub[[slot]])) {
      possible_cols <- intersect(
        c("COMPARISON", "comparison"),
        names(x_sub[[slot]])
      )

      if (length(possible_cols) > 0) {
        x_sub[[slot]] <- x_sub[[slot]][
          x_sub[[slot]][[possible_cols[1]]] == comparison,
          ,
          drop = FALSE
        ]
      }
    }
  }

  x_sub
}

brca_pa_similarity_er <- subset_jaccard_result_by_comparison(
  x                        = brca_pa_similarity,
  comparison               = target_network_comparison,
  remove_comparison_prefix = TRUE,
  keep_collection_prefix   = TRUE
)

message("\nER-only similarity matrix dimensions:")
print(dim(brca_pa_similarity_er$jaccard_sim))

message("\nFirst ER-only network node names:")
print(head(rownames(brca_pa_similarity_er$jaccard_sim), 10))

## Recompute clustering using only ERpositive_vs_ERnegative nodes.
brca_pa_clustering_er <- do_clust(brca_pa_similarity_er)

save_plot_safely(
  brca_pa_clustering_er$silhouette_plot,
  filename = "pathway_silhouette_plot_er_status.pdf",
  width    = 8,
  height   = 5
)

## Also overwrite the generic silhouette file used by the vignette.
save_plot_safely(
  brca_pa_clustering_er$silhouette_plot,
  filename = "pathway_silhouette_plot.pdf",
  width    = 8,
  height   = 5
)

## -----------------------------
## 7. ER-only network communities
## -----------------------------

net_results_er <- get_network_communities(
  x             = brca_pa_similarity_er,
  threshold     = 0.4,
  method        = "louvain",
  superterms    = TRUE,
  n_terms       = 4,
  remove_prefix = TRUE,
  seed          = 174
)

message("\nER-only network community summary:")
if (!is.null(net_results_er$superterms$summary)) {
  print(net_results_er$superterms$summary)
}

## -----------------------------
## 8. ER-only network plot variants
## -----------------------------

safe_network_plot_er <- function(type_value,
                                 superterms_value,
                                 filename,
                                 preferred_names,
                                 width = 9,
                                 height = 7) {
  message("\nTrying ER-only network plot: ", filename)

  plot_result <- tryCatch(
    {
      network_clust_gg(
        x                 = brca_pa_similarity_er,
        clust_result      = brca_pa_clustering_er,
        jaccard_threshold = 0.4,
        min_degree        = 4,
        superterms        = superterms_value,
        superterm_data    = if (superterms_value) net_results_er$superterms else NULL,
        type              = type_value,
        seed              = 174
      )
    },
    error = function(e) {
      message("Skipped ", filename, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(plot_result)) {
    return(invisible(NULL))
  }

  p <- tryCatch(
    extract_plot(plot_result, preferred_names = preferred_names),
    error = function(e) {
      message("Could not extract plot for ", filename, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(p)) {
    return(invisible(NULL))
  }

  save_plot_safely(
    p,
    filename = filename,
    width    = width,
    height   = height
  )

  invisible(p)
}

## Main ER-only network figure.
safe_network_plot_er(
  type_value       = "combined",
  superterms_value = TRUE,
  filename         = "pathway_network_communities_er_status.pdf",
  preferred_names  = c("combined", "superterms", "clean", "individual"),
  width            = 10,
  height           = 8
)

## Also overwrite the generic file used by the current vignette.
safe_network_plot_er(
  type_value       = "combined",
  superterms_value = TRUE,
  filename         = "pathway_network_communities.pdf",
  preferred_names  = c("combined", "superterms", "clean", "individual"),
  width            = 10,
  height           = 8
)

safe_network_plot_er(
  type_value       = "clean",
  superterms_value = FALSE,
  filename         = "pathway_network_er_clean_no_superterms.pdf",
  preferred_names  = c("clean", "plot", "network"),
  width            = 9,
  height           = 7
)

safe_network_plot_er(
  type_value       = "superterms",
  superterms_value = TRUE,
  filename         = "pathway_network_er_superterms.pdf",
  preferred_names  = c("superterms", "combined", "plot", "network"),
  width            = 9,
  height           = 7
)

safe_network_plot_er(
  type_value       = "individual",
  superterms_value = FALSE,
  filename         = "pathway_network_er_individual_no_superterms.pdf",
  preferred_names  = c("individual", "combined", "clean", "plot", "network"),
  width            = 11,
  height           = 9
)

safe_network_plot_er(
  type_value       = "individual",
  superterms_value = TRUE,
  filename         = "pathway_network_er_individual_with_superterms.pdf",
  preferred_names  = c("individual", "combined", "superterms", "plot", "network"),
  width            = 11,
  height           = 9
)

safe_network_plot_er(
  type_value       = "all",
  superterms_value = TRUE,
  filename         = "pathway_network_er_all_with_superterms.pdf",
  preferred_names  = c("combined", "superterms", "individual", "clean", "plot", "network"),
  width            = 10,
  height           = 8
)

message("\nDone. ER-only pathway network figures generated.")
