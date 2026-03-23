## data-raw/gsea_results.R
## Generates gsea_results a simulated merge_PA() output for three comparisons
## (TumorVsNormal, MetastasisVsNormal, and MetastasisVsTumor) across three
## MSigDB collections (HALLMARK, KEGG, GO). Uses geneset_list and deseq2_results
## already in the package.
## Run once to regenerate data/gsea_results.rda
## Requires: OmicsKit data objects already loaded


# Setup
# =============================================================================
set.seed(174)

data(geneset_list)
data(deseq2_results)

# ranked_genes: genes ordered by DESeq2 Wald stat (most positive = most upregulated)
ranked_genes <- deseq2_results$gene_id[
  order(deseq2_results$stat, decreasing = TRUE)
]

# Define gene sets per collection
# We use names from geneset_list split by prefix
# =============================================================================
gs_names   <- names(geneset_list)
hallmark   <- gs_names[grepl("^HALLMARK_",  gs_names)]
kegg       <- gs_names[grepl("^KEGG_",      gs_names)]
go         <- gs_names[grepl("^GO_",        gs_names)]

# Use all available from each collection
collection_map <- c(
  setNames(rep("HALLMARK", length(hallmark)), hallmark),
  setNames(rep("KEGG",     length(kegg)),     kegg),
  setNames(rep("GO",       length(go)),       go)
)

all_sets    <- names(collection_map)
n_sets      <- length(all_sets)

# Helper: simulate one comparison
# =============================================================================
simulate_comparison <- function(comparison_name, seed_offset = 0) {

  set.seed(174 + seed_offset)

  sizes <- vapply(all_sets, function(gs) length(geneset_list[[gs]]), integer(1))

  # Simulate NES: mix of up/down enrichment
  nes <- rnorm(n_sets, mean = 0, sd = 1.5)

  # Simulate FDR directly: ~60% significant (FDR < 0.05), rest non-significant
  n_sig    <- round(n_sets * 0.6)
  n_nonsig <- n_sets - n_sig
  fdr_sig    <- sort(runif(n_sig,    min = 0.001, max = 0.049))
  fdr_nonsig <- runif(n_nonsig, min = 0.06,  max = 0.95)
  rank_by_abs <- order(abs(nes), decreasing = TRUE)
  fdr           <- numeric(n_sets)
  fdr[rank_by_abs[seq_len(n_sig)]]             <- fdr_sig
  fdr[rank_by_abs[seq(n_sig + 1L, n_sets)]]    <- fdr_nonsig
  nom_pval <- pmin(1, fdr * runif(n_sets, 0.3, 0.9))
  fwer     <- pmin(1, fdr * 1.5)

  # Simulate rank at max
  rank_at_max <- sample(seq_len(length(ranked_genes)), n_sets, replace = TRUE)

  # Simulate leading edge: tags% ~ |NES| / 3, capped at 60%
  tags_pct   <- pmin(0.60, abs(nes) / 3 + runif(n_sets, 0, 0.1))
  list_pct   <- runif(n_sets, 0.20, 0.80)
  signal_pct <- tags_pct * (1 - list_pct) / (1 - tags_pct * list_pct + 1e-6)
  signal_pct <- pmin(signal_pct, 1)

  leading_edge <- sprintf(
    "tags=%.0f%%, list=%.0f%%, signal=%.0f%%",
    tags_pct   * 100,
    list_pct   * 100,
    signal_pct * 100
  )

  data.frame(
    NAME           = all_sets,
    SIZE           = sizes,
    ES             = nes * 0.7,
    NES            = nes,
    `NOM p-val`    = nom_pval,
    FDR            = fdr,
    `FWER p-val`   = fwer,
    `RANK AT MAX`  = rank_at_max,
    Log10FDR       = -log10(fdr),
    tags           = tags_pct,
    list           = list_pct,
    signal         = signal_pct,
    `LEADING EDGE` = leading_edge,
    COLLECTION     = unname(collection_map[all_sets]),
    COMPARISON     = comparison_name,
    stringsAsFactors = FALSE,
    check.names    = FALSE
  )
}

# Generate three comparisons
# =============================================================================
comp1 <- simulate_comparison("TumorVsNormal",       seed_offset = 0)
comp2 <- simulate_comparison("MetastasisVsNormal",   seed_offset = 7)
comp3 <- simulate_comparison("MetastasisVsTumor",    seed_offset = 14)

gsea_results <- rbind(comp1, comp2, comp3)

cat("gsea_results:\n")
cat("  Rows          :", nrow(gsea_results), "\n")
cat("  Comparisons   :", paste(unique(gsea_results$COMPARISON), collapse = ", "), "\n")
cat("  Collections   :", paste(unique(gsea_results$COLLECTION), collapse = ", "), "\n")
cat("  Gene sets     :", n_sets, "\n")
cat("  FDR < 0.05    :",
    sum(gsea_results$FDR < 0.05), "out of", nrow(gsea_results), "\n")
cat("  Columns       :", paste(colnames(gsea_results), collapse = ", "), "\n")

# Save
# =============================================================================
usethis::use_data(gsea_results, compress = "xz", overwrite = TRUE)

message("Done. Saved: data/gsea_results.rda")
