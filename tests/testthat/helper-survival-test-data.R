# Helper objects used by the survival plotting/modeling tests.
# Place this file in tests/testthat/.

make_cox_test_data <- function(n = 120, seed = 123) {
  set.seed(seed)

  gene_status <- factor(
    rep(c("WT", "Mut"), each = n / 2),
    levels = c("WT", "Mut")
  )
  age <- round(seq(42, 78, length.out = n) + rnorm(n, sd = 3), 1)
  x <- as.integer(gene_status == "Mut")

  true_time <- stats::rexp(n, rate = 0.004 * exp(0.45 * x + 0.015 * (age - 60)))
  censor_time <- stats::rexp(n, rate = 0.0025)

  data.frame(
    sample = paste0("S", seq_len(n)),
    PFI.time = round(pmin(true_time, censor_time) * 100, 2),
    PFI = as.integer(true_time <= censor_time),
    gene_status = gene_status,
    age = age,
    protein_score = stats::rnorm(n),
    sparse_var = factor(c(rep("Rare", 3), rep("Common", n - 3))),
    days_to_death = round(true_time * 100, 2),
    `_hidden_omic` = stats::rnorm(n),
    check.names = FALSE
  )
}

make_forest_test_data <- function() {
  data.frame(
    model = "adjusted",
    model_id = c("PAM50 adjusted", "DNAmethyl adjusted", "RPPA adjusted"),
    variable = c("PAM50Call_RNAseq", "_PANCAN_DNAMethyl_BRCA", "RPPA_Clusters_nature2012"),
    term_clean = c("Her2", "cluster 3", "Reactive"),
    reference = c("LumA", "cluster 4", "Basal"),
    HR = c(2.25, 0.43, 0.22),
    CI_low = c(1.25, 0.19, 0.06),
    CI_high = c(4.07, 0.98, 0.78),
    p.value = c(0.00715, 0.0435, 0.0186),
    n_used = c(100, 100, 100),
    n_events = c(30, 30, 30),
    adjusted_for = "age",
    stringsAsFactors = FALSE
  )
}

make_km_test_data <- function(n = 60) {
  data.frame(
    PFI.time = c(seq(25, 750, length.out = n / 2), seq(40, 900, length.out = n / 2)),
    PFI = rep(c(1, 0, 1, 1, 0, 0), length.out = n),
    GENE_muts = factor(rep(c("No", "Yes"), each = n / 2), levels = c("No", "Yes")),
    PAM50_group = factor(
      rep(c("LumA", "HER2", "Basal"), length.out = n),
      levels = c("LumA", "HER2", "Basal")
    )
  )
}
