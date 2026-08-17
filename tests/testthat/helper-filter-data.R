# Helper fixtures for detectability_filter() and trend_filter() tests.
# These fixtures are intentionally small and deterministic.

make_detect_fixture <- function() {
  norm_counts <- data.frame(
    N1 = c(10, 100, 100, 1, 20, 80, 10),
    N2 = c(20, 120, 110, 2, 25, 90, 12),
    B1 = c(100, 10, 100, 20, 5, 80, 15),
    B2 = c(120, 20, 110, 25, 10, 90, 16),
    C1 = c(15, 15, 100, 10, 10, 80, 100),
    C2 = c(20, 20, 110, 10, 10, 90, 120),
    D1 = c(15, 15, 100, 10, 10, 80, 5),
    D2 = c(20, 20, 110, 10, 10, 90, 10),
    check.names = FALSE
  )

  rownames(norm_counts) <- c(
    "g_up_detect",
    "g_down_detect",
    "g_low_baseMean",
    "g_up_low_expr",
    "g_down_low_baseline",
    "g_no_direction",
    "g_c_detect"
  )

  df_b <- data.frame(
    ensembl = rownames(norm_counts)[1:6],
    baseMean = c(100, 120, 20, 100, 100, 100),
    log2FoldChange = c(2, -2, 2, 2, -2, 0),
    padj = c(0.01, 0.01, 0.01, 0.02, 0.03, 0.04),
    stringsAsFactors = FALSE
  )

  df_c <- data.frame(
    ensembl = c("g_c_detect", "g_up_low_expr", "g_no_direction"),
    baseMean = c(100, 100, 100),
    log2FoldChange = c(2, 2, 0),
    padj = c(0.01, 0.02, 0.05),
    stringsAsFactors = FALSE
  )

  df_d <- data.frame(
    ensembl = c("g_down_detect", "g_c_detect"),
    baseMean = c(100, 100),
    log2FoldChange = c(-2, 2),
    padj = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  list(
    norm_counts = norm_counts,
    df_b = df_b,
    df_c = df_c,
    df_d = df_d,
    samples_baseline = c("N1", "N2"),
    samples_condition1 = c("B1", "B2"),
    samples_condition2 = c("C1", "C2"),
    samples_condition3 = c("D1", "D2")
  )
}

make_trend_fixture <- function(log_scale = FALSE) {
  expr_linear <- data.frame(
    N_P1 = c(10, 30, 30, 10, 50),
    T_P1 = c(20, 20, 10, 12, 60),
    N_P2 = c(12, 10, 40, 10, 50),
    T_P2 = c(24, 20, 15, 20, 55),
    check.names = FALSE
  )

  rownames(expr_linear) <- c(
    "g_up_pass",
    "g_up_fail",
    "g_down_pass",
    "g_down_fail",
    "g_no_direction"
  )

  expr <- if (log_scale) {
    as.data.frame(log2(expr_linear))
  } else {
    expr_linear
  }

  sampledata <- data.frame(
    sample_id = c("N_P1", "T_P1", "N_P2", "T_P2"),
    patient_id = c("P1", "P1", "P2", "P2"),
    sample_type = c("normal", "tumor", "normal", "tumor"),
    stringsAsFactors = FALSE
  )

  res <- data.frame(
    ensembl = c(
      "g_up_pass",
      "g_up_fail",
      "g_down_pass",
      "g_down_fail",
      "g_no_direction",
      "g_missing_expr"
    ),
    log2FoldChange = c(1.5, 1.2, -1.4, -1.1, 0, 2),
    padj = c(0.01, 0.01, 0.02, 0.03, 0.20, 0.01),
    stringsAsFactors = FALSE
  )

  list(
    expr = expr,
    sampledata = sampledata,
    res = res
  )
}

make_trend_unpaired_fixture <- function() {
  fixture <- make_trend_fixture()

  fixture$expr$N_P3 <- c(10, 10, 10, 10, 10)

  fixture$sampledata <- rbind(
    fixture$sampledata,
    data.frame(
      sample_id = "N_P3",
      patient_id = "P3",
      sample_type = "normal",
      stringsAsFactors = FALSE
    )
  )

  fixture
}
