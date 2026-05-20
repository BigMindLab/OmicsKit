#########################
# Function trend_filter #
#########################

#' Filter genes by patient-level trend consistency.
#'
#' This function filters differentially expressed genes according to the
#' consistency of their expression trend across paired patients. For each gene
#' \eqn{g} and patient \eqn{p}, the function calculates the mean expression in
#' baseline samples and condition samples:
#' \eqn{\bar{X}_{N}(g,p)} and \eqn{\bar{X}_{T}(g,p)}.
#'
#' For genes classified as UP-regulated at the group level, a gene is removed if
#' at least one paired patient shows:
#' \eqn{\bar{X}_{N}(g,p) >= \bar{X}_{T}(g,p) * ratio}.
#'
#' For genes classified as DOWN-regulated at the group level, a gene is removed
#' if at least one paired patient shows:
#' \eqn{\bar{X}_{T}(g,p) >= \bar{X}_{N}(g,p) * ratio}.
#'
#' The direction of regulation is defined from the group-level log2 fold-change.
#' The patient-level consistency check is performed using the expression matrix
#' supplied in `expr`.
#'
#' @param expr Numeric matrix or data frame of expression values with genes as
#'   rows and sample IDs as columns. Row names must contain gene IDs.
#' @param brca_rna_metadata_tumor_normal Data frame with sample metadata.
#' @param results Data frame or named list of data frames containing
#'   differential expression results. Each data frame must contain a gene ID
#'   column and a log2 fold-change column.
#' @param baseline Character vector identifying the baseline group or groups in
#'   `group.col`.
#' @param conditions Character vector or named list identifying the condition
#'   group for each comparison in `results`. If `results` is a named list,
#'   names in `conditions` should match names in `results`.
#' @param sample.col Column in `brca_rna_metadata_tumor_normal` containing sample IDs. Default is
#'   `"sample_id"`.
#' @param patient.col Column in `brca_rna_metadata_tumor_normal` containing patient IDs. Default is
#'   `"patient_id"`.
#' @param group.col Column in `brca_rna_metadata_tumor_normal` containing group labels. Default is
#'   `"sample_type"`.
#' @param gene.col Column in `results` containing gene IDs. Default is
#'   `"ensembl"`.
#' @param lfc.col Column in `results` containing log2 fold changes. Default is
#'   `"log2FoldChange"`.
#' @param ratio Numeric value greater than 1. Defines the tolerated reversal
#'   threshold. Default is `1.1`.
#' @param lfc.cutoff Non-negative numeric value used to define UP and DOWN
#'   regulation from `lfc.col`. Default is `0`.
#' @param scale Expression scale. Use `"linear"` for normalized counts, CPM, TPM
#'   or similar linear-scale values. Use `"log2"` for log2-transformed
#'   expression values. Default is `"linear"`.
#' @param require_complete_pairs Logical. If `TRUE`, the function stops when
#'   unpaired patients are found. If `FALSE`, only patients with both baseline
#'   and condition samples are used. Default is `FALSE`.
#' @param na.rm Logical. Should missing expression values be removed when
#'   calculating patient-level means? Default is `FALSE`.
#' @param return_removed Logical. Should the output include a vector of removed
#'   genes? Default is `TRUE`.
#'
#' @return A named list containing:
#'   \itemize{
#'     \item One filtered data frame per comparison.
#'     \item `TrendGenes`: unique genes passing the trend consistency filter.
#'     \item `Diagnostics`: gene-level filtering diagnostics.
#'     \item `Summary`: comparison-level filtering summary.
#'     \item `RemovedGenes`: unique genes removed by the filter, if
#'       `return_removed = TRUE`.
#'   }
#'
#' @details
#' The multiplicative rule `ratio = 1.1` is appropriate for linear-scale
#' expression values. If `expr` is log2-transformed, set `scale = "log2"`; the
#' function will use `log2(ratio)` as an additive threshold.
#'
#' @references Requena D. et al. Nat Commun 15, 10887 (2024).
#'
#' @examples
#' \dontrun{
#' trend_res <- trend_filter(
#'   expr = brca_rna_expr_tumor_normal_filtered,
#'   brca_rna_metadata_tumor_normal = brca_rna_metadata_tumor_normal,
#'   results = list(Tumor_vs_Normal = deseq_res),
#'   baseline = "normal",
#'   conditions = c(Tumor_vs_Normal = "tumor"),
#'   sample.col = "sample_id",
#'   patient.col = "patient_id",
#'   group.col = "sample_type"
#' )
#'
#' length(trend_res$TrendGenes)
#' head(trend_res$Tumor_vs_Normal)
#' head(trend_res$Diagnostics)
#' trend_res$Summary
#' }
#'
#' @seealso [detect_filter()]
#'
#' @export

trend_filter <- function(expr,
                         sampledata,
                         results,
                         baseline,
                         conditions,
                         sample.col = "sample_id",
                         patient.col = "patient_id",
                         group.col = "sample_type",
                         gene.col = "ensembl",
                         lfc.col = "log2FoldChange",
                         ratio = 1.1,
                         lfc.cutoff = 0,
                         scale = c("linear", "log2"),
                         require_complete_pairs = FALSE,
                         na.rm = FALSE,
                         return_removed = TRUE) {

  scale <- match.arg(scale)

  expr <- as.data.frame(expr, check.names = FALSE)
  sampledata <- as.data.frame(sampledata)

  if (is.null(rownames(expr))) {
    stop("expr must have gene IDs as row names.")
  }

  if (anyDuplicated(rownames(expr)) > 0) {
    stop("expr must have unique gene IDs as row names.")
  }

  required_cols <- c(sample.col, patient.col, group.col)
  missing_cols <- setdiff(required_cols, colnames(sampledata))

  if (length(missing_cols) > 0) {
    stop(
      "The following columns are missing from sampledata: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (anyDuplicated(sampledata[[sample.col]]) > 0) {
    stop("sampledata must contain one row per sample ID.")
  }

  if (!is.character(baseline) || length(baseline) < 1) {
    stop("baseline must be a character vector with at least one group label.")
  }

  if (!is.numeric(ratio) || length(ratio) != 1 || is.na(ratio) || ratio <= 1) {
    stop("ratio must be a single numeric value greater than 1.")
  }

  if (
    !is.numeric(lfc.cutoff) ||
    length(lfc.cutoff) != 1 ||
    is.na(lfc.cutoff) ||
    lfc.cutoff < 0
  ) {
    stop("lfc.cutoff must be a single numeric value greater than or equal to 0.")
  }

  if (!is.logical(require_complete_pairs) || length(require_complete_pairs) != 1) {
    stop("require_complete_pairs must be TRUE or FALSE.")
  }

  if (!is.logical(na.rm) || length(na.rm) != 1) {
    stop("na.rm must be TRUE or FALSE.")
  }

  if (!is.logical(return_removed) || length(return_removed) != 1) {
    stop("return_removed must be TRUE or FALSE.")
  }

  if (inherits(results, "data.frame")) {
    results <- list(Comparison1 = results)
  }

  if (!is.list(results) || length(results) == 0) {
    stop("results must be a data frame or a non-empty list of data frames.")
  }

  valid_results <- vapply(results, inherits, logical(1), what = "data.frame")

  if (!all(valid_results)) {
    stop("Every element in results must be a data frame.")
  }

  if (is.null(names(results)) || any(names(results) == "")) {
    names(results) <- paste0("Comparison", seq_along(results))
  }

  reserved_names <- c("TrendGenes", "Diagnostics", "Summary", "RemovedGenes")

  if (any(names(results) %in% reserved_names)) {
    stop(
      "The following names are reserved and cannot be used as comparison names: ",
      paste(reserved_names, collapse = ", ")
    )
  }

  if (is.character(conditions)) {
    condition_names <- names(conditions)
    conditions <- as.list(conditions)

    if (!is.null(condition_names)) {
      names(conditions) <- condition_names
    }
  }

  if (!is.list(conditions) || length(conditions) != length(results)) {
    stop("conditions must have the same length as results.")
  }

  valid_conditions <- vapply(
    conditions,
    function(x) is.character(x) && length(x) >= 1,
    logical(1)
  )

  if (!all(valid_conditions)) {
    stop("Every element in conditions must be a character vector.")
  }

  if (is.null(names(conditions)) || any(names(conditions) == "")) {
    names(conditions) <- names(results)
  }

  if (!all(names(results) %in% names(conditions))) {
    stop("Names in conditions must match names in results.")
  }

  conditions <- conditions[names(results)]

  sampledata[[sample.col]] <- as.character(sampledata[[sample.col]])
  sampledata[[patient.col]] <- as.character(sampledata[[patient.col]])
  sampledata[[group.col]] <- as.character(sampledata[[group.col]])

  get_patient_means <- function(condition_label, comparison_name) {
    baseline_data <- sampledata[
      sampledata[[group.col]] %in% baseline,
      ,
      drop = FALSE
    ]

    condition_data <- sampledata[
      sampledata[[group.col]] %in% condition_label,
      ,
      drop = FALSE
    ]

    if (nrow(baseline_data) == 0) {
      stop("No baseline samples found for comparison: ", comparison_name)
    }

    if (nrow(condition_data) == 0) {
      stop("No condition samples found for comparison: ", comparison_name)
    }

    baseline_patients <- unique(baseline_data[[patient.col]])
    condition_patients <- unique(condition_data[[patient.col]])

    paired_patients <- intersect(baseline_patients, condition_patients)
    all_patients <- union(baseline_patients, condition_patients)
    unpaired_patients <- setdiff(all_patients, paired_patients)

    if (length(paired_patients) == 0) {
      stop(
        "No paired patients found between baseline and condition for comparison: ",
        comparison_name
      )
    }

    if (length(unpaired_patients) > 0 && require_complete_pairs) {
      stop(
        "Unpaired patients found in comparison ",
        comparison_name,
        ": ",
        paste(unpaired_patients, collapse = ", ")
      )
    }

    if (length(unpaired_patients) > 0 && !require_complete_pairs) {
      warning(
        "Using only paired patients for comparison ",
        comparison_name,
        ". Unpaired patients were ignored: ",
        paste(unpaired_patients, collapse = ", "),
        call. = FALSE
      )
    }

    baseline_data <- baseline_data[
      baseline_data[[patient.col]] %in% paired_patients,
      ,
      drop = FALSE
    ]

    condition_data <- condition_data[
      condition_data[[patient.col]] %in% paired_patients,
      ,
      drop = FALSE
    ]

    used_samples <- unique(c(
      baseline_data[[sample.col]],
      condition_data[[sample.col]]
    ))

    missing_samples <- setdiff(used_samples, colnames(expr))

    if (length(missing_samples) > 0) {
      stop(
        "The following samples are missing from expr: ",
        paste(missing_samples, collapse = ", ")
      )
    }

    non_numeric_samples <- used_samples[
      !vapply(expr[, used_samples, drop = FALSE], is.numeric, logical(1))
    ]

    if (length(non_numeric_samples) > 0) {
      stop(
        "The following samples in expr are not numeric: ",
        paste(non_numeric_samples, collapse = ", ")
      )
    }

    mean_baseline <- vapply(
      paired_patients,
      function(patient) {
        samples <- baseline_data[
          baseline_data[[patient.col]] == patient,
          sample.col
        ]

        rowMeans(expr[, samples, drop = FALSE], na.rm = na.rm)
      },
      numeric(nrow(expr))
    )

    mean_condition <- vapply(
      paired_patients,
      function(patient) {
        samples <- condition_data[
          condition_data[[patient.col]] == patient,
          sample.col
        ]

        rowMeans(expr[, samples, drop = FALSE], na.rm = na.rm)
      },
      numeric(nrow(expr))
    )

    rownames(mean_baseline) <- rownames(expr)
    rownames(mean_condition) <- rownames(expr)

    colnames(mean_baseline) <- paired_patients
    colnames(mean_condition) <- paired_patients

    list(
      baseline = mean_baseline,
      condition = mean_condition,
      patients = paired_patients
    )
  }

  filter_one_comparison <- function(df, condition_label, comparison_name) {
    df <- as.data.frame(df)

    if (!gene.col %in% colnames(df)) {
      stop(
        "Column '",
        gene.col,
        "' was not found in comparison: ",
        comparison_name
      )
    }

    if (!lfc.col %in% colnames(df)) {
      stop(
        "Column '",
        lfc.col,
        "' was not found in comparison: ",
        comparison_name
      )
    }

    genes <- as.character(df[[gene.col]])

    lfc <- suppressWarnings(as.numeric(df[[lfc.col]]))

    invalid_lfc <- is.na(lfc) & !is.na(df[[lfc.col]])

    if (any(invalid_lfc)) {
      stop(
        "Column '",
        lfc.col,
        "' must be numeric in comparison: ",
        comparison_name
      )
    }

    direction <- rep("none", length(lfc))
    direction[!is.na(lfc) & lfc > lfc.cutoff] <- "up"
    direction[!is.na(lfc) & lfc < -lfc.cutoff] <- "down"

    missing_gene <- !genes %in% rownames(expr)

    diagnostics <- data.frame(
      gene = genes,
      comparison = comparison_name,
      log2FoldChange = lfc,
      direction = direction,
      passed_trend = FALSE,
      reason = NA_character_,
      failed_patients = NA_character_,
      stringsAsFactors = FALSE
    )

    diagnostics$reason[is.na(lfc)] <- "missing_lfc"
    diagnostics$reason[direction == "none" & !is.na(lfc)] <- "no_direction"
    diagnostics$reason[missing_gene] <- "missing_in_expr"

    evaluable <- which(!missing_gene & direction %in% c("up", "down"))

    if (length(evaluable) == 0) {
      filtered_df <- df[FALSE, , drop = FALSE]

      return(list(
        filtered = filtered_df,
        diagnostics = diagnostics
      ))
    }

    patient_means <- get_patient_means(condition_label, comparison_name)

    genes_eval <- genes[evaluable]

    mean_baseline <- patient_means$baseline[genes_eval, , drop = FALSE]
    mean_condition <- patient_means$condition[genes_eval, , drop = FALSE]
    patients <- patient_means$patients

    missing_expression <- is.na(mean_baseline) | is.na(mean_condition)
    has_missing_expression <- rowSums(missing_expression) > 0

    if (scale == "linear") {
      inconsistent_up <- mean_baseline >= mean_condition * ratio
      inconsistent_down <- mean_condition >= mean_baseline * ratio
    } else {
      log2_ratio <- log2(ratio)

      inconsistent_up <- mean_baseline >= mean_condition + log2_ratio
      inconsistent_down <- mean_condition >= mean_baseline + log2_ratio
    }

    direction_eval <- direction[evaluable]

    has_inconsistent_trend <- ifelse(
      direction_eval == "up",
      rowSums(inconsistent_up, na.rm = TRUE) > 0,
      rowSums(inconsistent_down, na.rm = TRUE) > 0
    )

    passed <- !has_missing_expression & !has_inconsistent_trend

    reason <- ifelse(
      passed,
      "passed",
      ifelse(
        has_missing_expression,
        "missing_expression",
        "inconsistent_trend"
      )
    )

    failed_patients <- vapply(
      seq_along(evaluable),
      function(i) {
        if (has_missing_expression[i]) {
          patient_ids <- patients[missing_expression[i, ]]
        } else if (direction_eval[i] == "up") {
          patient_ids <- patients[inconsistent_up[i, ]]
        } else {
          patient_ids <- patients[inconsistent_down[i, ]]
        }

        if (length(patient_ids) == 0) {
          NA_character_
        } else {
          paste(patient_ids, collapse = ";")
        }
      },
      character(1)
    )

    diagnostics$passed_trend[evaluable] <- passed
    diagnostics$reason[evaluable] <- reason
    diagnostics$failed_patients[evaluable] <- failed_patients

    filtered_df <- df[diagnostics$passed_trend, , drop = FALSE]

    list(
      filtered = filtered_df,
      diagnostics = diagnostics
    )
  }

  filtered_results <- list()
  diagnostics_list <- list()

  for (comparison_name in names(results)) {
    filtered <- filter_one_comparison(
      df = results[[comparison_name]],
      condition_label = conditions[[comparison_name]],
      comparison_name = comparison_name
    )

    filtered_results[[comparison_name]] <- filtered$filtered
    diagnostics_list[[comparison_name]] <- filtered$diagnostics
  }

  all_diagnostics <- do.call(rbind, diagnostics_list)
  rownames(all_diagnostics) <- NULL

  trend_genes <- unique(unlist(
    lapply(
      filtered_results,
      function(x) as.character(x[[gene.col]])
    ),
    use.names = FALSE
  ))

  trend_genes <- as.character(trend_genes)

  summary_table <- do.call(
    rbind,
    lapply(
      split(all_diagnostics, all_diagnostics$comparison),
      function(x) {
        data.frame(
          comparison = x$comparison[1],
          total_genes = nrow(x),
          kept_genes = sum(x$passed_trend, na.rm = TRUE),
          removed_genes = sum(!x$passed_trend, na.rm = TRUE),
          up_genes = sum(x$direction == "up", na.rm = TRUE),
          down_genes = sum(x$direction == "down", na.rm = TRUE),
          no_direction = sum(x$reason == "no_direction", na.rm = TRUE),
          missing_lfc = sum(x$reason == "missing_lfc", na.rm = TRUE),
          missing_in_expr = sum(x$reason == "missing_in_expr", na.rm = TRUE),
          inconsistent_trend = sum(x$reason == "inconsistent_trend", na.rm = TRUE),
          missing_expression = sum(x$reason == "missing_expression", na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(summary_table) <- NULL

  output <- filtered_results
  output$TrendGenes <- trend_genes
  output$Diagnostics <- all_diagnostics
  output$Summary <- summary_table

  if (return_removed) {
    output$RemovedGenes <- unique(as.character(
      all_diagnostics$gene[!all_diagnostics$passed_trend]
    ))
  }

  return(output)
}
