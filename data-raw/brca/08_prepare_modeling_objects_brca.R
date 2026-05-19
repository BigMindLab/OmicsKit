## 08_prepare_modeling_objects_brca.R
##
## Purpose:
##   Prepare real TCGA-BRCA clinical-only modeling objects for nice_KM(),
##   get_cox(), nice_forest(), get_glm(), and nice_ROC().
##
## Inputs:
##   - inst/extdata/brca/raw_xena/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix
##   - inst/extdata/brca/raw_xena/BRCA_survival.txt
##
## Outputs:
##   Package objects:
##   - brca_clinical_modeling_data
##   - brca_cox_univariable_clinical
##   - brca_cox_adjusted_clinical
##   - brca_clinical_level_summary_pfi
##   - brca_clinical_ref_template
##   - brca_glm_stage_clinical
##   - brca_glm_stage_train
##   - brca_glm_stage_test
##   - brca_glm_stage_predictions
##   - brca_roc_stage_clinical
##
## External files:
##   - inst/extdata/brca/intermediate/brca_clinical_modeling_data.tsv
##   - inst/extdata/brca/intermediate/brca_cox_univariable_clinical.tsv
##   - inst/extdata/brca/intermediate/brca_cox_adjusted_clinical.tsv
##   - inst/extdata/brca/intermediate/brca_glm_stage_predictions.tsv
##
## Rules:
##   Clinical data only. Do not compute or include omics_score, omics PCs,
##   PAM50, RPPA clusters, methylation clusters, CN clusters, miRNA clusters,
##   or any PANCAN molecular cluster variables.

options(stringsAsFactors = FALSE)

project_file <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    here::here(...)
  } else {
    file.path(...)
  }
}

raw_dir <- project_file("inst", "extdata", "brca", "raw_xena")
intermediate_dir <- project_file("inst", "extdata", "brca", "intermediate")

clinical_file <- file.path(raw_dir, "TCGA.BRCA.sampleMap_BRCA_clinicalMatrix")
survival_file <- file.path(raw_dir, "BRCA_survival.txt")

forbidden_variables <- c(
  "omics_score",
  "RNA PC1",
  "methylation PC1",
  "PAM50Call_RNAseq",
  "PAM50_mRNA_nature2012",
  "CN_Clusters_nature2012",
  "Integrated_Clusters_with_PAM50__nature2012",
  "Integrated_Clusters_no_exp__nature2012",
  "RPPA_Clusters_nature2012",
  "methylation_Clusters_nature2012",
  "miRNA_Clusters_nature2012",
  "_PANCAN_CNA_PANCAN_K8",
  "_PANCAN_DNAMethyl_BRCA",
  "_PANCAN_miRNA_PANCAN",
  "_PANCAN_UNC_RNAseq_PANCAN_K16"
)

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "Required raw file is missing for BRCA clinical modeling: ", label, "\n",
      "Expected location: ", normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }

  invisible(path)
}

require_file(clinical_file, "TCGA.BRCA.sampleMap_BRCA_clinicalMatrix")
require_file(survival_file, "BRCA_survival.txt")

clean_empty_as_na <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (!is.character(x)) {
    return(x)
  }

  y <- trimws(x)
  y[y == ""] <- NA_character_
  y[tolower(y) %in% c(
    "na", "n/a", "nan", "null", "none", "unknown", "not available",
    "not reported", "not applicable", "[not available]", "[not applicable]",
    "[unknown]", "--"
  )] <- NA_character_
  y
}

has_value <- function(x) {
  y <- clean_empty_as_na(as.character(x))
  !is.na(y) & nzchar(y)
}

clean_tcga_barcode <- function(x, level = c("sample", "patient")) {
  level <- match.arg(level)
  y <- toupper(trimws(as.character(x)))
  y <- gsub("\\.", "-", y)
  y[!has_value(y)] <- NA_character_

  if (identical(level, "patient")) {
    return(ifelse(!is.na(y) & nchar(y) >= 12L, substr(y, 1L, 12L), NA_character_))
  }

  ifelse(
    !is.na(y) & nchar(y) >= 16L,
    substr(y, 1L, 16L),
    ifelse(!is.na(y) & nchar(y) >= 15L, substr(y, 1L, 15L), NA_character_)
  )
}

make_age_group <- function(age) {
  age_num <- suppressWarnings(as.numeric(as.character(age)))
  out <- rep(NA_character_, length(age_num))
  out[is.finite(age_num) & age_num < 50] <- "<50"
  out[is.finite(age_num) & age_num >= 50] <- ">=50"
  factor(out, levels = c("<50", ">=50"))
}

make_stage_broad <- function(stage) {
  x <- toupper(clean_empty_as_na(as.character(stage)))
  x <- gsub("STAGE", "", x, fixed = TRUE)
  x <- trimws(x)

  out <- rep(NA_character_, length(x))
  out[grepl("^I($|[A-C])", x)] <- "Stage I"
  out[grepl("^II($|[A-C])", x)] <- "Stage II"
  out[grepl("^III($|[A-C])", x)] <- "Stage III"
  out[grepl("^IV($|[A-C])", x)] <- "Stage IV"
  factor(out, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
}

make_stage_advanced <- function(stage_broad) {
  x <- as.character(stage_broad)
  out <- rep(NA_integer_, length(x))
  out[x %in% c("Stage I", "Stage II")] <- 0L
  out[x %in% c("Stage III", "Stage IV")] <- 1L
  out
}

inspect_levels <- function(data, vars, event_col = "PFI", min_n = 10, min_events = 5) {
  rows <- lapply(vars, function(v) {
    if (!v %in% names(data)) {
      return(data.frame(
        variable = v,
        level = NA_character_,
        n = 0L,
        n_events = 0L,
        eligible = FALSE,
        reason = "missing variable",
        stringsAsFactors = FALSE
      ))
    }

    x <- data[[v]]
    event <- suppressWarnings(as.numeric(as.character(data[[event_col]])))
    keep <- !is.na(x) & is.finite(event)

    if (is.numeric(x) || is.integer(x)) {
      n <- sum(keep)
      n_events <- sum(event[keep] == 1, na.rm = TRUE)
      return(data.frame(
        variable = v,
        level = "continuous",
        n = n,
        n_events = n_events,
        eligible = n >= min_n && n_events >= min_events && length(unique(x[keep])) >= 2L,
        reason = if (n < min_n) "too few observations" else if (n_events < min_events) "too few events" else if (length(unique(x[keep])) < 2L) "no variation" else "ok",
        stringsAsFactors = FALSE
      ))
    }

    x <- droplevels(factor(x[keep]))
    event <- event[keep]

    if (length(x) == 0L) {
      return(data.frame(
        variable = v,
        level = NA_character_,
        n = 0L,
        n_events = 0L,
        eligible = FALSE,
        reason = "no complete observations",
        stringsAsFactors = FALSE
      ))
    }

    tab <- table(x)
    event_tab <- tapply(event == 1, x, sum, na.rm = TRUE)
    event_tab <- event_tab[names(tab)]
    event_tab[is.na(event_tab)] <- 0

    data.frame(
      variable = v,
      level = names(tab),
      n = as.integer(tab),
      n_events = as.integer(event_tab),
      eligible = as.integer(tab) >= min_n & as.integer(event_tab) >= min_events,
      reason = ifelse(as.integer(tab) < min_n, "too few observations",
                      ifelse(as.integer(event_tab) < min_events, "too few events", "ok")),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

make_reference_template <- function(data, vars) {
  rows <- lapply(vars, function(v) {
    if (!v %in% names(data)) {
      return(NULL)
    }

    x <- data[[v]]
    if (is.numeric(x) || is.integer(x)) {
      reference <- "continuous"
    } else {
      x <- droplevels(factor(x))
      reference <- if (nlevels(x) > 0L) levels(x)[1L] else NA_character_
    }

    data.frame(
      variable = v,
      reference = reference,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

require_columns <- function(data, columns, label = "data") {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stop(
      label, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

safe_factor <- function(x, ref = NULL, min_n = 5) {
  x <- clean_empty_as_na(as.character(x))
  tab <- sort(table(x), decreasing = TRUE)
  keep_levels <- names(tab)[tab >= min_n]
  x[!is.na(x) & !x %in% keep_levels] <- NA_character_
  out <- droplevels(factor(x))

  if (!is.null(ref) && ref %in% levels(out)) {
    out <- stats::relevel(out, ref = ref)
  }

  out
}

normalize_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

find_column <- function(data, candidates, label, required = TRUE) {
  data_names <- names(data)
  data_norm <- normalize_name(data_names)
  candidate_norm <- normalize_name(candidates)
  matched <- match(candidate_norm, data_norm)

  if (any(!is.na(matched))) {
    return(data_names[matched[which(!is.na(matched))[1L]]])
  }

  if (isTRUE(required)) {
    stop(
      "Could not find required column for ", label, ". Expected one of: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }

  NULL
}

add_optional_column <- function(out, source, output_name, candidates, transform = identity) {
  source_col <- find_column(source, candidates, output_name, required = FALSE)
  if (!is.null(source_col)) {
    out[[output_name]] <- transform(source[[source_col]])
  }
  out
}

read_xena_table <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
}

standardize_binary_event <- function(x) {
  y <- clean_empty_as_na(as.character(x))
  lower <- tolower(y)
  out <- suppressWarnings(as.numeric(y))
  out[lower %in% c("dead", "deceased", "event", "yes", "true", "progressed", "recurred")] <- 1
  out[lower %in% c("alive", "censored", "no", "false", "diseasefree", "disease_free")] <- 0
  out
}

drop_forbidden_columns <- function(data) {
  forbidden_norm <- normalize_name(forbidden_variables)
  data_norm <- normalize_name(names(data))
  forbidden_present <- names(data)[data_norm %in% forbidden_norm]

  if (length(forbidden_present) > 0L) {
    warning(
      "Dropping forbidden non-clinical or leakage variables: ",
      paste(forbidden_present, collapse = ", "),
      call. = FALSE
    )
    data <- data[, setdiff(names(data), forbidden_present), drop = FALSE]
  }

  data
}

source(project_file("R", "get_cox.R"))
source(project_file("R", "get_glm.R"))
source(project_file("R", "nice_ROC.R"))

for (pkg in c("survival", "broom", "dplyr", "tibble", "pROC", "ggplot2", "usethis")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package \"", pkg, "\" is required for BRCA clinical modeling.", call. = FALSE)
  }
}

phenotype <- read_xena_table(clinical_file)
survival <- read_xena_table(survival_file)

phenotype <- drop_forbidden_columns(phenotype)
survival <- drop_forbidden_columns(survival)

phenotype_id_col <- find_column(
  phenotype,
  c("sampleID", "sample", "Samples", "Sample", "ID"),
  "phenotype sampleID"
)
survival_sample_col <- find_column(
  survival,
  c("sample", "sampleID", "Samples", "Sample", "patient", "Patient", "_PATIENT"),
  "survival sample"
)

survival_os_col <- find_column(survival, c("OS"), "OS")
survival_os_time_col <- find_column(survival, c("OS.time", "OS_time", "OS Time"), "OS.time")
survival_pfi_col <- find_column(survival, c("PFI"), "PFI")
survival_pfi_time_col <- find_column(survival, c("PFI.time", "PFI_time", "PFI Time"), "PFI.time")

phenotype$sample <- clean_tcga_barcode(phenotype[[phenotype_id_col]], level = "sample")
phenotype$sampleID <- phenotype$sample
phenotype$patient <- clean_tcga_barcode(phenotype[[phenotype_id_col]], level = "patient")

survival$sample <- clean_tcga_barcode(survival[[survival_sample_col]], level = "sample")
survival$patient <- clean_tcga_barcode(survival[[survival_sample_col]], level = "patient")

merge_key <- if (any(!is.na(survival$sample)) && any(!is.na(phenotype$sample))) {
  "sample"
} else {
  "patient"
}

survival_small <- survival[, c(
  merge_key,
  survival_os_col,
  survival_os_time_col,
  survival_pfi_col,
  survival_pfi_time_col
), drop = FALSE]
names(survival_small) <- c(
  merge_key,
  ".survival_OS",
  ".survival_OS_time",
  ".survival_PFI",
  ".survival_PFI_time"
)
survival_small <- survival_small[!is.na(survival_small[[merge_key]]), , drop = FALSE]
survival_small <- survival_small[!duplicated(survival_small[[merge_key]]), , drop = FALSE]

merged <- merge(
  phenotype,
  survival_small,
  by = merge_key,
  all.x = FALSE,
  all.y = FALSE,
  sort = FALSE
)

if (nrow(merged) == 0L) {
  stop(
    "Clinical phenotype and survival data did not merge by ", merge_key, ".",
    call. = FALSE
  )
}

if (!"sample" %in% names(merged)) {
  merged$sample <- clean_tcga_barcode(merged[[phenotype_id_col]], level = "sample")
}
if (!"sampleID" %in% names(merged)) {
  merged$sampleID <- merged$sample
}
if (!"patient" %in% names(merged)) {
  merged$patient <- clean_tcga_barcode(merged[[phenotype_id_col]], level = "patient")
}

age_col <- find_column(
  merged,
  c(
    "age_at_initial_pathologic_diagnosis",
    "Age_at_Initial_Pathologic_Diagnosis_nature2012",
    "age_at_diagnosis",
    "age"
  ),
  "age_at_initial_pathologic_diagnosis"
)

brca_clinical_modeling_data <- data.frame(
  sample = merged$sample,
  sampleID = merged$sampleID,
  patient = merged$patient,
  OS = standardize_binary_event(merged$.survival_OS),
  OS.time = suppressWarnings(as.numeric(as.character(merged$.survival_OS_time))),
  PFI = standardize_binary_event(merged$.survival_PFI),
  PFI.time = suppressWarnings(as.numeric(as.character(merged$.survival_PFI_time))),
  age_at_initial_pathologic_diagnosis = suppressWarnings(as.numeric(as.character(merged[[age_col]]))),
  stringsAsFactors = FALSE
)

brca_clinical_modeling_data$age_group <- make_age_group(
  brca_clinical_modeling_data$age_at_initial_pathologic_diagnosis
)

optional_clinical <- list(
  ER_Status_nature2012 = c("ER_Status_nature2012", "ER_status", "ER Status", "estrogen_receptor_status"),
  PR_Status_nature2012 = c("PR_Status_nature2012", "PR_status", "PR Status", "progesterone_receptor_status"),
  HER2_Final_Status_nature2012 = c("HER2_Final_Status_nature2012", "HER2_Status_nature2012", "HER2_status", "HER2 Status", "her2_neu_status"),
  AJCC_Stage_nature2012 = c("AJCC_Stage_nature2012", "AJCC Stage", "pathologic_stage"),
  pathologic_T = c("pathologic_T", "pathologic T"),
  pathologic_N = c("pathologic_N", "pathologic N"),
  pathologic_M = c("pathologic_M", "pathologic M"),
  histological_type = c("histological_type", "histological type"),
  menopause_status = c("menopause_status", "menopause status"),
  radiation_therapy = c("radiation_therapy", "radiation therapy"),
  history_of_neoadjuvant_treatment = c("history_of_neoadjuvant_treatment", "history of neoadjuvant treatment")
)

for (nm in names(optional_clinical)) {
  brca_clinical_modeling_data <- add_optional_column(
    brca_clinical_modeling_data,
    merged,
    nm,
    optional_clinical[[nm]],
    transform = clean_empty_as_na
  )
}

if ("AJCC_Stage_nature2012" %in% names(brca_clinical_modeling_data)) {
  brca_clinical_modeling_data$AJCC_stage_broad <- make_stage_broad(
    brca_clinical_modeling_data$AJCC_Stage_nature2012
  )
} else {
  brca_clinical_modeling_data$AJCC_stage_broad <- factor(
    rep(NA_character_, nrow(brca_clinical_modeling_data)),
    levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
  )
}

brca_clinical_modeling_data$stage_advanced <- make_stage_advanced(
  brca_clinical_modeling_data$AJCC_stage_broad
)

for (v in intersect(c(
  "ER_Status_nature2012", "PR_Status_nature2012", "HER2_Final_Status_nature2012",
  "pathologic_T", "pathologic_N", "pathologic_M", "histological_type",
  "menopause_status", "radiation_therapy", "history_of_neoadjuvant_treatment"
), names(brca_clinical_modeling_data))) {
  brca_clinical_modeling_data[[v]] <- safe_factor(brca_clinical_modeling_data[[v]], min_n = 5)
}

if ("AJCC_stage_broad" %in% names(brca_clinical_modeling_data)) {
  brca_clinical_modeling_data$AJCC_stage_broad <- safe_factor(
    brca_clinical_modeling_data$AJCC_stage_broad,
    ref = "Stage I",
    min_n = 5
  )
}

brca_clinical_modeling_data <- drop_forbidden_columns(brca_clinical_modeling_data)

clinical_only_vars <- c(
  "age_at_initial_pathologic_diagnosis",
  "age_group",
  "ER_Status_nature2012",
  "PR_Status_nature2012",
  "HER2_Final_Status_nature2012",
  "AJCC_stage_broad",
  "pathologic_T",
  "pathologic_N",
  "pathologic_M",
  "histological_type",
  "menopause_status",
  "radiation_therapy",
  "history_of_neoadjuvant_treatment"
)
clinical_only_vars <- intersect(clinical_only_vars, names(brca_clinical_modeling_data))

brca_clinical_level_summary_pfi <- inspect_levels(
  brca_clinical_modeling_data,
  clinical_only_vars,
  event_col = "PFI",
  min_n = 5,
  min_events = 2
)

eligible_cox_vars <- unique(brca_clinical_level_summary_pfi$variable[
  brca_clinical_level_summary_pfi$eligible
])
eligible_cox_vars <- intersect(eligible_cox_vars, clinical_only_vars)

brca_clinical_ref_template <- make_reference_template(
  brca_clinical_modeling_data,
  eligible_cox_vars
)
ref_levels <- as.list(brca_clinical_ref_template$reference)
names(ref_levels) <- brca_clinical_ref_template$variable
ref_levels <- ref_levels[ref_levels != "continuous" & !is.na(ref_levels)]

if (length(eligible_cox_vars) == 0L) {
  warning("No clinical variables passed prefiltering for Cox models.", call. = FALSE)
  brca_cox_univariable_clinical <- data.frame()
  brca_cox_adjusted_clinical <- data.frame()
} else {
  brca_cox_univariable_clinical <- tryCatch(
    get_cox(
      data = brca_clinical_modeling_data,
      time_col = "PFI.time",
      event_col = "PFI",
      vars = eligible_cox_vars,
      model = "univariable",
      ref_levels = ref_levels,
      min_n = 10,
      min_n_per_level = 5,
      min_events = 5,
      min_events_per_level = 2,
      verbose = TRUE
    ),
    error = function(e) {
      warning("Univariable Cox modeling skipped: ", conditionMessage(e), call. = FALSE)
      data.frame(skip_reason = conditionMessage(e), stringsAsFactors = FALSE)
    }
  )

  adjust_vars <- intersect(
    c("age_at_initial_pathologic_diagnosis", "AJCC_stage_broad"),
    eligible_cox_vars
  )
  adjusted_vars <- setdiff(eligible_cox_vars, adjust_vars)

  if (length(adjust_vars) == 0L || length(adjusted_vars) == 0L) {
    brca_cox_adjusted_clinical <- data.frame(
      skip_reason = "No eligible adjustment variables or target variables for adjusted Cox models.",
      stringsAsFactors = FALSE
    )
  } else {
    brca_cox_adjusted_clinical <- tryCatch(
      get_cox(
        data = brca_clinical_modeling_data,
        time_col = "PFI.time",
        event_col = "PFI",
        vars = adjusted_vars,
        model = "adjusted",
        adjust_vars = adjust_vars,
        keep_adjust_terms = FALSE,
        ref_levels = ref_levels,
        min_n = 10,
        min_n_per_level = 5,
        min_events = 5,
        min_events_per_level = 2,
        verbose = TRUE
      ),
      error = function(e) {
        warning("Adjusted Cox modeling skipped: ", conditionMessage(e), call. = FALSE)
        data.frame(skip_reason = conditionMessage(e), stringsAsFactors = FALSE)
      }
    )
  }
}

set.seed(2025)

glm_predictor_candidates <- c(
  "age_at_initial_pathologic_diagnosis",
  "ER_Status_nature2012",
  "PR_Status_nature2012",
  "HER2_Final_Status_nature2012",
  "histological_type",
  "menopause_status"
)
glm_predictors <- intersect(glm_predictor_candidates, names(brca_clinical_modeling_data))

glm_data <- brca_clinical_modeling_data[, c("stage_advanced", glm_predictors), drop = FALSE]
glm_data <- glm_data[stats::complete.cases(glm_data), , drop = FALSE]
glm_data <- glm_data[glm_data$stage_advanced %in% c(0, 1), , drop = FALSE]

for (v in glm_predictors) {
  if (!is.numeric(glm_data[[v]]) && !is.integer(glm_data[[v]])) {
    glm_data[[v]] <- droplevels(factor(glm_data[[v]]))
  }
}

valid_glm_predictors <- glm_predictors[vapply(glm_predictors, function(v) {
  x <- glm_data[[v]]
  if (is.numeric(x) || is.integer(x)) {
    length(unique(x)) >= 2L
  } else {
    nlevels(droplevels(factor(x))) >= 2L
  }
}, logical(1L))]

glm_data <- glm_data[, c("stage_advanced", valid_glm_predictors), drop = FALSE]

if (nrow(glm_data) < 20L ||
    length(unique(glm_data$stage_advanced)) < 2L ||
    length(valid_glm_predictors) == 0L) {
  warning("Clinical GLM/ROC skipped because data are too sparse after filtering.", call. = FALSE)
  brca_glm_stage_train <- data.frame()
  brca_glm_stage_test <- data.frame()
  brca_glm_stage_clinical <- data.frame(skip_reason = "Too sparse for clinical stage_advanced GLM.", stringsAsFactors = FALSE)
  brca_glm_stage_predictions <- data.frame()
  brca_roc_stage_clinical <- list(
    skip_reason = "Too sparse for clinical stage_advanced ROC.",
    auc_table = data.frame()
  )
} else {
  idx0 <- which(glm_data$stage_advanced == 0)
  idx1 <- which(glm_data$stage_advanced == 1)

  train_idx <- c(
    sample(idx0, size = max(1L, floor(0.60 * length(idx0)))),
    sample(idx1, size = max(1L, floor(0.60 * length(idx1))))
  )

  brca_glm_stage_train <- glm_data[sort(train_idx), , drop = FALSE]
  brca_glm_stage_test <- glm_data[setdiff(seq_len(nrow(glm_data)), sort(train_idx)), , drop = FALSE]

  if (nrow(brca_glm_stage_test) == 0L ||
      length(unique(brca_glm_stage_train$stage_advanced)) < 2L ||
      length(unique(brca_glm_stage_test$stage_advanced)) < 2L) {
    warning("Clinical GLM/ROC skipped because stratified split lacked both outcome classes.", call. = FALSE)
    brca_glm_stage_clinical <- data.frame(skip_reason = "Stratified split lacked both outcome classes.", stringsAsFactors = FALSE)
    brca_glm_stage_predictions <- data.frame()
    brca_roc_stage_clinical <- list(
      skip_reason = "Stratified split lacked both outcome classes.",
      auc_table = data.frame()
    )
  } else {
    brca_glm_stage_clinical <- get_glm(
      data = brca_glm_stage_train,
      outcome = "stage_advanced",
      predictors = valid_glm_predictors,
      family = "binomial",
      adjust_method = "BH",
      verbose = FALSE
    )

    clinical_glm_fit <- attr(brca_glm_stage_clinical, "model")
    predicted_probability <- stats::predict(
      clinical_glm_fit,
      newdata = brca_glm_stage_test,
      type = "response"
    )

    brca_glm_stage_predictions <- data.frame(
      row_id = rownames(brca_glm_stage_test),
      stage_advanced = brca_glm_stage_test$stage_advanced,
      predicted_probability = as.numeric(predicted_probability),
      predicted_class = as.integer(predicted_probability >= 0.5),
      stringsAsFactors = FALSE
    )

    brca_roc_stage_clinical <- nice_ROC(
      models = list("Clinical only" = clinical_glm_fit),
      data = brca_glm_stage_test,
      outcome = "stage_advanced",
      show_delong = FALSE,
      plot_title = "TCGA-BRCA clinical stage model",
      return_data = TRUE
    )
  }
}

require_columns(
  brca_clinical_modeling_data,
  c("sample", "sampleID", "patient", "OS", "OS.time", "PFI", "PFI.time",
    "age_at_initial_pathologic_diagnosis", "age_group", "AJCC_stage_broad",
    "stage_advanced"),
  label = "brca_clinical_modeling_data"
)

dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.table(
  brca_clinical_modeling_data,
  file = file.path(intermediate_dir, "brca_clinical_modeling_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

utils::write.table(
  brca_cox_univariable_clinical,
  file = file.path(intermediate_dir, "brca_cox_univariable_clinical.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

utils::write.table(
  brca_cox_adjusted_clinical,
  file = file.path(intermediate_dir, "brca_cox_adjusted_clinical.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

utils::write.table(
  brca_glm_stage_predictions,
  file = file.path(intermediate_dir, "brca_glm_stage_predictions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

usethis::use_data(brca_clinical_modeling_data, overwrite = TRUE)
usethis::use_data(brca_cox_univariable_clinical, overwrite = TRUE)
usethis::use_data(brca_cox_adjusted_clinical, overwrite = TRUE)
usethis::use_data(brca_glm_stage_clinical, overwrite = TRUE)
usethis::use_data(brca_roc_stage_clinical, overwrite = TRUE)
usethis::use_data(brca_clinical_level_summary_pfi, overwrite = TRUE)
usethis::use_data(brca_clinical_ref_template, overwrite = TRUE)
usethis::use_data(brca_glm_stage_train, overwrite = TRUE)
usethis::use_data(brca_glm_stage_test, overwrite = TRUE)
usethis::use_data(brca_glm_stage_predictions, overwrite = TRUE)

message("Saved clinical-only BRCA modeling objects to data/.")
message("Wrote clinical modeling TSV outputs to inst/extdata/brca/intermediate/.")
