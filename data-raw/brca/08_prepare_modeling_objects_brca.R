## 08_prepare_modeling_objects_brca.R
##
## Purpose:
##   Prepare real TCGA-BRCA clinical-only modeling objects for OmicsKit examples:
##   - nice_KM()
##   - get_cox()
##   - nice_forest()
##   - get_glm()
##   - nice_ROC()
##
## Key correction:
##   Survival endpoints are merged explicitly from BRCA_survival.txt using
##   sample-level barcodes first and patient-level fallback second. This prevents
##   PFI from becoming all NA when PFI.time is present.
##
## Inputs:
##   - data/brca_metadata.rda
##   - inst/extdata/brca/raw_xena/BRCA_survival.txt
##
## Outputs:
##   Package data objects in data/:
##   - brca_clinical_modeling_data
##   - brca_cox_univariable_clinical
##   - brca_cox_adjusted_clinical
##   - brca_glm_stage_clinical
##   - brca_roc_stage_clinical
##   - brca_clinical_level_summary_pfi
##   - brca_clinical_ref_template
##   - brca_glm_stage_train
##   - brca_glm_stage_test
##   - brca_glm_stage_predictions
##   - brca_clinical_modeling_log
##
## Intermediate TSV outputs in inst/extdata/brca/intermediate/.

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
data_dir <- project_file("data")

metadata_file <- file.path(data_dir, "brca_metadata.rda")
survival_file <- file.path(raw_dir, "BRCA_survival.txt")

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "Required input is missing for BRCA clinical modeling: ", label, "\n",
      "Expected location: ", normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }
}

require_file(metadata_file, "data/brca_metadata.rda")
require_file(survival_file, "BRCA_survival.txt")

suppressPackageStartupMessages({
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.", call. = FALSE)
  }
  if (!requireNamespace("usethis", quietly = TRUE)) {
    stop("Package 'usethis' is required.", call. = FALSE)
  }
})

has_value <- function(x) {
  y <- trimws(as.character(x))
  !is.na(y) & nzchar(y) &
    !(tolower(y) %in% c(
      "na", "n/a", "nan", "null", "none", "unknown", "not available",
      "not reported", "not applicable", "[not available]",
      "[not applicable]", "[unknown]", "--", ""
    ))
}

clean_text <- function(x) {
  y <- trimws(as.character(x))
  y[!has_value(y)] <- NA_character_
  y
}

clean_tcga_barcode <- function(x, level = c("sample", "patient")) {
  level <- match.arg(level)
  y <- toupper(trimws(as.character(x)))
  y <- gsub("\\.", "-", y)
  y[!has_value(y)] <- NA_character_

  out <- rep(NA_character_, length(y))

  if (identical(level, "patient")) {
    hit <- regexpr("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", y, perl = TRUE)
    ok <- !is.na(y) & hit > 0L
    out[ok] <- regmatches(y, hit)[ok]
    return(out)
  }

  hit <- regexpr("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-[0-9]{2}", y, perl = TRUE)
  ok <- !is.na(y) & hit > 0L
  out[ok] <- regmatches(y, hit)[ok]
  out
}

as_numeric_clean <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

safe_factor <- function(x) {
  y <- clean_text(x)
  factor(y)
}

clean_empty_as_na <- function(data) {
  data[] <- lapply(data, function(x) {
    if (is.character(x) || is.factor(x)) clean_text(x) else x
  })
  data
}

make_age_group <- function(age) {
  age_num <- as_numeric_clean(age)
  out <- ifelse(is.na(age_num), NA_character_, ifelse(age_num >= 60, "Age >= 60", "Age < 60"))
  factor(out, levels = c("Age < 60", "Age >= 60"))
}

make_stage_broad <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[!has_value(y)] <- NA_character_
  out <- rep(NA_character_, length(y))
  out[grepl("STAGE I($|[A-Z]| )", y)] <- "Stage I"
  out[grepl("STAGE II($|[A-Z]| )", y)] <- "Stage II"
  out[grepl("STAGE III($|[A-Z]| )", y)] <- "Stage III"
  out[grepl("STAGE IV($|[A-Z]| )", y)] <- "Stage IV"
  factor(out, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
}

make_stage_advanced <- function(stage_broad) {
  y <- as.character(stage_broad)
  out <- rep(NA_integer_, length(y))
  out[y %in% c("Stage I", "Stage II")] <- 0L
  out[y %in% c("Stage III", "Stage IV")] <- 1L
  out
}

collapse_rare_levels <- function(x, min_n = 10, other_label = "Other") {
  y <- clean_text(x)
  tab <- table(y, useNA = "no")
  rare <- names(tab)[tab < min_n]
  y[!is.na(y) & y %in% rare] <- other_label
  factor(y)
}

require_columns <- function(data, cols, label = "data") {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", label, ": ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

inspect_levels <- function(data, vars, event_col = NULL, include_missing = TRUE) {
  vars <- intersect(vars, names(data))
  if (length(vars) == 0L) {
    return(data.frame())
  }

  out <- lapply(vars, function(v) {
    x <- data[[v]]
    if (is.factor(x)) x <- as.character(x)
    x <- clean_text(x)

    x_label <- if (include_missing) ifelse(is.na(x), "<missing>", x) else x[!is.na(x)]
    tab <- as.data.frame(table(x_label, useNA = "no"), stringsAsFactors = FALSE)
    names(tab) <- c("level", "n")
    tab$variable <- v
    tab$percent <- round(100 * tab$n / sum(tab$n), 2)

    if (!is.null(event_col) && event_col %in% names(data)) {
      event <- as_numeric_clean(data[[event_col]])
      if (!include_missing) event <- event[!is.na(x)]
      event_count <- tapply(event == 1, x_label, sum, na.rm = TRUE)
      tab$n_events <- as.integer(event_count[tab$level])
      tab$n_events[is.na(tab$n_events)] <- 0L
    }

    tab[, c("variable", "level", "n", "percent", if ("n_events" %in% names(tab)) "n_events" else character(0)), drop = FALSE]
  })

  do.call(rbind, out)
}

make_reference_template <- function(level_summary) {
  if (!nrow(level_summary)) return(list())
  x <- subset(level_summary, level != "<missing>")
  split(x, x$variable) |>
    lapply(function(z) z$level[which.max(z$n)])
}

is_forbidden_modeling_var <- function(v) {
  grepl(
    paste(
      c(
        "PAM50", "RPPA", "methyl", "Methyl", "miRNA", "PANCAN", "CN_Clusters",
        "Integrated_Clusters", "GENOMIC_ID", "gistic", "mutation", "omics_score",
        "RNAseq_PANCAN", "DNAMethyl", "CNA"
      ),
      collapse = "|"
    ),
    v
  )
}

prefilter_cox_variable <- function(data, var, time_col, event_col,
                                   min_n = 30,
                                   min_events = 10,
                                   min_n_per_level = 10,
                                   min_events_per_level = 2) {
  if (!var %in% names(data)) {
    return(list(retained = FALSE, reason = "missing variable"))
  }

  d <- data[, c(time_col, event_col, var), drop = FALSE]
  names(d) <- c("time", "event", "x")
  d$time <- as_numeric_clean(d$time)
  d$event <- as_numeric_clean(d$event)
  d <- d[!is.na(d$time) & !is.na(d$event) & !is.na(d$x), , drop = FALSE]

  if (nrow(d) < min_n) return(list(retained = FALSE, reason = "too few total observations"))
  if (sum(d$event == 1, na.rm = TRUE) < min_events) return(list(retained = FALSE, reason = "too few total events"))

  if (is.character(d$x) || is.factor(d$x)) {
    x <- factor(as.character(d$x))
    if (nlevels(droplevels(x)) < 2L) return(list(retained = FALSE, reason = "fewer than two levels"))
    level_n <- table(x)
    if (any(level_n < min_n_per_level)) return(list(retained = FALSE, reason = "too few observations per level"))
    event_by_level <- tapply(d$event == 1, x, sum, na.rm = TRUE)
    if (any(event_by_level < min_events_per_level)) return(list(retained = FALSE, reason = "too few events per level"))
  } else {
    if (stats::var(as_numeric_clean(d$x), na.rm = TRUE) == 0) return(list(retained = FALSE, reason = "zero variance"))
  }

  list(retained = TRUE, reason = "retained")
}

fit_cox_univariable <- function(data, vars, time_col, event_col) {
  if (length(vars) == 0L) {
    return(data.frame(
      variable = character(), term = character(), HR = numeric(), conf.low = numeric(),
      conf.high = numeric(), p.value = numeric(), n = integer(), n_events = integer(),
      stringsAsFactors = FALSE
    ))
  }

  res <- lapply(vars, function(v) {
    d <- data[, c(time_col, event_col, v), drop = FALSE]
    names(d) <- c("time", "event", "x")
    d$time <- as_numeric_clean(d$time)
    d$event <- as_numeric_clean(d$event)
    d <- d[!is.na(d$time) & !is.na(d$event) & !is.na(d$x), , drop = FALSE]

    if (is.character(d$x)) d$x <- factor(d$x)
    if (is.factor(d$x)) d$x <- droplevels(d$x)

    fit <- tryCatch(
      survival::coxph(survival::Surv(time, event) ~ x, data = d),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      return(data.frame(
        variable = v, term = NA_character_, HR = NA_real_, conf.low = NA_real_,
        conf.high = NA_real_, p.value = NA_real_, n = nrow(d),
        n_events = sum(d$event == 1, na.rm = TRUE), note = fit$message,
        stringsAsFactors = FALSE
      ))
    }

    sm <- summary(fit)
    ci <- sm$conf.int
    co <- sm$coefficients
    data.frame(
      variable = v,
      term = rownames(co),
      HR = ci[, "exp(coef)"],
      conf.low = ci[, "lower .95"],
      conf.high = ci[, "upper .95"],
      p.value = co[, "Pr(>|z|)"],
      n = stats::nobs(fit),
      n_events = fit$nevent,
      note = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, res)
}

fit_cox_adjusted <- function(data, vars, time_col, event_col) {
  if (length(vars) == 0L) {
    return(data.frame(
      variable = character(), term = character(), HR = numeric(), conf.low = numeric(),
      conf.high = numeric(), p.value = numeric(), n = integer(), n_events = integer(),
      stringsAsFactors = FALSE
    ))
  }

  d <- data[, c(time_col, event_col, vars), drop = FALSE]
  d[[time_col]] <- as_numeric_clean(d[[time_col]])
  d[[event_col]] <- as_numeric_clean(d[[event_col]])

  for (v in vars) {
    if (is.character(d[[v]])) d[[v]] <- factor(d[[v]])
    if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  }

  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) < 30 || sum(d[[event_col]] == 1, na.rm = TRUE) < 10) {
    return(data.frame(
      variable = vars, term = NA_character_, HR = NA_real_, conf.low = NA_real_,
      conf.high = NA_real_, p.value = NA_real_, n = nrow(d),
      n_events = sum(d[[event_col]] == 1, na.rm = TRUE),
      note = "not enough complete observations/events for adjusted model",
      stringsAsFactors = FALSE
    ))
  }

  form <- stats::as.formula(paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", paste(vars, collapse = " + ")))
  fit <- tryCatch(survival::coxph(form, data = d), error = function(e) e)
  if (inherits(fit, "error")) {
    return(data.frame(
      variable = vars, term = NA_character_, HR = NA_real_, conf.low = NA_real_,
      conf.high = NA_real_, p.value = NA_real_, n = nrow(d),
      n_events = sum(d[[event_col]] == 1, na.rm = TRUE), note = fit$message,
      stringsAsFactors = FALSE
    ))
  }

  sm <- summary(fit)
  ci <- sm$conf.int
  co <- sm$coefficients
  data.frame(
    variable = sub("^`?([^`]+)`?.*$", "\\1", rownames(co)),
    term = rownames(co),
    HR = ci[, "exp(coef)"],
    conf.low = ci[, "lower .95"],
    conf.high = ci[, "upper .95"],
    p.value = co[, "Pr(>|z|)"],
    n = stats::nobs(fit),
    n_events = fit$nevent,
    note = NA_character_,
    stringsAsFactors = FALSE
  )
}

stratified_train_test_split <- function(y, train_prop = 0.6, seed = 2025) {
  set.seed(seed)
  idx0 <- which(y == 0)
  idx1 <- which(y == 1)
  train0 <- sample(idx0, size = floor(length(idx0) * train_prop))
  train1 <- sample(idx1, size = floor(length(idx1) * train_prop))
  sort(c(train0, train1))
}

align_factor_levels <- function(train, test, vars) {
  for (v in vars) {
    if (!v %in% names(train) || !v %in% names(test)) next
    if (is.factor(train[[v]]) || is.character(train[[v]])) {
      train[[v]] <- droplevels(factor(as.character(train[[v]])))
      test_raw <- as.character(test[[v]])
      unseen <- !is.na(test_raw) & !(test_raw %in% levels(train[[v]]))
      if (any(unseen) && "Other" %in% levels(train[[v]])) {
        test_raw[unseen] <- "Other"
      } else if (any(unseen)) {
        test_raw[unseen] <- NA_character_
      }
      test[[v]] <- factor(test_raw, levels = levels(train[[v]]))
    }
  }
  list(train = train, test = test)
}

has_perfect_prediction <- function(data, predictor, outcome) {
  if (!predictor %in% names(data)) return(FALSE)
  x <- data[[predictor]]
  y <- data[[outcome]]
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  if (length(unique(y)) < 2L) return(TRUE)
  if (is.factor(x) || is.character(x)) {
    tab <- table(x, y)
    any(rowSums(tab > 0) == 1)
  } else {
    FALSE
  }
}

# ── Load metadata and survival ------------------------------------------------
load(metadata_file)
if (!exists("brca_metadata")) {
  stop("data/brca_metadata.rda does not contain object brca_metadata.", call. = FALSE)
}

survival_raw <- utils::read.delim(
  survival_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE
)

require_columns(survival_raw, c("sample", "_PATIENT", "OS", "OS.time", "PFI", "PFI.time"), "BRCA_survival.txt")

survival_raw$survival_sample16 <- clean_tcga_barcode(survival_raw$sample, level = "sample")
survival_raw$survival_patient <- clean_tcga_barcode(survival_raw$`_PATIENT`, level = "patient")

for (v in c("OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time")) {
  if (v %in% names(survival_raw)) survival_raw[[v]] <- as_numeric_clean(survival_raw[[v]])
}

# Prefer exact sample-level match, then fallback to patient.
sample_match <- match(brca_metadata$sample16, survival_raw$survival_sample16)
patient_match <- match(brca_metadata$patient, survival_raw$survival_patient)
matched_index <- sample_match
fallback <- is.na(matched_index) & !is.na(patient_match)
matched_index[fallback] <- patient_match[fallback]

message("Survival merge diagnostics:")
message("  Clinical rows: ", nrow(brca_metadata))
message("  Matched by sample16: ", sum(!is.na(sample_match)))
message("  Matched by patient fallback: ", sum(fallback))
message("  Unmatched: ", sum(is.na(matched_index)))

# ── Build clinical-only modeling dataset -------------------------------------
clinical_vars <- c(
  "sampleID", "sample16", "patient", "sample_code", "sample_type", "tumor_normal",
  "ER_status_raw", "ER_group", "PR_status_raw", "HER2_status_raw",
  "age_at_initial_pathologic_diagnosis",
  "ER_Status_nature2012", "PR_Status_nature2012", "HER2_Final_Status_nature2012",
  "AJCC_Stage_nature2012", "pathologic_T", "pathologic_N", "pathologic_M",
  "histological_type", "menopause_status", "radiation_therapy",
  "history_of_neoadjuvant_treatment"
)
clinical_vars <- intersect(clinical_vars, names(brca_metadata))

brca_clinical_modeling_data <- brca_metadata[, clinical_vars, drop = FALSE]
brca_clinical_modeling_data <- clean_empty_as_na(brca_clinical_modeling_data)

# Explicitly assign survival endpoints from matched survival rows.
for (v in c("OS", "OS.time", "DSS", "DSS.time", "DFI", "DFI.time", "PFI", "PFI.time")) {
  if (v %in% names(survival_raw)) {
    brca_clinical_modeling_data[[v]] <- survival_raw[[v]][matched_index]
  }
}

brca_clinical_modeling_data$age_at_initial_pathologic_diagnosis <- as_numeric_clean(
  brca_clinical_modeling_data$age_at_initial_pathologic_diagnosis
)
brca_clinical_modeling_data$age_group <- make_age_group(
  brca_clinical_modeling_data$age_at_initial_pathologic_diagnosis
)

if ("AJCC_Stage_nature2012" %in% names(brca_clinical_modeling_data)) {
  brca_clinical_modeling_data$AJCC_stage_broad <- make_stage_broad(
    brca_clinical_modeling_data$AJCC_Stage_nature2012
  )
} else {
  brca_clinical_modeling_data$AJCC_stage_broad <- factor(NA_character_)
}
brca_clinical_modeling_data$stage_advanced <- make_stage_advanced(
  brca_clinical_modeling_data$AJCC_stage_broad
)

categorical_vars <- intersect(
  c(
    "ER_Status_nature2012", "PR_Status_nature2012", "HER2_Final_Status_nature2012",
    "histological_type", "menopause_status", "radiation_therapy",
    "history_of_neoadjuvant_treatment", "pathologic_T", "pathologic_N", "pathologic_M",
    "AJCC_stage_broad"
  ),
  names(brca_clinical_modeling_data)
)

for (v in categorical_vars) {
  brca_clinical_modeling_data[[v]] <- collapse_rare_levels(brca_clinical_modeling_data[[v]], min_n = 10)
}

message("Survival endpoint diagnostics after merge:")
message("PFI distribution:")
print(table(brca_clinical_modeling_data$PFI, useNA = "ifany"))
message("OS distribution:")
print(table(brca_clinical_modeling_data$OS, useNA = "ifany"))
message("PFI.time summary:")
print(summary(brca_clinical_modeling_data$PFI.time))
message("OS.time summary:")
print(summary(brca_clinical_modeling_data$OS.time))

# ── Select Cox endpoint -------------------------------------------------------
cox_endpoint <- list(time_col = "OS.time", event_col = "OS", endpoint = "OS")

pfi_events <- if ("PFI" %in% names(brca_clinical_modeling_data)) {
  sum(brca_clinical_modeling_data$PFI == 1, na.rm = TRUE)
} else 0L
pfi_times <- if ("PFI.time" %in% names(brca_clinical_modeling_data)) {
  sum(!is.na(brca_clinical_modeling_data$PFI.time))
} else 0L

if (pfi_events >= 20 && pfi_times >= 100) {
  cox_endpoint <- list(time_col = "PFI.time", event_col = "PFI", endpoint = "PFI")
}

message("Selected Cox endpoint: ", cox_endpoint$endpoint)
message("  time_col: ", cox_endpoint$time_col)
message("  event_col: ", cox_endpoint$event_col)
message("  events: ", sum(brca_clinical_modeling_data[[cox_endpoint$event_col]] == 1, na.rm = TRUE))

# ── Cox clinical-only ---------------------------------------------------------
cox_candidates <- c(
  "age_at_initial_pathologic_diagnosis", "age_group",
  "ER_Status_nature2012", "PR_Status_nature2012", "HER2_Final_Status_nature2012",
  "AJCC_stage_broad", "histological_type", "menopause_status", "radiation_therapy",
  "history_of_neoadjuvant_treatment", "pathologic_T", "pathologic_N", "pathologic_M"
)
cox_candidates <- intersect(cox_candidates, names(brca_clinical_modeling_data))
forbidden <- cox_candidates[is_forbidden_modeling_var(cox_candidates)]
if (length(forbidden) > 0L) {
  warning("Dropping forbidden non-clinical or leakage variables: ", paste(forbidden, collapse = ", "), call. = FALSE)
}
cox_candidates <- setdiff(cox_candidates, forbidden)

cox_prefilter <- do.call(rbind, lapply(cox_candidates, function(v) {
  chk <- prefilter_cox_variable(
    brca_clinical_modeling_data,
    v,
    time_col = cox_endpoint$time_col,
    event_col = cox_endpoint$event_col,
    min_n = 30,
    min_events = 10,
    min_n_per_level = 10,
    min_events_per_level = 2
  )
  data.frame(variable = v, retained = chk$retained, reason = chk$reason, stringsAsFactors = FALSE)
}))

cox_retained <- cox_prefilter$variable[cox_prefilter$retained]
message("Cox candidates retained: ", paste(cox_retained, collapse = ", "))
message("Cox candidates dropped and reason:")
print(cox_prefilter[!cox_prefilter$retained, , drop = FALSE])

brca_cox_univariable_clinical <- fit_cox_univariable(
  brca_clinical_modeling_data,
  vars = cox_retained,
  time_col = cox_endpoint$time_col,
  event_col = cox_endpoint$event_col
)

adjusted_vars <- intersect(c("age_at_initial_pathologic_diagnosis", "AJCC_stage_broad", "ER_Status_nature2012"), cox_retained)
brca_cox_adjusted_clinical <- fit_cox_adjusted(
  brca_clinical_modeling_data,
  vars = adjusted_vars,
  time_col = cox_endpoint$time_col,
  event_col = cox_endpoint$event_col
)

brca_clinical_level_summary_pfi <- inspect_levels(
  brca_clinical_modeling_data,
  vars = cox_candidates,
  event_col = cox_endpoint$event_col,
  include_missing = TRUE
)
brca_clinical_ref_template <- make_reference_template(brca_clinical_level_summary_pfi)

# ── GLM / ROC clinical-only ---------------------------------------------------
glm_predictors <- intersect(
  c(
    "age_at_initial_pathologic_diagnosis",
    "ER_Status_nature2012", "PR_Status_nature2012", "HER2_Final_Status_nature2012",
    "histological_type", "menopause_status", "radiation_therapy",
    "history_of_neoadjuvant_treatment"
  ),
  names(brca_clinical_modeling_data)
)

glm_predictors <- setdiff(glm_predictors, glm_predictors[is_forbidden_modeling_var(glm_predictors)])

# Never use variables that define/encode stage_advanced.
glm_predictors <- setdiff(glm_predictors, c(
  "AJCC_Stage_nature2012", "AJCC_stage_broad", "pathologic_T", "pathologic_N", "pathologic_M"
))

glm_data <- brca_clinical_modeling_data[, c("stage_advanced", glm_predictors), drop = FALSE]
glm_data <- glm_data[!is.na(glm_data$stage_advanced), , drop = FALSE]

for (v in glm_predictors) {
  if (is.character(glm_data[[v]]) || is.factor(glm_data[[v]])) {
    glm_data[[v]] <- collapse_rare_levels(glm_data[[v]], min_n = 10)
  }
}

glm_data <- glm_data[stats::complete.cases(glm_data), , drop = FALSE]

if (length(unique(glm_data$stage_advanced)) < 2L || nrow(glm_data) < 50) {
  warning("Not enough complete data to fit stage_advanced GLM. Saving empty GLM objects.", call. = FALSE)
  brca_glm_stage_train <- glm_data[0, , drop = FALSE]
  brca_glm_stage_test <- glm_data[0, , drop = FALSE]
  brca_glm_stage_clinical <- NULL
  brca_glm_stage_predictions <- data.frame()
  brca_roc_stage_clinical <- list(auc = NA_real_, roc = NULL, note = "not enough data")
  glm_retained <- character(0)
  glm_dropped <- data.frame(variable = glm_predictors, reason = "not enough complete data", stringsAsFactors = FALSE)
} else {
  train_idx <- stratified_train_test_split(glm_data$stage_advanced, train_prop = 0.6, seed = 2025)
  train <- glm_data[train_idx, , drop = FALSE]
  test <- glm_data[-train_idx, , drop = FALSE]

  aligned <- align_factor_levels(train, test, glm_predictors)
  train <- aligned$train
  test <- aligned$test

  glm_dropped <- data.frame(variable = character(), reason = character(), stringsAsFactors = FALSE)
  glm_retained <- glm_predictors

  # Drop predictors with fewer than two levels or clear perfect prediction in train.
  for (v in glm_predictors) {
    if (!v %in% names(train)) next
    if (is.factor(train[[v]]) || is.character(train[[v]])) {
      if (nlevels(droplevels(factor(train[[v]]))) < 2L) {
        glm_retained <- setdiff(glm_retained, v)
        glm_dropped <- rbind(glm_dropped, data.frame(variable = v, reason = "fewer than two levels in train", stringsAsFactors = FALSE))
      } else if (has_perfect_prediction(train, v, "stage_advanced")) {
        glm_retained <- setdiff(glm_retained, v)
        glm_dropped <- rbind(glm_dropped, data.frame(variable = v, reason = "perfect prediction in train", stringsAsFactors = FALSE))
      }
    } else {
      if (stats::var(train[[v]], na.rm = TRUE) == 0) {
        glm_retained <- setdiff(glm_retained, v)
        glm_dropped <- rbind(glm_dropped, data.frame(variable = v, reason = "zero variance", stringsAsFactors = FALSE))
      }
    }
  }

  train <- train[, c("stage_advanced", glm_retained), drop = FALSE]
  test <- test[, c("stage_advanced", glm_retained), drop = FALSE]
  train <- train[stats::complete.cases(train), , drop = FALSE]
  test <- test[stats::complete.cases(test), , drop = FALSE]

  message("GLM train outcome distribution:")
  print(table(train$stage_advanced, useNA = "ifany"))
  message("GLM test outcome distribution:")
  print(table(test$stage_advanced, useNA = "ifany"))
  message("GLM predictors retained: ", paste(glm_retained, collapse = ", "))
  message("GLM predictors dropped and reason:")
  print(glm_dropped)

  if (length(glm_retained) == 0L) {
    warning("No GLM predictors retained. Saving empty GLM objects.", call. = FALSE)
    brca_glm_stage_clinical <- NULL
    brca_glm_stage_predictions <- data.frame()
    brca_roc_stage_clinical <- list(auc = NA_real_, roc = NULL, note = "no predictors retained")
  } else {
    glm_formula <- stats::as.formula(paste("stage_advanced ~", paste(glm_retained, collapse = " + ")))
    separation_warning <- FALSE
    brca_glm_stage_clinical <- withCallingHandlers(
      stats::glm(glm_formula, data = train, family = stats::binomial()),
      warning = function(w) {
        if (grepl("fitted probabilities numerically 0 or 1 occurred", conditionMessage(w))) {
          separation_warning <<- TRUE
          invokeRestart("muffleWarning")
        }
      }
    )

    pred <- stats::predict(brca_glm_stage_clinical, newdata = test, type = "response")
    brca_glm_stage_predictions <- data.frame(
      observed = test$stage_advanced,
      predicted_probability = as.numeric(pred),
      stringsAsFactors = FALSE
    )

    if (requireNamespace("pROC", quietly = TRUE) && length(unique(test$stage_advanced)) == 2L) {
      roc_obj <- pROC::roc(
        response = brca_glm_stage_predictions$observed,
        predictor = brca_glm_stage_predictions$predicted_probability,
        quiet = TRUE
      )
      brca_roc_stage_clinical <- list(
        roc = roc_obj,
        auc = as.numeric(pROC::auc(roc_obj)),
        outcome = "stage_advanced",
        predictors = glm_retained,
        note = if (separation_warning) "glm separation warning occurred" else NA_character_
      )
    } else {
      brca_roc_stage_clinical <- list(
        roc = NULL,
        auc = NA_real_,
        outcome = "stage_advanced",
        predictors = glm_retained,
        note = "pROC unavailable or test outcome has fewer than two classes"
      )
    }
  }

  brca_glm_stage_train <- train
  brca_glm_stage_test <- test
}

# ── Log -----------------------------------------------------------------------
brca_clinical_modeling_log <- list(
  created = as.character(Sys.time()),
  survival_merge = list(
    n_clinical_rows = nrow(brca_metadata),
    n_matched_by_sample16 = sum(!is.na(sample_match)),
    n_matched_by_patient_fallback = sum(fallback),
    n_unmatched = sum(is.na(matched_index))
  ),
  cox_endpoint = cox_endpoint,
  cox_prefilter = cox_prefilter,
  cox_retained = cox_retained,
  glm_predictors_retained = glm_retained,
  glm_predictors_dropped = glm_dropped,
  note = "Clinical-only modeling. Omics/PANCAN cluster variables are excluded."
)

# ── Save outputs --------------------------------------------------------------
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.table(
  brca_clinical_modeling_data,
  file = file.path(intermediate_dir, "brca_clinical_modeling_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
utils::write.table(
  brca_cox_univariable_clinical,
  file = file.path(intermediate_dir, "brca_cox_univariable_clinical.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
utils::write.table(
  brca_cox_adjusted_clinical,
  file = file.path(intermediate_dir, "brca_cox_adjusted_clinical.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
utils::write.table(
  brca_glm_stage_predictions,
  file = file.path(intermediate_dir, "brca_glm_stage_predictions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

usethis::use_data(brca_clinical_modeling_data, compress = "xz", overwrite = TRUE)
usethis::use_data(brca_cox_univariable_clinical, compress = "xz", overwrite = TRUE)
usethis::use_data(brca_cox_adjusted_clinical, compress = "xz", overwrite = TRUE)
usethis::use_data(brca_glm_stage_clinical, compress = "xz", overwrite = TRUE)
usethis::use_data(brca_roc_stage_clinical, compress = "xz", overwrite = TRUE)
usethis::use_data(brca_clinical_level_summary_pfi, compress = "xz", overwrite = TRUE)
usethis::use_data(brca_clinical_ref_template, compress = "xz", overwrite = TRUE)
usethis::use_data(brca_glm_stage_train, compress = "xz", overwrite = TRUE)
usethis::use_data(brca_glm_stage_test, compress = "xz", overwrite = TRUE)
usethis::use_data(brca_glm_stage_predictions, compress = "xz", overwrite = TRUE)
usethis::use_data(brca_clinical_modeling_log, compress = "xz", overwrite = TRUE)

message("Saved clinical-only BRCA modeling objects to data/.")
message("Wrote clinical modeling TSV outputs to inst/extdata/brca/intermediate/.")
