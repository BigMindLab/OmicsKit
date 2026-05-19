## 01_prepare_metadata_brca.R
##
## Purpose:
##   Prepare a harmonized TCGA-BRCA metadata object for DEA, proteomics,
##   pathway, omics-layer, concordance, and clinical modeling scripts.
##
## Inputs:
##   - inst/extdata/brca/raw_xena/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix
##   - inst/extdata/brca/raw_xena/BRCA_survival.txt
##
## Outputs:
##   - data/brca_metadata.rda
##   - inst/extdata/brca/intermediate/brca_metadata.tsv
##
## Expected file locations:
##   - Raw files: inst/extdata/brca/raw_xena/
##   - Intermediate TSV: inst/extdata/brca/intermediate/
##   - Package data object: data/

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

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "Required raw file is missing for BRCA metadata: ", label, "\n",
      "Expected location: ", normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }
}

require_file(clinical_file, "TCGA.BRCA.sampleMap_BRCA_clinicalMatrix")
require_file(survival_file, "BRCA_survival.txt")

has_value <- function(x) {
  y <- trimws(as.character(x))
  !is.na(y) &
    nzchar(y) &
    !(tolower(y) %in% c(
      "na", "n/a", "nan", "null", "none", "unknown", "not available",
      "not reported", "not applicable", "[not available]",
      "[not applicable]", "[unknown]", "--"
    ))
}

clean_tcga_barcode <- function(x, level = c("sample", "patient")) {
  level <- match.arg(level)
  y <- toupper(trimws(as.character(x)))
  y <- gsub("\\.", "-", y)
  y[!has_value(y)] <- NA_character_

  n_chars <- if (identical(level, "sample")) 16L else 12L
  ifelse(!is.na(y) & nchar(y) >= n_chars, substr(y, 1L, n_chars), NA_character_)
}

tcga_sample_code <- function(x) {
  y <- clean_tcga_barcode(x, level = "sample")
  ifelse(!is.na(y) & nchar(y) >= 15L, substr(y, 14L, 15L), NA_character_)
}

make_er_group <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))

  valid <- has_value(y)
  positive <- valid & grepl("positive|\\bpos\\b|\\+", y)
  negative <- valid & grepl("negative|\\bneg\\b", y)

  out[positive] <- "ER_positive"
  out[negative & !positive] <- "ER_negative"
  out
}

normalize_name <- function(x) {
  tolower(gsub("[^a-z0-9]+", "", x))
}

find_column <- function(data, candidates) {
  data_names <- names(data)
  data_norm <- normalize_name(data_names)
  candidate_norm <- normalize_name(candidates)

  matched <- match(candidate_norm, data_norm)
  if (any(!is.na(matched))) {
    return(data_names[matched[which(!is.na(matched))[1L]]])
  }

  NULL
}

get_required_column <- function(data, candidates, label, source_file) {
  column <- find_column(data, candidates)
  if (is.null(column)) {
    stop(
      "Required column is missing for brca_metadata: ", label, "\n",
      "Expected one of: ", paste(candidates, collapse = ", "), "\n",
      "Source file: ", normalizePath(source_file, mustWork = FALSE),
      call. = FALSE
    )
  }
  column
}

get_optional_column <- function(data, candidates, label) {
  column <- find_column(data, candidates)
  if (is.null(column)) {
    warning("Optional clinical column is missing and will be omitted: ", label,
            call. = FALSE)
  }
  column
}

clean_text <- function(x) {
  y <- trimws(as.character(x))
  y[!has_value(y)] <- NA_character_
  y
}

tcga_sample_type <- function(sample_code) {
  labels <- c(
    "01" = "Primary Tumor",
    "02" = "Recurrent Solid Tumor",
    "03" = "Primary Blood Derived Cancer - Peripheral Blood",
    "05" = "Additional New Primary",
    "06" = "Metastatic",
    "07" = "Additional Metastatic",
    "10" = "Blood Derived Normal",
    "11" = "Solid Tissue Normal",
    "12" = "Buccal Cell Normal",
    "13" = "EBV Immortalized Normal",
    "14" = "Bone Marrow Normal",
    "20" = "Control Analyte",
    "40" = "Recurrent Blood Derived Cancer - Bone Marrow",
    "50" = "Cell Lines",
    "60" = "Primary Xenograft Tissue",
    "61" = "Cell Line Derived Xenograft Tissue"
  )

  unname(labels[sample_code])
}

tcga_tumor_normal <- function(sample_code) {
  out <- rep(NA_character_, length(sample_code))
  tumor_codes <- sprintf("%02d", 1:9)
  normal_codes <- sprintf("%02d", 10:19)

  out[sample_code %in% tumor_codes] <- "Tumor"
  out[sample_code %in% normal_codes] <- "Normal"
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

clinical <- read_xena_table(clinical_file)
survival <- read_xena_table(survival_file)

clinical_id_col <- get_required_column(
  clinical,
  candidates = c("sampleID", "sample", "Samples", "Sample", "ID"),
  label = "clinical sample identifier",
  source_file = clinical_file
)

survival_id_col <- get_required_column(
  survival,
  candidates = c("sample", "sampleID", "Samples", "Sample", "patient", "Patient", "_PATIENT"),
  label = "survival sample or patient identifier",
  source_file = survival_file
)

er_col <- get_required_column(
  clinical,
  candidates = c(
    "ER_Status_nature2012",
    "ER_status",
    "ER Status",
    "breast_carcinoma_estrogen_receptor_status",
    "estrogen_receptor_status"
  ),
  label = "ER status",
  source_file = clinical_file
)

pr_col <- get_required_column(
  clinical,
  candidates = c(
    "PR_Status_nature2012",
    "PR_status",
    "PR Status",
    "breast_carcinoma_progesterone_receptor_status",
    "progesterone_receptor_status"
  ),
  label = "PR status",
  source_file = clinical_file
)

her2_col <- get_required_column(
  clinical,
  candidates = c(
    "HER2_Final_Status_nature2012",
    "HER2_Status_nature2012",
    "HER2_status",
    "HER2 Status",
    "lab_proc_her2_neu_immunohistochemistry_receptor_status",
    "her2_neu_status"
  ),
  label = "HER2 status",
  source_file = clinical_file
)

age_col <- get_required_column(
  clinical,
  candidates = c(
    "age_at_initial_pathologic_diagnosis",
    "Age_at_Initial_Pathologic_Diagnosis_nature2012",
    "age_at_diagnosis",
    "age"
  ),
  label = "age_at_initial_pathologic_diagnosis",
  source_file = clinical_file
)

rna_genomic_candidates <- c("RNA_genomic_id", "RNA genomic id", "RNAseq_genomic_id")
rppa_genomic_candidates <- c("RPPA_genomic_id", "RPPA genomic id")

rna_genomic_col <- find_column(clinical, rna_genomic_candidates)
rppa_genomic_col <- find_column(clinical, rppa_genomic_candidates)

os_col <- get_required_column(
  survival,
  candidates = c("OS"),
  label = "OS",
  source_file = survival_file
)

os_time_col <- get_required_column(
  survival,
  candidates = c("OS.time", "OS_time", "OS Time"),
  label = "OS.time",
  source_file = survival_file
)

pfi_col <- get_required_column(
  survival,
  candidates = c("PFI"),
  label = "PFI",
  source_file = survival_file
)

pfi_time_col <- get_required_column(
  survival,
  candidates = c("PFI.time", "PFI_time", "PFI Time"),
  label = "PFI.time",
  source_file = survival_file
)

sample16 <- clean_tcga_barcode(clinical[[clinical_id_col]], level = "sample")

if (!is.null(rna_genomic_col)) {
  rna_sample16 <- clean_tcga_barcode(clinical[[rna_genomic_col]], level = "sample")
  sample16 <- ifelse(is.na(sample16), rna_sample16, sample16)
}

if (!is.null(rppa_genomic_col)) {
  rppa_sample16 <- clean_tcga_barcode(clinical[[rppa_genomic_col]], level = "sample")
  sample16 <- ifelse(is.na(sample16), rppa_sample16, sample16)
}

patient <- clean_tcga_barcode(clinical[[clinical_id_col]], level = "patient")
patient <- ifelse(is.na(patient), clean_tcga_barcode(sample16, level = "patient"), patient)
sample_code <- tcga_sample_code(sample16)

brca_metadata <- data.frame(
  sampleID = sample16,
  sample16 = sample16,
  patient = patient,
  sample_code = sample_code,
  sample_type = tcga_sample_type(sample_code),
  tumor_normal = tcga_tumor_normal(sample_code),
  ER_status_raw = clean_text(clinical[[er_col]]),
  ER_group = make_er_group(clinical[[er_col]]),
  PR_status_raw = clean_text(clinical[[pr_col]]),
  HER2_status_raw = clean_text(clinical[[her2_col]]),
  age_at_initial_pathologic_diagnosis = clean_text(clinical[[age_col]]),
  stringsAsFactors = FALSE
)

survival$patient <- clean_tcga_barcode(survival[[survival_id_col]], level = "patient")
survival <- survival[!is.na(survival$patient), , drop = FALSE]
survival <- survival[!duplicated(survival$patient), , drop = FALSE]

survival_match <- match(brca_metadata$patient, survival$patient)
brca_metadata$OS <- survival[[os_col]][survival_match]
brca_metadata$OS.time <- survival[[os_time_col]][survival_match]
brca_metadata$PFI <- survival[[pfi_col]][survival_match]
brca_metadata$PFI.time <- survival[[pfi_time_col]][survival_match]

optional_columns <- list(
  AJCC_Stage_nature2012 = c("AJCC_Stage_nature2012", "AJCC Stage", "pathologic_stage"),
  pathologic_T = c("pathologic_T", "pathologic T"),
  pathologic_N = c("pathologic_N", "pathologic N"),
  pathologic_M = c("pathologic_M", "pathologic M"),
  histological_type = c("histological_type", "histological type"),
  menopause_status = c("menopause_status", "menopause status"),
  radiation_therapy = c("radiation_therapy", "radiation therapy"),
  history_of_neoadjuvant_treatment = c(
    "history_of_neoadjuvant_treatment",
    "history of neoadjuvant treatment"
  ),
  RNA_genomic_id = rna_genomic_candidates,
  RPPA_genomic_id = rppa_genomic_candidates
)

for (output_name in names(optional_columns)) {
  source_col <- get_optional_column(clinical, optional_columns[[output_name]], output_name)
  if (!is.null(source_col)) {
    brca_metadata[[output_name]] <- clean_text(clinical[[source_col]])
  }
}

if (any(is.na(brca_metadata$sampleID))) {
  stop(
    "Some clinical rows do not contain a valid TCGA sample barcode in column: ",
    clinical_id_col,
    call. = FALSE
  )
}

if (any(is.na(brca_metadata$patient))) {
  stop("Some clinical rows do not contain a valid TCGA patient barcode.", call. = FALSE)
}

dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.table(
  brca_metadata,
  file = file.path(intermediate_dir, "brca_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

if (!requireNamespace("usethis", quietly = TRUE)) {
  stop(
    "Package \"usethis\" is required to save brca_metadata with usethis::use_data().",
    call. = FALSE
  )
}

usethis::use_data(brca_metadata, overwrite = TRUE)

message("Saved brca_metadata to data/brca_metadata.rda")
message("Wrote TSV copy to inst/extdata/brca/intermediate/brca_metadata.tsv")
