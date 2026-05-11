# =============================================================================
# OmicsKit — Ensembl annotation utilities
# Functions: list_ensembl_versions, list_ensembl_species, get_annotations
# Internal helpers: .resolve_ensembl_host, .get_symbol_attribute,
#                   .resolve_symbol_column, .handle_biomart_connection_error,
#                   .handle_biomart_query_error
# =============================================================================


############################
# Function get_annotations #
############################

#' Get gene or transcript annotations from Ensembl BioMart.
#'
#' This function annotates a column of transcripts or gene IDs (ENSEMBL) with information of the Biomart.
#'
#' The Gene information added include:
#' - Gene ENSEMBL ID, gene Symbol, Description, Biotype and Chromosome.
#' - Gene start, end and length
#'
#' Annotates a vector of Ensembl gene or transcript IDs using BioMart. If
#' transcript IDs are provided, they are also annotated with information
#' of the genes to which they belong. Works for any species available in
#' Ensembl — use [list_ensembl_species()] to find the correct `species` value
#' for your organism, and [list_ensembl_versions()] to verify version availability.
#'
#' @details
#' **Workflow:**
#' \enumerate{
#'   \item `list_ensembl_species()` → find dataset name.
#'   \item `list_ensembl_versions()` → confirm version exists.
#'   \item `get_annotations(ensembl_ids, species = "drerio_gene_ensembl")`.
#' }
#'
#' **Symbol resolution by species:**
#' \itemize{
#'   \item Human (`hsapiens`): uses `hgnc_symbol`, falling back to
#'         `external_gene_name` when HGNC symbol is absent.
#'   \item Mouse (`mmusculus`): uses `mgi_symbol`.
#'   \item Rat (`rnorvegicus`): uses `rgd_symbol`.
#'   \item Zebrafish (`drerio`): uses `zfin_id_symbol`.
#'   \item Drosophila (`dmelanogaster`): uses `flybasename_gene`.
#'   \item All other species: uses `external_gene_name` (universal fallback).
#' }
#'
#' @param ensembl_ids A character vector of Ensembl gene or transcript IDs as inputs.
#' @param species The BioMart dataset identifier. Default =
#'   `"hsapiens_gene_ensembl"` (human). Use [list_ensembl_species()] to find
#'   the identifier for other organisms.
#' @param mode Either `"genes"` or `"transcripts"`. Default = `"genes"`.
#' @param version Ensembl release version as a string (e.g. `"112"`, `"114"`),
#'   or `"current"`. Use [list_ensembl_versions()] to see available versions,
#'   and [list_ensembl_species()] to confirm your species exists in that
#'   version. Default = `"current"`.
#' @param filename Name of the output file without extension. If `NULL`, the
#'   result is returned but **not** saved to disk. Default =
#'   `"gene_annotations"`.
#' @param format Output format: `"csv"` or `"xlsx"`. Ignored when
#'   `filename = NULL`. Default = `"csv"`.
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @return A data frame with one row per input ID and columns: `geneID`,
#'   `symbol`, `biotype`, `chromosome`, `gene_start`, `gene_end`,
#'   `gene_length`, `description`. For `mode = "transcripts"`, a `transcriptID`
#'   column is prepended. If `filename` is not `NULL`, the data frame is also
#'   written to disk.
#'
#' @examples
#' \dontrun{
#' # Step 1 — find your species dataset
#' list_ensembl_species(version = "112", filter = "zebrafish")
#' # -> "drerio_gene_ensembl"
#'
#' # Step 2 — annotate (human, default)
#' annotations <- get_annotations(
#'   ensembl_ids = c("ENSG00000141510", "ENSG00000012048"),
#'   species     = "hsapiens_gene_ensembl",
#'   version     = "112",
#'   filename    = NULL
#' )
#'
#' # Mouse
#' annotations_mouse <- get_annotations(
#'   ensembl_ids = c("ENSMUSG00000059552", "ENSMUSG00000024610"),
#'   species     = "mmusculus_gene_ensembl",
#'   version     = "112",
#'   filename    = NULL
#' )
#'
#' # Zebrafish
#' annotations_zf <- get_annotations(
#'   ensembl_ids = c("ENSDARG00000002333"),
#'   species     = "drerio_gene_ensembl",
#'   version     = "112",
#'   filename    = NULL
#' )
#' }
#'
#' @note Requires an active internet connection to query the Ensembl BioMart.
#'   `gene_length` is computed as `gene_end - gene_start + 1` (genomic span).
#'   For TPM calculation with [tpm()], transcript-level lengths are more
#'   accurate.
#'
#' @seealso [list_ensembl_versions()], [list_ensembl_species()] to find
#'   the identifier for other organisms, [add_annotations()], [tpm()]
#'
#' @export

get_annotations <- function(ensembl_ids,
                            species  = "hsapiens_gene_ensembl",
                            mode     = "genes",
                            version  = "current",
                            filename = "gene_annotations",
                            format   = "csv") {

  # --- Dependency check -------------------------------------------------------
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "Package \"biomaRt\" must be installed to use this function.\n",
      "Install it with: BiocManager::install(\"biomaRt\")",
      call. = FALSE
    )
  }

  # --- Resolve host and connect -----------------------------------------------
  host <- .resolve_ensembl_host(version)

  message("Connecting to Ensembl BioMart  [version: ", version,
          " | species: ", species, "] ...")

  ensembl <- tryCatch(
    biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                     dataset = species,
                     host    = host),
    error = function(e) .handle_biomart_connection_error(e, host, species, version)
  )

  # --- Detect appropriate symbol attribute for this species -------------------
  symbol_info <- .get_symbol_attribute(ensembl, species)
  is_human    <- grepl("^hsapiens", species)

  # --- Define BioMart attributes ----------------------------------------------
  if (is_human) {
    # Human: fetch both hgnc_symbol + external_gene_name, merge later
    annotations_attrs <- c(
      "ensembl_gene_id", "hgnc_symbol",    "external_gene_name",
      "gene_biotype",    "chromosome_name",
      "start_position",  "end_position",   "description"
    )
  } else {
    annotations_attrs <- c(
      "ensembl_gene_id", symbol_info$attr,
      "gene_biotype",    "chromosome_name",
      "start_position",  "end_position",   "description"
    )
  }

  final_col_names <- c(
    "geneID", "symbol", "biotype", "chromosome",
    "gene_start", "gene_end", "description"
  )

  # --- Query BioMart ----------------------------------------------------------
  if (mode == "transcripts") {

    df <- data.frame(transcriptID = ensembl_ids)

    genemap <- tryCatch(
      biomaRt::getBM(
        attributes = c("ensembl_transcript_id_version", annotations_attrs),
        filters    = "ensembl_transcript_id_version",
        values     = df$transcriptID,
        mart       = ensembl
      ),
      error = function(e) .handle_biomart_query_error(e, host)
    )

    genemap <- .resolve_symbol_column(genemap, is_human, symbol_info$attr)

    keep_cols <- c("ensembl_transcript_id_version", "ensembl_gene_id",
                   "symbol", "gene_biotype", "chromosome_name",
                   "start_position", "end_position", "description")

    idx <- match(df$transcriptID, genemap$ensembl_transcript_id_version)
    df  <- merge(
      df,
      genemap[idx, keep_cols],
      by.x = "transcriptID",
      by.y = "ensembl_transcript_id_version"
    )
    names(df) <- c("transcriptID", final_col_names)

  } else {

    df <- data.frame(geneID = ensembl_ids)

    genemap <- tryCatch(
      biomaRt::getBM(
        attributes = annotations_attrs,
        filters    = "ensembl_gene_id",
        values     = df$geneID,
        mart       = ensembl
      ),
      error = function(e) .handle_biomart_query_error(e, host)
    )

    genemap <- .resolve_symbol_column(genemap, is_human, symbol_info$attr)

    keep_cols <- c("ensembl_gene_id", "symbol", "gene_biotype",
                   "chromosome_name", "start_position",
                   "end_position",    "description")

    idx <- match(df$geneID, genemap$ensembl_gene_id)
    df  <- merge(
      df,
      genemap[idx, keep_cols],
      by.x = "geneID",
      by.y = "ensembl_gene_id"
    )
    names(df) <- final_col_names
  }

  # --- Compute gene length ----------------------------------------------------
  df$gene_length <- df$gene_end - df$gene_start + 1
  df <- df %>% dplyr::relocate("gene_length", .before = "description")

  # --- Save to disk (optional) ------------------------------------------------
  if (!is.null(filename)) {

    if (format == "xlsx") {
      if (!requireNamespace("openxlsx", quietly = TRUE)) {
        stop(
          "Package \"openxlsx\" must be installed to save in xlsx format.\n",
          "Install it with: install.packages(\"openxlsx\")",
          call. = FALSE
        )
      }
      openxlsx::write.xlsx(df,
                           file     = paste0(filename, ".xlsx"),
                           colNames = TRUE,
                           rowNames = FALSE,
                           append   = FALSE)
    } else {
      utils::write.csv(df, row.names = FALSE, file = paste0(filename, ".csv"))
    }

    message("Annotations saved to: ", paste0(filename, ".", format))
  }

  return(df)
}

############################
# Function list_ensembl_versions #
############################

#' List available Ensembl releases and their BioMart hosts.
#'
#' A convenience wrapper around [biomaRt::listEnsemblArchives()] that returns
#' all available Ensembl release versions with their dates and host URLs.
#'
#' @details
#' Ensembl uses a **universal release numbering system**: every release (e.g.
#' v112, v113) covers *all* species simultaneously. What varies between species
#' is whether the species is present in a given release and which genome
#' assembly is used. To check whether your species of interest is available in
#' a particular version, use [list_ensembl_species()].
#'
#' The recommended workflow is:
#' \enumerate{
#'   \item Find the dataset name for your species: `list_ensembl_species()`.
#'   \item Find available versions and confirm species presence:
#'         `list_ensembl_versions()` + `list_ensembl_species(version = "X")`.
#'   \item Annotate: `get_annotations(species = "...", version = "...")`.
#' }
#'
#' @return A data frame with columns `name`, `date`, `url`, `version`, and
#'   `current_release` (marked with `*` for the active release).
#'
#' @examples
#' \dontrun{
#' # See all available Ensembl versions
#' list_ensembl_versions()
#'
#' # Then check if your species is in a specific version
#' list_ensembl_species(version = "112")
#' }
#'
#' @seealso [list_ensembl_species()], [get_annotations()]
#' @export

list_ensembl_versions <- function() {

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "Package \"biomaRt\" must be installed to use this function.\n",
      "Install it with: BiocManager::install(\"biomaRt\")",
      call. = FALSE
    )
  }

  archives <- tryCatch(
    biomaRt::listEnsemblArchives(),
    error = function(e) {
      stop(
        "Could not retrieve Ensembl archive list.\n",
        "This is likely a connectivity issue with the Ensembl server.\n\n",
        "Suggestions:\n",
        "  - Check your internet connection.\n",
        "  - Try again later.\n",
        "  - Visit https://www.ensembl.org/info/website/archives/index.html",
        call. = FALSE
      )
    }
  )

  df <- archives[, c("name", "date", "url", "version", "current_release")]
  # Exclude the legacy GRCh37 entry — it is a fixed human-only site,
  # not a versioned multi-species release.
  df <- df[df$version != "GRCh37", ]
  rownames(df) <- NULL

  message(
    "NOTE: Ensembl release numbers are UNIVERSAL across all species.\n",
    "A species may not exist in every version. Use list_ensembl_species(version)\n",
    "to confirm that your species dataset is available in a given release."
  )

  return(df)
}


############################
# Function list_ensembl_species #
############################

#' List all species datasets available in a given Ensembl version.
#'
#' Queries BioMart at the specified Ensembl version and returns a data frame
#' of all available gene datasets. Use this to:
#' \itemize{
#'   \item Find the exact `species` string to pass to [get_annotations()].
#'   \item Confirm that a species is available in a specific Ensembl version.
#' }
#'
#' @details
#' Because Ensembl releases are universal (same version number for all
#' species), this function simply lists which species datasets exist in the
#' BioMart of the requested version. Species added in later releases will not
#' appear in older versions.
#'
#' @param version Ensembl release version as a string (e.g. `"112"`, `"114"`),
#'   or `"current"` for the latest release. Default = `"current"`.
#' @param filter Optional. A search string to filter results by dataset name
#'   or description (case-insensitive). Useful for quickly finding a species
#'   without scrolling through hundreds of rows. Default = `NULL` (no filter).
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{dataset}{The dataset identifier to pass to `species` in
#'       [get_annotations()] (e.g. `"hsapiens_gene_ensembl"`).}
#'     \item{description}{Human-readable species name and assembly
#'       (e.g. `"Human genes (GRCh38.p14)"`).}
#'     \item{version}{Assembly version string.}
#'   }
#'
#' @examples
#' \dontrun{
#' # All species in version 112
#' list_ensembl_species(version = "112")
#'
#' # Find zebrafish dataset in version 112
#' list_ensembl_species(version = "112", filter = "zebrafish")
#'
#' # Find mouse dataset
#' list_ensembl_species(filter = "mouse")
#'
#' # Find all fish species
#' list_ensembl_species(filter = "fish")
#' }
#'
#' @seealso [list_ensembl_versions()], [get_annotations()]
#' @export

list_ensembl_species <- function(version = "current", filter = NULL) {

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "Package \"biomaRt\" must be installed to use this function.\n",
      "Install it with: BiocManager::install(\"biomaRt\")",
      call. = FALSE
    )
  }

  host <- .resolve_ensembl_host(version)

  message("Fetching species list from Ensembl version ", version, "...")

  mart <- tryCatch(
    biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host = host),
    error = function(e) .handle_biomart_connection_error(e, host)
  )

  datasets <- tryCatch(
    biomaRt::listDatasets(mart),
    error = function(e) .handle_biomart_query_error(e, host)
  )

  # Keep only gene datasets
  df <- datasets[grepl("_gene_ensembl$", datasets$dataset), ]
  rownames(df) <- NULL

  # Apply optional text filter across dataset name and description
  if (!is.null(filter) && nchar(trimws(filter)) > 0) {
    mask <- grepl(filter, df$dataset,     ignore.case = TRUE) |
      grepl(filter, df$description, ignore.case = TRUE)
    df <- df[mask, ]

    if (nrow(df) == 0) {
      message(
        "No species found matching \"", filter, "\" in Ensembl version ",
        version, ".\n",
        "Try a broader search term or check list_ensembl_species(version = \"",
        version, "\") without a filter."
      )
      return(invisible(df))
    }
  }

  message(
    nrow(df), " species datasets found",
    if (!is.null(filter)) paste0(" matching \"", filter, "\"") else "",
    " (Ensembl v", version, ").\n",
    "Pass the `dataset` column value to the `species` argument of get_annotations()."
  )

  return(df)
}


# =============================================================================
# Internal helpers
# =============================================================================

# Resolve BioMart host URL from version string
.resolve_ensembl_host <- function(version) {

  if (tolower(as.character(version)) == "current") {
    return("https://www.ensembl.org")
  }

  archives <- tryCatch(
    biomaRt::listEnsemblArchives(),
    error = function(e) {
      stop(
        "Could not retrieve the list of Ensembl archives.\n",
        "This may be due to a network or server issue.\n\n",
        "Suggestions:\n",
        "  - Check your internet connection.\n",
        "  - Use version = \"current\" to try the main server.\n",
        "  - Visit https://www.ensembl.org/info/website/archives/index.html",
        call. = FALSE
      )
    }
  )

  matched <- archives[archives$version == as.character(version), ]

  if (nrow(matched) == 0) {
    available <- paste(
      sort(archives$version[archives$version != "GRCh37"]),
      collapse = ", "
    )
    stop(
      "Ensembl version \"", version, "\" was not found.\n",
      "Available versions: ", available, "\n",
      "Use list_ensembl_versions() to see the full list.",
      call. = FALSE
    )
  }

  return(matched$url[1])
}


# Detect the best symbol attribute for a species
# Returns a list: list(attr = "attribute_name", source = "label for messages")
.get_symbol_attribute <- function(mart, species) {

  if (grepl("^hsapiens", species)) {
    return(list(attr = "hgnc_symbol", source = "HGNC"))
  }

  # Curated map: species prefix → preferred symbol attribute + source label
  curated_map <- list(
    mmusculus     = list(attr = "mgi_symbol",        source = "MGI"),
    rnorvegicus   = list(attr = "rgd_symbol",        source = "RGD"),
    drerio        = list(attr = "zfin_id_symbol",    source = "ZFIN"),
    dmelanogaster = list(attr = "flybasename_gene",  source = "FlyBase"),
    celegans      = list(attr = "wormbase_gene",     source = "WormBase")
  )

  # Check if this species matches any curated prefix
  for (prefix in names(curated_map)) {
    if (grepl(paste0("^", prefix), species)) {
      entry <- curated_map[[prefix]]

      # Verify the attribute actually exists in this mart/version
      available_attrs <- tryCatch(
        biomaRt::listAttributes(mart)$name,
        error = function(e) character(0)
      )

      if (entry$attr %in% available_attrs) {
        return(entry)
      } else {
        message(
          "Note: preferred symbol attribute '", entry$attr, "' (",
          entry$source, ") is not available for this species/version.\n",
          "Falling back to 'external_gene_name'."
        )
        return(list(attr = "external_gene_name", source = "External"))
      }
    }
  }

  # Default fallback for non-curated species
  return(list(attr = "external_gene_name", source = "External"))
}


# Resolve symbol column (human: HGNC fallback to external; others: rename attr)
.resolve_symbol_column <- function(genemap, is_human, symbol_attr) {

  if (is_human) {
    # Ensure both columns exist
    if (!"hgnc_symbol" %in% names(genemap)) {
      genemap$hgnc_symbol <- NA_character_
    }
    if (!"external_gene_name" %in% names(genemap)) {
      genemap$external_gene_name <- NA_character_
    }

    # Prefer HGNC, fallback to external gene name
    genemap$symbol <- ifelse(
      is.na(genemap$hgnc_symbol) | genemap$hgnc_symbol == "",
      genemap$external_gene_name,
      genemap$hgnc_symbol
    )
  } else {
    # Non-human: rename preferred attribute to "symbol" if present
    if (symbol_attr %in% names(genemap)) {
      genemap$symbol <- genemap[[symbol_attr]]
    } else if ("external_gene_name" %in% names(genemap)) {
      genemap$symbol <- genemap[["external_gene_name"]]
    } else {
      genemap$symbol <- NA_character_
    }
  }

  return(genemap)
}


# Standardized connection error handler
.handle_biomart_connection_error <- function(e, host, species = NULL, version = NULL) {

  ctx <- ""
  if (!is.null(version)) {
    ctx <- paste0(ctx, "version: ", version)
  }
  if (!is.null(species)) {
    ctx <- paste0(ctx, if (nchar(ctx) > 0) " | " else "", "species: ", species)
  }
  if (nchar(ctx) > 0) ctx <- paste0(" [", ctx, "]")

  stop(
    "Could not connect to Ensembl BioMart", ctx, ".\n",
    "Host: ", host, "\n\n",
    "Suggestions:\n",
    "  - Check your internet connection.\n",
    "  - Try again later.\n",
    "  - Use version = \"current\" to try the main server.\n",
    "  - Visit https://www.ensembl.org/info/website/archives/index.html",
    call. = FALSE
  )
}


# Standardized query error handler
.handle_biomart_query_error <- function(e, host) {
  stop(
    "BioMart query failed.\n",
    "Host: ", host, "\n\n",
    "Suggestions:\n",
    "  - Check your internet connection.\n",
    "  - Try again later.\n",
    "  - Confirm that the species dataset exists in this version.\n",
    "  - Use list_ensembl_species(version = \"X\") to verify.",
    call. = FALSE
  )
}
