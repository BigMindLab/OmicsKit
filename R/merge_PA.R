#######################
# Function merge_PA   #
#######################

utils::globalVariables(c(
  "LEADING EDGE", "tags", "signal",
  "FDR q-val", "Log10FDR", "FWER p-val", "Comparison", "...12"
))

#' Merge GSEA result files into a single data frame
#'
#' Reads all `.tsv` files produced by `GSEA_merge.sh` (from the GSEA.sh
#' pipeline) from a directory, standardizes numeric columns, parses the
#' leading edge string, computes `-log10(FDR)`, and returns a single merged
#' data frame ready for downstream use with [splot_PA()], [multiplot_PA()] , [getgenesPA()], and
#' [heatmap_PA()].
#'
#' **Input file format:** Each `.tsv` file corresponds to one MSigDB collection
#' (e.g., `H.tsv`, `C2.tsv`) and must follow the standard GSEA output
#' format with the following columns:
#'   * `NAME`: gene set name.
#'   * `SIZE`: number of genes in the gene set.
#'   * `ES`: enrichment score.
#'   * `NES`: normalized enrichment score.
#'   * `NOM p-val`: nominal p-value.
#'   * `FDR q-val`: false discovery rate. Values of exactly `0` indicate
#'     that no permutation produced an equally extreme NES (i.e., the true
#'     FDR is below the permutation resolution `1 / n_permutations`). These
#'     are replaced by `fdr_replace` to avoid `-Inf` in log-transforms.
#'   * `FWER p-val`: family-wise error rate.
#'   * `RANK AT MAX`: gene rank at maximum enrichment score.
#'   * `LEADING EDGE`: string encoding the leading edge statistics in the
#'     format `"tags=XX%, list=XX%, signal=XX%"`. Parsed into three numeric
#'     columns: `tags` (fraction of gene set in leading edge), `list`
#'     (fraction of ranked list used), and `signal` (enrichment signal
#'     strength).
#'   * `Comparison`: name of the comparison (e.g., `"TumorVsNormal"`).
#'     Renamed to `COMPARISON` in the output. Required for visualization
#'     with [splot_PA()] or [multiplot_PA()] .
#'   * `GS<br> follow link to MSigDB` and `GS DETAILS`: removed automatically.
#'
#' **Output columns:** All input columns (minus the two removed above) plus:
#'   * `COLLECTION`: name of the MSigDB collection, derived from the file name
#'     by removing the `.tsv` suffix.
#'   * `tags`, `list`, `signal`: numeric leading edge components (0-1 scale).
#'   * `Log10FDR`: `-log10(FDR)` computed after applying `fdr_replace`.
#'   * `FDR`: renamed from `FDR q-val`.
#'   * `COMPARISON`: renamed from `Comparison`.
#'
#' @param input_directory Character. Path to the directory containing one or
#'   more GSEA collection result files in `.tsv` format (e.g., output of
#'   `GSEA_merge.sh`). Each file must end in `.tsv`.
#' @param fdr_replace Numeric. Value used to replace `FDR q-val = 0`. This
#'   occurs when no permutation produced an NES as extreme as the observed
#'   one, meaning the true FDR is below `1 / n_permutations`. With the
#'   standard 1,000 permutations, the recommended value is `0.001`. Adjust
#'   to `1 / n_permutations` if a different number of permutations was used.
#'   Default: `0.001`.
#'
#' @note The input `.tsv` files must contain a `Comparison` column identifying
#'   each comparison (e.g., `"TumorVsNormal"`). This column is renamed to
#'   `COMPARISON` in the output and is required by [splot_PA()] or [multiplot_PA()]  to operate in
#'   multi-comparison mode. If your files come from a single comparison and
#'   do not have this column, add it manually to each file before merging:
#'   `your_data$Comparison <- "YourComparisonName"`.
#'
#' @return A data frame (`gsea_data`) containing all merged and processed GSEA
#'   results with standardized column names.
#'
#' @examples
#' \dontrun{
#' # Merge all GSEA collection TSV files from a directory
#' gsea_data <- merge_PA(
#'   input_directory = "path/to/gsea_results/",
#'   fdr_replace     = 0.001   # standard for 1000 permutations
#' )
#'
#' # Inspect result
#' head(gsea_data)
#' colnames(gsea_data)
#'
#' # Use directly in downstream functions
#' gsl        <- list_gmts("path/to/gmt_folder/")
#' ranked     <- deseq2_results$gene_id[order(deseq2_results$stat,
#'                                            decreasing = TRUE)]
#' gene_lists <- getgenesPA(gsea_data, gsl, ranked)
#' pa_annot   <- addgenesPA(gsea_data, gene_lists)
#'
#' plot_PA(gsea_data, comparison_col = "COMPARISON")
#' }
#'
#' @seealso [splot_PA()] or [multiplot_PA()]  for visualization of merged results;
#'   [getgenesPA()] and [addgenesPA()] for gene-level annotation;
#'   [heatmap_PA()] for leading edge heatmaps;
#'   [list_gmts()] to load gene sets.
#'
#' @importFrom magrittr %>%
#' @export

merge_PA <- function(input_directory,
                     fdr_replace = 0.001) {

  if (!requireNamespace("readr",      quietly = TRUE)) {
    stop("Package \"readr\" must be installed to use this function.",      call. = FALSE)
  }
  if (!requireNamespace("tidyr",      quietly = TRUE)) {
    stop("Package \"tidyr\" must be installed to use this function.",      call. = FALSE)
  }
  if (!requireNamespace("tidyselect", quietly = TRUE)) {
    stop("Package \"tidyselect\" must be installed to use this function.", call. = FALSE)
  }

  # --- Input validation --
  if (!dir.exists(input_directory)) {
    stop("The specified directory does not exist: ", input_directory, call. = FALSE)
  }

  files <- list.files(path = input_directory, pattern = "\\.tsv$",
                      full.names = TRUE)
  if (length(files) == 0) {
    stop("No TSV files found in: ", input_directory, call. = FALSE)
  }

  if (!is.numeric(fdr_replace) || fdr_replace <= 0 || fdr_replace >= 1) {
    stop("`fdr_replace` must be a numeric value between 0 and 1.", call. = FALSE)
  }

  # --- Numeric columns present in every GSEA output file
  numeric_cols <- c("SIZE", "ES", "NES", "NOM p-val",
                    "FDR q-val", "FWER p-val", "RANK AT MAX")

  # --- Read and combine all TSV files
  read_one <- function(file) {
    data      <- readr::read_tsv(file, show_col_types = FALSE)
    coll_name <- sub("\\.tsv$", "", basename(file))
    data <- data %>%
      dplyr::mutate(
        dplyr::across(tidyselect::all_of(numeric_cols), as.numeric)
      )
    data$COLLECTION <- coll_name
    return(data)
  }

  gsea_data <- lapply(files, read_one) %>%
    dplyr::bind_rows()

  # Remove empty trailing column if present
  if ("...12" %in% colnames(gsea_data)) {
    gsea_data <- dplyr::select(gsea_data, -`...12`)
  }

  # --- Report any non-numeric values in numeric columns --
  problematic <- gsea_data %>%
    dplyr::filter(
      dplyr::if_any(tidyselect::all_of(numeric_cols),
                    ~ !grepl("^-?[0-9.]+([eE][+-]?[0-9]+)?$",
                             as.character(.)))
    )
  if (nrow(problematic) > 0) {
    message("Warning: ", nrow(problematic),
            " row(s) with non-numeric values in numeric columns:")
    print(problematic)
  }

  # --- Check Comparison column ---
  if (!"Comparison" %in% colnames(gsea_data)) {
    stop(
      "No 'Comparison' column found in the TSV files. ",
      "Add it manually to each file before merging: ",
      "your_data$Comparison <- 'YourComparisonName'.",
      call. = FALSE
    )
  }

  # --- Process and standardize ---
  gsea_data <- gsea_data %>%
    dplyr::select(-dplyr::any_of(c("GS<br> follow link to MSigDB",
                                   "GS DETAILS"))) %>%
    tidyr::separate(
      col    = `LEADING EDGE`,
      into   = c("tags", "list", "signal"),
      sep    = ",",
      remove = FALSE
    ) %>%
    dplyr::mutate(
      tags    = 0.01 * as.numeric(sub("%", "", sub("tags=",   "", tags))),
      list    = 0.01 * as.numeric(sub("%", "", sub("list=",   "", list))),
      signal  = 0.01 * as.numeric(sub("%", "", sub("signal=", "", signal))),
      # Replace FDR = 0 to avoid -Inf in log-transforms.
      # FDR = 0 means no permutation matched the observed NES,
      # so the true FDR < 1/n_permutations (typically 0.001 for 1000 perms).
      `FDR q-val` = ifelse(`FDR q-val` == 0, fdr_replace, `FDR q-val`),
      Log10FDR    = -log10(`FDR q-val`)
    ) %>%
    dplyr::relocate(Log10FDR, .after = `FWER p-val`) %>%
    dplyr::rename(COMPARISON = Comparison, FDR = `FDR q-val`)

  # --- Message and return
  message("Merged GSEA data processed: ",
          nrow(gsea_data), " gene sets from ",
          length(unique(gsea_data$COLLECTION)), " collection(s) and ",
          length(unique(gsea_data$COMPARISON)), " comparison(s).")

  return(gsea_data)
}
