#######################
# Function merge_GSEA #
#######################

#' Merge GSEA results data frames.
#'
#' After running GSEA_all.sh from GSEA.sh, merge_GSEA function joins .tsv files to a single file
#'
#' @param input_directory The directory containing the GSEA collection results in TSV format.
#' @param output_file The output file to save the merged data. If not provided, the file will be saved in the input directory.
#' @importFrom magrittr %>%
#' @export


merge_GSEA <- function(input_directory, output_file = "collections_merged_gsea_data.tsv") {

  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package \"dplyr\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("readr", quietly = TRUE)) stop("Package \"readr\" must be installed to use this function.", call. = FALSE)
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Package \"tidyr\" must be installed to use this function.", call. = FALSE)

  # Validate input directory and check for TSV files
  if (!dir.exists(input_directory)) {
    stop("The specified directory does not exist: ", input_directory)
  }
  files <- list.files(path = input_directory, pattern = "\\.tsv$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No TSV files found in ", input_directory)
  }

  # Function to read each file and add a column with the modified file name
  read_file <- function(file) {
    data <- readr::read_tsv(file)
    file_name <- basename(file)
    file_name <- sub("_all.tsv$", "", file_name)  # Change the pattern if necessary
    numeric_cols <- c("SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "RANK AT MAX")
    data <- data %>%
      dplyr::mutate(dplyr::across(tidyselect::all_of(numeric_cols), as.numeric))
    data$COLLECTION <- file_name
    return(data)
  }

  # Read and combine all TSV files into a single data frame, remove empty columns
  gsea_data <- lapply(files, read_file) %>% dplyr::bind_rows() %>% dplyr::select(-`...12`)

 # Find problematic values in numeric columns
  gsea_data %>%
    dplyr::filter(dplyr::if_any(tidyselect::all_of(numeric_cols), ~ !grepl("^-?[0-9.]+$", .))) %>%
    print()

  # Data processing: selection, separation, mutation, and renaming of columns
  gsea_data <- gsea_data %>%
    dplyr::select(-"GS<br> follow link to MSigDB", -"GS DETAILS") %>%
    tidyr::separate(col = `LEADING EDGE`, into = c("tags", "list", "signal"), sep = ",", remove = FALSE) %>%
    dplyr::mutate(
      tags = 0.01 * as.numeric(sub("%", "", sub("tags=", "", tags))),
      list = 0.01 * as.numeric(sub("%", "", sub("list=", "", list))),
      signal = 0.01 * as.numeric(sub("%", "", sub("signal=", "", signal))),
      `FDR q-val` = ifelse(`FDR q-val` == 0, 0.001, `FDR q-val`),
      `Log10FDR` = -log10(`FDR q-val`)
    ) %>%
    dplyr::relocate(`Log10FDR`, .after = `FWER p-val`) %>%
    dplyr::rename(COMPARISON = Comparison, FDR = `FDR q-val`)

  # Save the processed data to a TSV file
  readr::write_tsv(gsea_data, output_file)
  message("GSEA data saved to:", output_file, "\n")

  return(TRUE)
}
