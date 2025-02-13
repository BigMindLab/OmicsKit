
# Function merge_GSEA #  

#' Merge GSEA results data frames.
#' 
#' After run GSEA_all.sh from GSEA.sh, move the .tsv files to a single directory 
#' 
#' @param input_directory The directory containing the GSEA collection results in TSV format.
#' @param output_file The output file to save the merged data. If not provided, the file will be saved in the input directory.
#' @importFrom magrittr %>%
#' @export
 

merge_GSEA <- function(input_directory, output_file) {
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }
 
  # Read command-line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: merge_GSEA.R <input_directory> [output_file]\n")
  }
  input_directory <- args[1]
  # If a second argument is provided, use it as the output file; otherwise, save the file in the input directory.
  output_file <- ifelse(length(args) >= 2, args[2], file.path(input_directory, "collections_merged_gsea_data.tsv"))
  
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
    # Convert all numeric columns explicitly
    numeric_cols <- c("SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "RANK AT MAX")
    data <- data %>%
      dplyr::mutate(across(tidyselect::all_of(numeric_cols), as.numeric))
    data$COLLECTION <- file_name
    return(data)
  }

  # Read and combine all TSV files into a single data frame
  gsea_data <- lapply(files, read_file) %>% dplyr::bind_rows()
  
  gsea_data <- gsea_data %>%
    dplyr::select(-`...12`)
  
  # Define numeric columns
  numeric_cols <- c("SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "RANK AT MAX")
  
  # Find problematic values in numeric columns
  gsea_data %>%
    dplyr::filter(dplyr::if_any(all_of(numeric_cols), ~ !grepl("^-?[0-9.]+$", .))) %>%
    print()
  
  # Data processing: selection, separation, mutation, and renaming of columns
  gsea_data <- gsea_data %>%
    dplyr::select(-"GS<br> follow link to MSigDB", -"GS DETAILS") %>%  #, -"...12"
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
  
  cat("GSEA data saved to:", output_file, "\n")
}



