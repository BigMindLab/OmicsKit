######################
# Function list_gmts #
######################

#' Read GMT files from a directory into a named gene set list
#'
#' Scans a directory for `.gmt` files, parses them, and returns a single named
#' list where each element is a character vector of gene symbols for one gene
#' set. The output is ready to be passed directly to [calc_jaccard()].
#'
#' **GMT format:** each row contains the gene set name in column 1, an optional
#' description in column 2, and gene symbols from column 3 onward. Empty fields
#' are automatically removed. This is the standard format used by MSigDB
#' (Molecular Signatures Database) and other gene set annotation resources.
#'
#' @param dir Character. Path to the directory containing one or more `.gmt`
#'   files. The function searches the directory non-recursively.
#'
#' @return A named list where each element is a character vector of gene symbols
#'   for one gene set. Names correspond to gene set names as defined in column 1
#'   of the GMT files. If the same gene set name appears in multiple files, the
#'   last occurrence overwrites the earlier one.
#'
#' @examples
#' \dontrun{
#' # Read all GMT files from a directory
#' geneset_list <- list_gmts("path/to/gmt_files/")
#'
#' # Inspect output
#' length(geneset_list)              # number of gene sets
#' names(geneset_list)[1:5]          # first five gene set names
#' geneset_list[["KEGG_APOPTOSIS"]]  # genes in a specific set
#'
#' # Pass directly to calc_jaccard
#' jac <- calc_jaccard(geneset_list, results_df, fdr_th = 0.05)
#' }
#'
#' @seealso [calc_jaccard()]
#' @export

list_gmts <- function(dir) {

  if (!is.character(dir) || length(dir) != 1) {
    stop("`dir` must be a single character string with the path to a directory.",
         call. = FALSE)
  }
  if (!dir.exists(dir)) {
    stop("Directory not found: ", dir, call. = FALSE)
  }

  gmt_files <- list.files(dir, pattern = "\\.gmt$", full.names = TRUE)

  if (length(gmt_files) == 0) {
    stop("No .gmt files found in: ", dir, call. = FALSE)
  }

  geneset_list <- list()

  for (f in gmt_files) {
    gmt <- utils::read.delim(f, header = FALSE, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(gmt))) {
      name  <- gmt$V1[i]
      genes <- as.character(gmt[i, 3:ncol(gmt)])
      genes <- genes[genes != ""]
      geneset_list[[name]] <- genes
    }
  }

  if (length(geneset_list) == 0) {
    stop("No gene sets could be parsed from the .gmt files in: ", dir,
         call. = FALSE)
  }

  message("Loaded ", length(geneset_list), " gene sets from ",
          length(gmt_files), " GMT file(s).")

  return(geneset_list)
}
