############################
# Function get_annotations #
############################

#' Get annotations from Ensembl.
#'
#' This function annotates a column of transcripts or gene IDs (ENSEMBL) with information of the Biomart.
#' If transcript IDs are provided, they are also annotated with information of the genes to which they belong.
#'
#' The Gene information added include:
#' - Gene ENSEMBL ID, HGNC Symbol, Description, Biotype and Chromosome.
#' - Gene start, end and length
#'
#' @param ensembl_ids The column of transcripts to be used as input.
#' @param mode To specify the IDs provided, between "transcripts" or "genes". Default = genes.
#' @param version This function can use the versions 102, 103, and 112 of Ensembl. Default = "Current".
#' @param filename The name of the output file, which is table. Default = gene_annotations.
#' @param format The output is saved in .csv or .xlsx formats. Default = csv.
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export

get_annotations <- function(ensembl_ids, mode = "genes", filename = "gene_annotations", version = "Current", format = "csv") {

  if (!requireNamespace("biomaRt", quietly = TRUE)) stop("Package \"biomaRt\" must be installed to use this function.", call. = FALSE)

  if (version == "102") {
    ensembl = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                               dataset = "hsapiens_gene_ensembl",
                               host = "https://nov2020.archive.ensembl.org")
  } else if (version == "103") {
    ensembl = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                               dataset = "hsapiens_gene_ensembl",
                               host = "https://feb2021.archive.ensembl.org")
  } else if (version == "112") {
    ensembl = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                               dataset = "hsapiens_gene_ensembl",
                               host = "https://may2024.archive.ensembl.org")
  } else if (version == "Current") {
    ensembl = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                               dataset = "hsapiens_gene_ensembl",
                               host = "https://useast.ensembl.org")
  }

  annotations <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype",
                   "chromosome_name", "start_position", "end_position",
                   "description")

  new_names <- c("geneID", "symbol", "biotype", "chromosome",
                 "gene_start", "gene_end", "description")

  # There are more annotations available in the biomaRt, check listAttributes(mart = ensembl)
  # The terms "go_id" and "name_1006" can be added in a future release.

  if (mode == "transcripts") {

    df <- data.frame(transcriptID = ensembl_ids)

    genemap <- biomaRt::getBM(attributes = c("ensembl_transcript_id_version", annotations),
                              filters = "ensembl_transcript_id_version",
                              values = df$transcriptID,
                              mart = ensembl)

    idx <- match(df$transcriptID, genemap$ensembl_transcript_id_version)
    df <- merge(df, genemap[idx, c("ensembl_transcript_id_version", annotations)],
                by.x = "transcriptID", by.y = "ensembl_transcript_id_version")
    names(df) <- c("transcriptID", new_names)

  } else {

    df <- data.frame(geneID = ensembl_ids)

    genemap <- biomaRt::getBM(attributes = annotations,
                              filters = "ensembl_gene_id",
                              values = df$geneID,
                              mart = ensembl)

    idx <- match(df$geneID, genemap$ensembl_gene_id)
    df <- merge(df, genemap[idx, annotations], by.x = "geneID", by.y = "ensembl_gene_id")
    names(df) <- new_names

  }

  df$gene_length <- df$gene_end - df$gene_start + 1
  df <- df %>% dplyr::relocate(.data$gene_length, .before = "description")

  if (format == "xlsx") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop(
        "Package \"openxlsx\" must be installed to use this function.",
        call. = FALSE
        )
    }

    openxlsx::write.xlsx(df, file = paste0(filename, ".xlsx"), colNames = T, rowNames = F, append = F)
  } else {
    utils::write.csv(df, rowNames = F, file = paste0(filename, ".csv"))
  }

  return(df)
}
