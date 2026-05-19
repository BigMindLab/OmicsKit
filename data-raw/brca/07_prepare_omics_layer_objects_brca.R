## 07_prepare_omics_layer_objects_brca.R
##
## Purpose:
##   Prepare real sample-specific TCGA-BRCA omics-layer track files and compact
##   metadata objects for nice_GenomeTrack() and nice_circos().
##
## Inputs:
##   - inst/extdata/brca/raw_xena/TCGA.BRCA.sampleMap_HiSeqV2/HiSeqV2
##   - inst/extdata/brca/raw_xena/TCGA.BRCA.sampleMap_SNP6_nocnv_genomicSegment/SNP6_nocnv_genomicSegment
##   - inst/extdata/brca/raw_xena/TCGA.BRCA.sampleMap_HumanMethylation450/HumanMethylation450
##   - inst/extdata/brca/raw_xena/mc3_BRCA_mc3.txt/BRCA_mc3.txt
##   - inst/extdata/brca/raw_xena/hg19ToHg38.over.chain/hg19ToHg38.over.chain
##   - inst/extdata/brca/raw_xena/hugo_probemap.tsv
##   - inst/extdata/brca/raw_xena/probemap_450k.tsv
##
## Outputs:
##   - Sample-specific BigWig and BED files under inst/extdata/brca/omics_layers/
##   - data/brca_omics_layer_tracks.rda
##
## Notes:
##   This script prepares data only. Figure rendering belongs in
##   data-raw/brca/09_render_brca_figures.R.

options(stringsAsFactors = FALSE)

sample_id <- "TCGA-A1-A0SH-01"

standard_chrs <- paste0("chr", c(1:22, "X", "Y"))
brca1_chr <- "chr17"
brca1_start_hg19 <- 41196312 - 500000
brca1_end_hg19 <- 41277500 + 500000
brca1_region_hg38_for_plot <- c(chr = "chr17", start = 42944295, end = 43610338)

hg38_lengths <- c(
  chr1 = 248956422, chr2 = 242193529, chr3 = 198295559,
  chr4 = 190214555, chr5 = 181538259, chr6 = 170805979,
  chr7 = 159345973, chr8 = 145138636, chr9 = 138394717,
  chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
  chr13 = 114364328, chr14 = 107043718, chr15 = 101991189,
  chr16 = 90338345, chr17 = 83257441, chr18 = 80373285,
  chr19 = 58617616, chr20 = 64444167, chr21 = 46709983,
  chr22 = 50818468, chrX = 156040895, chrY = 57227415
)

project_file <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    here::here(...)
  } else {
    file.path(...)
  }
}

raw_dir <- project_file("inst", "extdata", "brca", "raw_xena")
omics_dir <- project_file("inst", "extdata", "brca", "omics_layers")

rna_file <- file.path(raw_dir, "TCGA.BRCA.sampleMap_HiSeqV2", "HiSeqV2")
cnv_file <- file.path(
  raw_dir,
  "TCGA.BRCA.sampleMap_SNP6_nocnv_genomicSegment",
  "SNP6_nocnv_genomicSegment"
)
methyl_file <- file.path(
  raw_dir,
  "TCGA.BRCA.sampleMap_HumanMethylation450",
  "HumanMethylation450"
)
mutation_file <- file.path(raw_dir, "mc3_BRCA_mc3.txt", "BRCA_mc3.txt")
chain_file <- file.path(raw_dir, "hg19ToHg38.over.chain", "hg19ToHg38.over.chain")
hugo_probemap_file <- file.path(raw_dir, "hugo_probemap.tsv")
probemap_450k_file <- file.path(raw_dir, "probemap_450k.tsv")

rna_bw_file <- file.path(omics_dir, paste0(sample_id, "_rna_expr_hg38.bw"))
rna_1mb_bw_file <- file.path(omics_dir, paste0(sample_id, "_rna_expr_1Mb_hg38.bw"))
rna_1mb_bed_file <- file.path(omics_dir, paste0(sample_id, "_rna_expr_1Mb_hg38.bed"))
cnv_bw_file <- file.path(omics_dir, paste0(sample_id, "_cnv_hg38.bw"))
cnv_bed_file <- file.path(omics_dir, paste0(sample_id, "_cnv_hg38.bed"))
methyl_brca1_bw_file <- file.path(omics_dir, paste0(sample_id, "_methyl_brca1_hg38.bw"))
methyl_1mb_bw_file <- file.path(omics_dir, paste0(sample_id, "_methyl_450k_meanBeta_1Mb_hg38.bw"))
methyl_1mb_bed_file <- file.path(omics_dir, paste0(sample_id, "_methyl_450k_meanBeta_1Mb_hg38.bed"))
mutation_bed_file <- file.path(omics_dir, paste0(sample_id, "_mutations_brca1_hg38.bed"))

require_file <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required raw file is missing for BRCA omics-layer preparation:\n",
      normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }

  invisible(path)
}

for (path in c(
  rna_file, cnv_file, methyl_file, mutation_file, chain_file,
  hugo_probemap_file, probemap_450k_file
)) {
  require_file(path)
}

required_packages <- c(
  "GenomicRanges", "GenomeInfoDb", "IRanges", "S4Vectors",
  "rtracklayer", "usethis"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package \"", pkg, "\" is required for BRCA omics-layer preparation.",
      call. = FALSE
    )
  }
}

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

clean_chr <- function(x) {
  y <- trimws(as.character(x))
  y <- sub("^chr", "", y, ignore.case = TRUE)
  y <- toupper(y)
  paste0("chr", y)
}

clean_tcga_sample <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y <- gsub("\\.", "-", y)
  y[!has_value(y)] <- NA_character_

  ifelse(
    !is.na(y) & nchar(y) >= 16L,
    substr(y, 1L, 16L),
    ifelse(!is.na(y) & nchar(y) >= 15L, substr(y, 1L, 15L), NA_character_)
  )
}

sample_key <- function(x) {
  y <- clean_tcga_sample(x)
  ifelse(!is.na(y) & nchar(y) >= 15L, substr(y, 1L, 15L), NA_character_)
}

target_sample_key <- sample_key(sample_id)

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

as_numeric_clean <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

assign_hg38_seqinfo <- function(gr, standard_chrs) {
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  gr <- GenomeInfoDb::keepSeqlevels(
    gr,
    value = intersect(standard_chrs, GenomeInfoDb::seqlevels(gr)),
    pruning.mode = "coarse"
  )
  gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")

  seqinfo <- GenomeInfoDb::Seqinfo(
    seqnames = standard_chrs,
    seqlengths = hg38_lengths[standard_chrs],
    genome = "hg38"
  )

  GenomeInfoDb::seqinfo(gr) <- seqinfo[GenomeInfoDb::seqlevels(gr)]
  gr
}

validate_no_overlaps_for_bigwig <- function(gr) {
  gr <- GenomicRanges::sort(gr)
  overlaps <- GenomicRanges::findOverlaps(
    gr,
    drop.self = TRUE,
    drop.redundant = TRUE,
    ignore.strand = TRUE
  )

  if (length(overlaps) > 0L) {
    stop(
      "BigWig export requires non-overlapping intervals; overlaps remain after processing.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

resolve_overlaps_for_bigwig <- function(gr, agg = c("max", "mean")) {
  agg <- match.arg(agg)

  if (!"score" %in% names(S4Vectors::mcols(gr))) {
    if ("value" %in% names(S4Vectors::mcols(gr))) {
      gr$score <- gr$value
    } else {
      stop("GRanges must contain a numeric `score` column.", call. = FALSE)
    }
  }

  gr$score <- as_numeric_clean(gr$score)
  gr <- gr[is.finite(gr$score)]
  gr <- GenomicRanges::sort(gr)

  if (length(gr) == 0L) {
    stop("No finite signal values remain for overlap resolution.", call. = FALSE)
  }

  overlaps <- GenomicRanges::findOverlaps(
    gr,
    drop.self = TRUE,
    drop.redundant = TRUE,
    ignore.strand = TRUE
  )

  if (length(overlaps) == 0L) {
    return(gr)
  }

  disjoined <- GenomicRanges::disjoin(gr, ignore.strand = TRUE)
  hits <- GenomicRanges::findOverlaps(disjoined, gr, ignore.strand = TRUE)
  scores_by_query <- split(gr$score[S4Vectors::subjectHits(hits)], S4Vectors::queryHits(hits))

  score <- rep(NA_real_, length(disjoined))
  idx <- as.integer(names(scores_by_query))

  if (identical(agg, "max")) {
    score[idx] <- vapply(scores_by_query, max, numeric(1L), na.rm = TRUE)
  } else {
    score[idx] <- vapply(scores_by_query, mean, numeric(1L), na.rm = TRUE)
  }

  disjoined$score <- score
  disjoined <- disjoined[is.finite(disjoined$score)]
  GenomicRanges::sort(disjoined)
}

write_bed_from_granges <- function(gr, bed_file, mutation = FALSE) {
  name <- if ("name" %in% names(S4Vectors::mcols(gr))) {
    as.character(gr$name)
  } else {
    rep(".", length(gr))
  }

  score <- if ("score" %in% names(S4Vectors::mcols(gr))) {
    as_numeric_clean(gr$score)
  } else if ("value" %in% names(S4Vectors::mcols(gr))) {
    as_numeric_clean(gr$value)
  } else {
    rep(1, length(gr))
  }

  if (isTRUE(mutation)) {
    score[!is.finite(score)] <- 1
  }

  bed <- data.frame(
    chr = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr) - 1L,
    end = GenomicRanges::end(gr),
    name = name,
    score = score,
    stringsAsFactors = FALSE
  )

  utils::write.table(
    bed,
    file = bed_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    na = "."
  )

  invisible(bed_file)
}

export_signal_files <- function(gr, bed_file = NULL, bw_file = NULL) {
  if (!"score" %in% names(S4Vectors::mcols(gr))) {
    stop("Signal GRanges must contain a numeric `score` column.", call. = FALSE)
  }

  gr$score <- as_numeric_clean(gr$score)
  gr <- gr[is.finite(gr$score)]
  gr <- assign_hg38_seqinfo(gr, standard_chrs)
  gr <- GenomicRanges::sort(gr)

  if (length(gr) == 0L) {
    stop("No finite signal intervals remain for export.", call. = FALSE)
  }

  if (!is.null(bw_file)) {
    validate_no_overlaps_for_bigwig(gr)
    rtracklayer::export(gr, bw_file, format = "BigWig")
  }

  if (!is.null(bed_file)) {
    write_bed_from_granges(gr, bed_file, mutation = FALSE)
  }

  invisible(gr)
}

make_binned_signal <- function(gr, bin_size = 1e6, fun = "mean", min_n = 1) {
  if (!"score" %in% names(S4Vectors::mcols(gr))) {
    stop("GRanges must contain a numeric `score` column for binning.", call. = FALSE)
  }

  fun <- match.arg(fun, c("mean", "max"))
  gr <- assign_hg38_seqinfo(gr, standard_chrs)
  gr$score <- as_numeric_clean(gr$score)
  gr <- gr[is.finite(gr$score)]

  if (length(gr) == 0L) {
    stop("No finite intervals remain for binned signal generation.", call. = FALSE)
  }

  midpoint <- floor((GenomicRanges::start(gr) + GenomicRanges::end(gr)) / 2)
  chr <- as.character(GenomicRanges::seqnames(gr))
  bin_index <- floor((midpoint - 1) / bin_size)
  key <- paste(chr, bin_index, sep = "\t")
  split_scores <- split(gr$score, key)

  bin_score <- if (identical(fun, "mean")) {
    vapply(split_scores, mean, numeric(1L), na.rm = TRUE)
  } else {
    vapply(split_scores, max, numeric(1L), na.rm = TRUE)
  }

  bin_n <- vapply(split_scores, length, integer(1L))
  keep <- bin_n >= min_n & is.finite(bin_score)
  parts <- strsplit(names(split_scores)[keep], "\t", fixed = TRUE)

  out_chr <- vapply(parts, `[`, character(1L), 1L)
  out_bin <- as.integer(vapply(parts, `[`, character(1L), 2L))
  out_start <- out_bin * bin_size + 1
  out_end <- pmin((out_bin + 1) * bin_size, hg38_lengths[out_chr])

  out <- GenomicRanges::GRanges(
    seqnames = out_chr,
    ranges = IRanges::IRanges(start = out_start, end = out_end),
    score = as.numeric(bin_score[keep]),
    n = as.integer(bin_n[keep]),
    name = paste0(out_chr, ":", out_start, "-", out_end)
  )

  assign_hg38_seqinfo(out, standard_chrs)
}

load_chain_file <- function(chain_file) {
  require_file(chain_file)
  rtracklayer::import.chain(chain_file)
}

read_xena_sample_vector <- function(path, sample_id, value_name) {
  require_file(path)
  header <- strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]

  if (length(header) < 2L) {
    stop("Xena matrix has no sample columns: ", path, call. = FALSE)
  }

  matched <- which(sample_key(header) == target_sample_key)
  matched <- matched[matched != 1L]

  if (length(matched) == 0L) {
    stop(
      "Sample ", sample_id, " was not found in Xena matrix:\n",
      normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }

  col_classes <- rep("NULL", length(header))
  col_classes[1L] <- "character"
  col_classes[matched] <- "numeric"

  tbl <- utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    colClasses = col_classes
  )

  values <- as.matrix(tbl[, -1L, drop = FALSE])
  mode(values) <- "numeric"

  data.frame(
    id = trimws(as.character(tbl[[1L]])),
    value = rowMeans(values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )[c("id", "value")]
}

read_tsv <- function(path) {
  require_file(path)
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
}

read_probemap <- function(path, label) {
  tbl <- read_tsv(path)
  if ("#id" %in% names(tbl) && !"id" %in% names(tbl)) {
    names(tbl)[names(tbl) == "#id"] <- "id"
  }

  id_col <- find_column(tbl, c("id", "gene", "symbol", "name"), paste(label, "id"))
  chr_col <- find_column(tbl, c("chrom", "chr", "chromosome"), paste(label, "chromosome"))
  start_col <- find_column(
    tbl,
    c("chromStart", "start", "txStart", "probe_start", "position", "pos"),
    paste(label, "start")
  )
  end_col <- find_column(
    tbl,
    c("chromEnd", "end", "txEnd", "probe_end"),
    paste(label, "end"),
    required = FALSE
  )

  out <- data.frame(
    id = trimws(as.character(tbl[[id_col]])),
    chr = clean_chr(tbl[[chr_col]]),
    start = as_numeric_clean(tbl[[start_col]]),
    end = if (!is.null(end_col)) as_numeric_clean(tbl[[end_col]]) else as_numeric_clean(tbl[[start_col]]),
    stringsAsFactors = FALSE
  )

  out <- out[has_value(out$id) & out$chr %in% standard_chrs &
               is.finite(out$start) & is.finite(out$end) &
               out$start <= out$end, , drop = FALSE]
  out
}

make_granges_hg19 <- function(chr, start, end, score, name = NULL, extra = NULL) {
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = as.integer(start), end = as.integer(end)),
    score = as_numeric_clean(score)
  )

  if (!is.null(name)) {
    gr$name <- as.character(name)
  }

  if (!is.null(extra)) {
    for (nm in names(extra)) {
      S4Vectors::mcols(gr)[[nm]] <- extra[[nm]]
    }
  }

  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  gr <- gr[as.character(GenomicRanges::seqnames(gr)) %in% standard_chrs &
             is.finite(gr$score)]
  GenomicRanges::sort(gr)
}

lift_first <- function(gr, chain) {
  gr$..row_id <- seq_along(gr)
  lifted <- rtracklayer::liftOver(gr, chain)
  keep <- lengths(lifted) > 0L

  if (!any(keep)) {
    stop("No intervals mapped during hg19 to hg38 liftOver.", call. = FALSE)
  }

  first <- lapply(lifted[keep], function(x) x[1L])
  out <- unlist(GenomicRanges::GRangesList(first), use.names = FALSE)
  S4Vectors::mcols(out)$..row_id <- NULL
  assign_hg38_seqinfo(out, standard_chrs)
}

lift_all <- function(gr, chain) {
  lifted <- rtracklayer::liftOver(gr, chain)
  keep <- lengths(lifted) > 0L

  if (!any(keep)) {
    stop("No intervals mapped during hg19 to hg38 liftOver.", call. = FALSE)
  }

  out <- unlist(lifted[keep], use.names = FALSE)
  assign_hg38_seqinfo(out, standard_chrs)
}

filter_region <- function(df, chr, start, end) {
  df[df$chr == chr & df$start <= end & df$end >= start, , drop = FALSE]
}

dir.create(omics_dir, recursive = TRUE, showWarnings = FALSE)

chain <- load_chain_file(chain_file)

## RNA-seq ---------------------------------------------------------------------
message("Preparing RNA-seq signal track for ", sample_id, "...")
rna_sample <- read_xena_sample_vector(rna_file, sample_id, "rna_expr")
names(rna_sample) <- c("id", "score")
rna_sample$id <- sub("\\|.*$", "", rna_sample$id)

hugo_probemap <- read_probemap(hugo_probemap_file, "hugo_probemap")
rna_merged <- merge(rna_sample, hugo_probemap, by = "id", all = FALSE)

if (nrow(rna_merged) == 0L) {
  stop("RNA expression rows did not match hugo_probemap.tsv gene IDs.", call. = FALSE)
}

rna_hg19 <- make_granges_hg19(
  chr = rna_merged$chr,
  start = rna_merged$start,
  end = rna_merged$end,
  score = rna_merged$score,
  name = rna_merged$id
)
rna_hg38 <- lift_first(rna_hg19, chain)
rna_hg38 <- resolve_overlaps_for_bigwig(rna_hg38, agg = "max")
export_signal_files(rna_hg38, bw_file = rna_bw_file)

rna_1mb_hg38 <- make_binned_signal(rna_hg38, bin_size = 1e6, fun = "mean", min_n = 1)
rna_1mb_hg38 <- resolve_overlaps_for_bigwig(rna_1mb_hg38, agg = "mean")
export_signal_files(rna_1mb_hg38, bed_file = rna_1mb_bed_file, bw_file = rna_1mb_bw_file)

## CNV -------------------------------------------------------------------------
message("Preparing CNV signal track for ", sample_id, "...")
cnv <- read_tsv(cnv_file)

cnv_sample_col <- find_column(cnv, c("sample", "sampleID", "Sample", "ID"), "CNV sample")
cnv_chr_col <- find_column(cnv, c("chrom", "chr", "chromosome"), "CNV chromosome")
cnv_start_col <- find_column(cnv, c("start", "chromStart", "loc.start"), "CNV start")
cnv_end_col <- find_column(cnv, c("end", "chromEnd", "loc.end"), "CNV end")
cnv_score_col <- find_column(
  cnv,
  c("value", "score", "seg.mean", "Segment_Mean", "segmean", "log2", "log2ratio"),
  "CNV score"
)

cnv$sample_key <- sample_key(cnv[[cnv_sample_col]])
cnv_keep <- cnv$sample_key == target_sample_key
cnv <- cnv[cnv_keep, , drop = FALSE]

if (nrow(cnv) == 0L) {
  stop("No CNV segments found for sample ", sample_id, ".", call. = FALSE)
}

cnv_chr <- clean_chr(cnv[[cnv_chr_col]])
cnv_start <- as_numeric_clean(cnv[[cnv_start_col]])
cnv_end <- as_numeric_clean(cnv[[cnv_end_col]])
cnv_score <- as_numeric_clean(cnv[[cnv_score_col]])

cnv_keep <- cnv_chr %in% standard_chrs &
  is.finite(cnv_start) & is.finite(cnv_end) &
  cnv_start < cnv_end & is.finite(cnv_score)

cnv_hg19 <- make_granges_hg19(
  chr = cnv_chr[cnv_keep],
  start = cnv_start[cnv_keep],
  end = cnv_end[cnv_keep],
  score = cnv_score[cnv_keep],
  name = paste0("CNV_", seq_len(sum(cnv_keep)))
)
cnv_hg38 <- lift_all(cnv_hg19, chain)
cnv_hg38 <- resolve_overlaps_for_bigwig(cnv_hg38, agg = "max")
export_signal_files(cnv_hg38, bed_file = cnv_bed_file, bw_file = cnv_bw_file)

## Methylation 450k ------------------------------------------------------------
message("Preparing methylation 450k signal tracks for ", sample_id, "...")
methyl_sample <- read_xena_sample_vector(methyl_file, sample_id, "beta")
names(methyl_sample) <- c("id", "score")

probemap_450k <- read_probemap(probemap_450k_file, "probemap_450k")
methyl_merged <- merge(methyl_sample, probemap_450k, by = "id", all = FALSE)

if (nrow(methyl_merged) == 0L) {
  stop("Methylation probe rows did not match probemap_450k.tsv probe IDs.", call. = FALSE)
}

methyl_brca1 <- filter_region(
  methyl_merged,
  chr = brca1_chr,
  start = brca1_start_hg19,
  end = brca1_end_hg19
)

if (nrow(methyl_brca1) == 0L) {
  stop("No methylation 450k probes found in the BRCA1 hg19 plus/minus 500 kb region.", call. = FALSE)
}

methyl_brca1_hg19 <- make_granges_hg19(
  chr = methyl_brca1$chr,
  start = methyl_brca1$start,
  end = methyl_brca1$end,
  score = methyl_brca1$score,
  name = methyl_brca1$id
)
methyl_brca1_hg38 <- lift_first(methyl_brca1_hg19, chain)
methyl_brca1_hg38 <- resolve_overlaps_for_bigwig(methyl_brca1_hg38, agg = "mean")
export_signal_files(methyl_brca1_hg38, bw_file = methyl_brca1_bw_file)

methyl_hg19 <- make_granges_hg19(
  chr = methyl_merged$chr,
  start = methyl_merged$start,
  end = methyl_merged$end,
  score = methyl_merged$score,
  name = methyl_merged$id
)
methyl_hg38 <- lift_first(methyl_hg19, chain)
methyl_1mb_hg38 <- make_binned_signal(methyl_hg38, bin_size = 1e6, fun = "mean", min_n = 1)
methyl_1mb_hg38 <- resolve_overlaps_for_bigwig(methyl_1mb_hg38, agg = "mean")
export_signal_files(methyl_1mb_hg38, bed_file = methyl_1mb_bed_file, bw_file = methyl_1mb_bw_file)

## Mutation --------------------------------------------------------------------
message("Preparing mutation BED track for ", sample_id, "...")
mut <- read_tsv(mutation_file)

mut_sample_col <- find_column(
  mut,
  c("Tumor_Sample_Barcode", "sample", "sampleID", "Sample"),
  "mutation sample"
)
mut_chr_col <- find_column(mut, c("Chromosome", "chrom", "chr"), "mutation chromosome")
mut_start_col <- find_column(mut, c("Start_Position", "start", "chromStart"), "mutation start")
mut_end_col <- find_column(mut, c("End_Position", "end", "chromEnd"), "mutation end")
mut_effect_col <- find_column(
  mut,
  c("Variant_Classification", "effect", "variant_classification"),
  "mutation functional effect"
)
mut_gene_col <- find_column(
  mut,
  c("Hugo_Symbol", "gene", "Gene", "gene_symbol"),
  "mutation gene",
  required = FALSE
)
mut_protein_col <- find_column(
  mut,
  c("HGVSp_Short", "Protein_Change", "protein_change", "HGVSp"),
  "mutation protein change",
  required = FALSE
)

functional_effects <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
  "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del", "In_Frame_Ins",
  "Translation_Start_Site", "Nonstop_Mutation"
)

mut$sample_key <- sample_key(mut[[mut_sample_col]])
mut_chr <- clean_chr(mut[[mut_chr_col]])
mut_start <- as_numeric_clean(mut[[mut_start_col]])
mut_end <- as_numeric_clean(mut[[mut_end_col]])
mut_effect <- as.character(mut[[mut_effect_col]])

mut_keep <- mut$sample_key == target_sample_key &
  mut_effect %in% functional_effects &
  mut_chr == brca1_chr &
  is.finite(mut_start) & is.finite(mut_end) &
  mut_start <= brca1_end_hg19 & mut_end >= brca1_start_hg19

mut <- mut[mut_keep, , drop = FALSE]

if (nrow(mut) == 0L) {
  stop(
    "No functional mutations found for ", sample_id,
    " in the BRCA1 hg19 plus/minus 500 kb region.",
    call. = FALSE
  )
}

mut_chr <- clean_chr(mut[[mut_chr_col]])
mut_start <- as_numeric_clean(mut[[mut_start_col]])
mut_end <- as_numeric_clean(mut[[mut_end_col]])
mut_end[!is.finite(mut_end) | mut_end < mut_start] <- mut_start[!is.finite(mut_end) | mut_end < mut_start]

mut_gene <- if (!is.null(mut_gene_col)) as.character(mut[[mut_gene_col]]) else rep("Mutation", nrow(mut))
mut_protein <- if (!is.null(mut_protein_col)) as.character(mut[[mut_protein_col]]) else rep("", nrow(mut))
mut_effect <- as.character(mut[[mut_effect_col]])
mut_name <- trimws(paste(mut_gene, mut_protein, mut_effect))

if (!any(grepl("BRCA1", mut_gene, ignore.case = FALSE) & grepl("Q934\\*", mut_protein))) {
  warning(
    "BRCA1 p.Q934* was not detected in the mutation annotations after filtering; ",
    "check MC3 protein-change column naming if this is unexpected.",
    call. = FALSE
  )
}

mut_hg19 <- make_granges_hg19(
  chr = mut_chr,
  start = mut_start,
  end = mut_end,
  score = rep(1, nrow(mut)),
  name = mut_name,
  extra = list(
    gene = mut_gene,
    protein_change = mut_protein,
    variant_classification = mut_effect
  )
)
mut_hg38 <- lift_first(mut_hg19, chain)
mut_hg38 <- assign_hg38_seqinfo(mut_hg38, standard_chrs)
write_bed_from_granges(mut_hg38, mutation_bed_file, mutation = TRUE)

## Package object --------------------------------------------------------------
rel_path <- function(path) {
  root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  full <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(root, "/")

  if (startsWith(full, prefix)) {
    substring(full, nchar(prefix) + 1L)
  } else {
    full
  }
}

brca_omics_layer_tracks <- list(
  sample_id = sample_id,
  genome_build = "hg38",
  brca1_region_hg38 = brca1_region_hg38_for_plot,
  track_files = list(
    rna_bw = rel_path(rna_bw_file),
    rna_1Mb_bw = rel_path(rna_1mb_bw_file),
    rna_1Mb_bed = rel_path(rna_1mb_bed_file),
    cnv_bw = rel_path(cnv_bw_file),
    cnv_bed = rel_path(cnv_bed_file),
    methyl_brca1_bw = rel_path(methyl_brca1_bw_file),
    methyl_1Mb_bw = rel_path(methyl_1mb_bw_file),
    methyl_1Mb_bed = rel_path(methyl_1mb_bed_file),
    mutation_bed = rel_path(mutation_bed_file)
  ),
  track_labels = c(
    rna = "RNA-seq",
    cnv = "CNV",
    methylation = "Methylation 450k",
    mutations = "Mutations"
  ),
  notes = c(
    sample_provenance = paste0(
      sample_id,
      " selected because it is present in RNA, CNV, methylation 450k, and mutation layers."
    ),
    brca1_mutation = "Sample carries a BRCA1 nonsense mutation p.Q934* in MC3.",
    liftover = "Input genomic coordinates are hg19 and were lifted to hg38 with hg19ToHg38.over.chain.",
    file_generation_date = as.character(Sys.Date())
  )
)

usethis::use_data(brca_omics_layer_tracks, overwrite = TRUE)

message("Saved brca_omics_layer_tracks to data/brca_omics_layer_tracks.rda")
message("Wrote sample-specific omics-layer files to inst/extdata/brca/omics_layers/")
