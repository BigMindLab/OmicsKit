########################
# Function nice_circos #
########################

#' Required cytoband columns for circos plots
#'
#' @keywords internal
#' @noRd
.circos_required_cytoband_cols <- c("chr", "start", "end", "name", "gieStain")


#' Default qualitative colors for circos tracks
#'
#' @keywords internal
#' @noRd
.circos_default_colors <- c(
  "#7F77DD", "#1D9E75", "#EF9F27", "#D85A30", "#60A5FA",
  "#A855F7", "#14B8A6", "#F97316", "#64748B", "#84CC16",
  "#EC4899", "#06B6D4"
)


#' Safely coerce values to numeric
#'
#' @param x Vector-like object.
#'
#' @return Numeric vector.
#'
#' @keywords internal
#' @noRd
.circos_as_numeric <- function(x) {
  nm <- names(x)

  if (is.numeric(x)) {
    out <- as.numeric(x)
  } else {
    out <- suppressWarnings(as.numeric(as.character(x)))
  }

  names(out) <- nm
  out
}


#' Normalize column names for fuzzy matching
#'
#' @param x Character vector of names.
#'
#' @return Lower-case character vector with non-alphanumeric characters removed.
#'
#' @keywords internal
#' @noRd
.circos_clean_names <- function(x) {
  gsub("[^a-z0-9]", "", tolower(x))
}


#' Find a column index by aliases
#'
#' @param nms Character vector of observed names.
#' @param aliases Character vector of accepted aliases.
#'
#' @return Integer index of the first matching alias, or `NA_integer_`.
#'
#' @keywords internal
#' @noRd
.circos_find_alias <- function(nms, aliases) {
  nms_clean <- .circos_clean_names(nms)
  aliases_clean <- .circos_clean_names(aliases)

  idx <- which(nms_clean %in% aliases_clean)

  if (length(idx) == 0) {
    return(NA_integer_)
  }

  idx[1]
}


#' Ensure data frame column names exist and are unique
#'
#' @param x A data frame.
#'
#' @return A data frame with non-empty unique names.
#'
#' @keywords internal
#' @noRd
.circos_make_unique_names <- function(x) {
  if (is.null(names(x))) {
    names(x) <- paste0("V", seq_len(ncol(x)))
  }

  empty <- is.na(names(x)) | !nzchar(names(x))
  names(x)[empty] <- paste0("V", which(empty))
  names(x) <- make.unique(names(x), sep = "_")

  x
}


#' Construct a circos genome descriptor
#'
#' @param genome_build Character scalar. Genome build name.
#' @param source Character scalar describing genome source.
#' @param cytoband Cytoband data frame.
#' @param chromosomes Character vector of chromosomes in plotting order.
#' @param plot_type Character vector passed to `circlize`.
#' @param chr_lengths Optional named numeric vector of chromosome lengths.
#'
#' @return Object of class `"circos_genome"`.
#'
#' @keywords internal
#' @noRd
.make_circos_genome <- function(genome_build, source, cytoband, chromosomes,
                                plot_type, chr_lengths = NULL) {
  structure(
    list(
      genome_build = genome_build,
      source       = source,
      cytoband     = cytoband,
      chromosomes  = chromosomes,
      plot_type    = plot_type,
      chr_lengths  = chr_lengths
    ),
    class = "circos_genome"
  )
}


#' Validate and standardize cytoband data for circos plots
#'
#' Accepts UCSC/circlize-like cytoband data and standardizes it to the column
#' names required by `circlize::circos.initializeWithIdeogram()`.
#'
#' @param cytoband Data frame containing cytoband information.
#' @param arg Character scalar used in error messages.
#'
#' @return A data frame with columns `chr`, `start`, `end`, `name`,
#'   and `gieStain`.
#'
#' @section Notes:
#' The function accepts both explicitly named cytoband columns and unnamed
#' five-column UCSC-style cytobands in the order:
#' `chr`, `start`, `end`, `name`, `gieStain`.
#'
#' @keywords internal
#' @noRd
.validate_circos_cytoband <- function(cytoband, arg = "cytoband") {
  required <- .circos_required_cytoband_cols

  if (!is.data.frame(cytoband)) {
    stop(
      "`", arg, "` must be a data frame with columns: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }

  cytoband <- .circos_make_unique_names(cytoband)

  aliases <- list(
    chr      = c("chr", "chrom", "chromosome", "seqnames", "#chrom"),
    start    = c("start", "chromStart", "chrom_start", "chromStart0"),
    end      = c("end", "chromEnd", "chrom_end", "chromEnd1"),
    name     = c("name", "band", "cytoband", "id"),
    gieStain = c("gieStain", "gie_stain", "stain", "gie")
  )

  for (target in names(aliases)) {
    if (!target %in% names(cytoband)) {
      idx <- .circos_find_alias(names(cytoband), aliases[[target]])

      if (!is.na(idx)) {
        names(cytoband)[idx] <- target
      }
    }
  }

  if (!all(required %in% names(cytoband))) {
    if (ncol(cytoband) >= 5) {
      cytoband <- cytoband[, seq_len(5), drop = FALSE]
      names(cytoband) <- required
    } else {
      if (identical(arg, "custom_cytoband")) {
        stop(
          "custom_cytoband must have columns: chr, start, end, name, ",
          "gieStain or at least 5 columns in UCSC cytoband order.",
          call. = FALSE
        )
      }

      stop(
        "`", arg, "` must have columns: ",
        paste(required, collapse = ", "),
        " or at least 5 columns in UCSC cytoband order.",
        call. = FALSE
      )
    }
  } else {
    cytoband <- cytoband[, required, drop = FALSE]
  }

  cytoband$chr      <- as.character(cytoband$chr)
  cytoband$start    <- .circos_as_numeric(cytoband$start)
  cytoband$end      <- .circos_as_numeric(cytoband$end)
  cytoband$name     <- as.character(cytoband$name)
  cytoband$gieStain <- as.character(cytoband$gieStain)

  bad <- is.na(cytoband$chr) |
    !is.finite(cytoband$start) |
    !is.finite(cytoband$end) |
    cytoband$start > cytoband$end

  if (any(bad)) {
    stop(
      "`", arg, "` has invalid chr/start/end values.",
      call. = FALSE
    )
  }

  cytoband
}


#' Filter cytoband data by chromosome order
#'
#' @param cytoband Cytoband data frame.
#' @param chromosome_index Optional character vector of chromosomes to retain.
#' @param genome_build Character scalar used in error messages.
#'
#' @return A list with `cytoband` and `chromosomes`.
#'
#' @keywords internal
#' @noRd
.filter_circos_cytoband <- function(cytoband, chromosome_index,
                                    genome_build = "genome") {
  cytoband <- .validate_circos_cytoband(cytoband)
  cytoband$chr <- as.character(cytoband$chr)

  available <- unique(cytoband$chr)

  if (is.null(chromosome_index)) {
    chromosomes <- available

    cytoband <- cytoband[
      order(
        match(cytoband$chr, chromosomes),
        cytoband$start,
        cytoband$end
      ),
      ,
      drop = FALSE
    ]

    return(list(cytoband = cytoband, chromosomes = chromosomes))
  }

  chromosome_index <- as.character(chromosome_index)

  if (anyDuplicated(chromosome_index)) {
    stop("`chromosome_index` contains duplicated chromosomes.", call. = FALSE)
  }

  invalid <- setdiff(chromosome_index, available)

  if (length(invalid) > 0) {
    stop(
      "`chromosome_index` contains chromosomes not found in ",
      genome_build,
      ": ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }

  cytoband <- cytoband[cytoband$chr %in% chromosome_index, , drop = FALSE]
  cytoband$chr <- factor(cytoband$chr, levels = chromosome_index)

  cytoband <- cytoband[
    order(cytoband$chr, cytoband$start, cytoband$end),
    ,
    drop = FALSE
  ]

  cytoband$chr <- as.character(cytoband$chr)

  list(cytoband = cytoband, chromosomes = chromosome_index)
}


#' Create a synthetic cytoband from chromosome lengths
#'
#' @param chromosome_lengths Named numeric vector of chromosome lengths.
#'
#' @return A cytoband-like data frame with one `gneg` band per chromosome.
#'
#' @keywords internal
#' @noRd
.synthetic_cytoband <- function(chromosome_lengths) {
  if (is.null(names(chromosome_lengths)) ||
      any(!nzchar(names(chromosome_lengths)))) {
    stop("`chromosome_lengths` must be a named numeric vector.", call. = FALSE)
  }

  chromosome_lengths <- .circos_as_numeric(chromosome_lengths)

  if (any(!is.finite(chromosome_lengths)) || any(chromosome_lengths <= 0)) {
    stop(
      "`chromosome_lengths` must contain positive finite lengths.",
      call. = FALSE
    )
  }

  chromosomes <- names(chromosome_lengths)

  data.frame(
    chr      = chromosomes,
    start    = rep(0, length(chromosomes)),
    end      = as.numeric(chromosome_lengths),
    name     = chromosomes,
    gieStain = rep("gneg", length(chromosomes)),
    stringsAsFactors = FALSE
  )
}


#' Infer chromosome lengths from cytoband data
#'
#' @param cytoband Standardized cytoband data frame.
#'
#' @return Named numeric vector of chromosome lengths.
#'
#' @keywords internal
#' @noRd
.circos_chr_lengths_from_cytoband <- function(cytoband) {
  cytoband <- .validate_circos_cytoband(cytoband)

  x <- tapply(cytoband$end, cytoband$chr, max, na.rm = TRUE)
  stats::setNames(as.numeric(x), names(x))
}


#' Resolve cytoband and chromosome order for circos plots
#'
#' Builds a reusable genome descriptor for [nice_circos()]. For `ce10`, `ce11`,
#' and `dm6`, chromosome sizes are provided internally, so no UCSC download is
#' needed. Other genome builds are resolved with
#' `circlize::read.cytoband(species = genome_build)`.
#'
#' This function is tolerant of different cytoband column names. It standardizes
#' cytoband data to `chr`, `start`, `end`, `name`, and `gieStain` before
#' plotting.
#'
#' @param genome_build Character string. Genome build name, such as `"ce10"`,
#'   `"ce11"`, `"dm6"`, `"hg38"`, `"hg19"`, or `"mm10"`.
#' @param chromosome_index Optional character vector of chromosomes to keep and
#'   order in the circos plot. If `NULL`, all available chromosomes are used.
#'
#' @return An object of class `"circos_genome"` with elements:
#' \describe{
#'   \item{`genome_build`}{Requested genome build.}
#'   \item{`source`}{Source of the genome information: `"built-in"`,
#'   `"UCSC"`, or `"UCSC-chr.len"`.}
#'   \item{`cytoband`}{Data frame with columns `chr`, `start`, `end`,
#'   `name`, and `gieStain`.}
#'   \item{`chromosomes`}{Filtered chromosome vector in plotting order.}
#'   \item{`plot_type`}{Default `circlize` plot type vector.}
#'   \item{`chr_lengths`}{Named numeric vector of chromosome lengths.}
#' }
#'
#' @section Notes:
#' For `ce10`, `ce11`, and `dm6`, synthetic cytobands are generated from
#' built-in chromosome lengths. These do not display real cytogenetic bands;
#' they provide chromosome axes and labels.
#'
#' For UCSC-supported builds, the cytoband is obtained through
#' `circlize::read.cytoband()`. If only chromosome lengths are available,
#' synthetic cytobands are generated and the default plot type omits the
#' ideogram.
#'
#' @references
#' Gu Z, Gu L, Eils R, Schlesner M, Brors B. circlize implements and enhances
#' circular visualization in R. Bioinformatics. 2014;30(19):2811-2812.
#'
#' @seealso [nice_circos()], [is_circos_genome()]
#'
#' @examples
#' ce10 <- resolve_genome("ce10")
#' ce10
#'
#' resolve_genome("dm6", chromosome_index = c("chrX", "chr2L"))
#'
#' \dontrun{
#' hg38 <- resolve_genome("hg38", chromosome_index = "chr17")
#' print(hg38)
#' head(hg38$cytoband)
#' }
#'
#' @export
resolve_genome <- function(genome_build, chromosome_index = NULL) {
  if (!is.character(genome_build) || length(genome_build) != 1 ||
      is.na(genome_build) || !nzchar(genome_build)) {
    stop(
      "`genome_build` must be a single non-empty character string.",
      call. = FALSE
    )
  }

  genome_key <- tolower(genome_build)

  if (genome_key %in% c("ce10", "ce11")) {
    lengths <- c(
      chrI   = 15072434,
      chrII  = 15279421,
      chrIII = 13783801,
      chrIV  = 17493829,
      chrV   = 20924180,
      chrX   = 17718942
    )

    cytoband <- .synthetic_cytoband(lengths)

    filtered <- .filter_circos_cytoband(
      cytoband,
      chromosome_index,
      genome_build = genome_build
    )

    return(.make_circos_genome(
      genome_build = genome_build,
      source       = "built-in",
      cytoband     = filtered$cytoband,
      chromosomes  = filtered$chromosomes,
      plot_type    = c("axis", "labels"),
      chr_lengths  = lengths[filtered$chromosomes]
    ))
  }

  if (identical(genome_key, "dm6")) {
    lengths <- c(
      chr2L = 23513712,
      chr2R = 25286936,
      chr3L = 28110227,
      chr3R = 32079331,
      chr4  = 1348131,
      chrX  = 23542271
    )

    cytoband <- .synthetic_cytoband(lengths)

    filtered <- .filter_circos_cytoband(
      cytoband,
      chromosome_index,
      genome_build = genome_build
    )

    return(.make_circos_genome(
      genome_build = genome_build,
      source       = "built-in",
      cytoband     = filtered$cytoband,
      chromosomes  = filtered$chromosomes,
      plot_type    = c("axis", "labels"),
      chr_lengths  = lengths[filtered$chromosomes]
    ))
  }

  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop(
      "Package \"circlize\" must be installed to resolve UCSC genome builds.",
      call. = FALSE
    )
  }

  cytoband_result <- tryCatch(
    circlize::read.cytoband(species = genome_build),
    error = function(e) {
      stop(
        "Could not resolve genome build `", genome_build,
        "` with circlize::read.cytoband(). Original error: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  used_chr_len_only <- FALSE

  if (is.data.frame(cytoband_result)) {
    cytoband <- cytoband_result
    source <- "UCSC"
  } else if (!is.null(cytoband_result$df) &&
             is.data.frame(cytoband_result$df) &&
             nrow(cytoband_result$df) > 0) {
    cytoband <- cytoband_result$df
    source <- "UCSC"
  } else if (!is.null(cytoband_result$chr.len)) {
    cytoband <- .synthetic_cytoband(cytoband_result$chr.len)
    source <- "UCSC-chr.len"
    used_chr_len_only <- TRUE
  } else {
    stop(
      "`circlize::read.cytoband()` returned an unsupported object for `",
      genome_build,
      "`.",
      call. = FALSE
    )
  }

  cytoband <- .validate_circos_cytoband(cytoband)

  filtered <- .filter_circos_cytoband(
    cytoband,
    chromosome_index,
    genome_build = genome_build
  )

  chr_lengths <- .circos_chr_lengths_from_cytoband(filtered$cytoband)
  chr_lengths <- chr_lengths[filtered$chromosomes]

  .make_circos_genome(
    genome_build = genome_build,
    source       = source,
    cytoband     = filtered$cytoband,
    chromosomes  = filtered$chromosomes,
    plot_type    = if (used_chr_len_only) {
      c("axis", "labels")
    } else {
      c("ideogram", "axis", "labels")
    },
    chr_lengths  = chr_lengths
  )
}


#' Test whether an object is a circos genome descriptor
#'
#' @param x An R object.
#'
#' @return Logical scalar. `TRUE` if `x` inherits from class
#'   `"circos_genome"`, otherwise `FALSE`.
#'
#' @seealso [resolve_genome()], [nice_circos()]
#'
#' @examples
#' ce10 <- resolve_genome("ce10")
#' is_circos_genome(ce10)
#' is_circos_genome(data.frame())
#'
#' @export
is_circos_genome <- function(x) {
  inherits(x, "circos_genome")
}


#' Print a circos genome descriptor
#'
#' @param x A `"circos_genome"` object returned by [resolve_genome()].
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @seealso [resolve_genome()]
#'
#' @examples
#' ce10 <- resolve_genome("ce10")
#' print(ce10)
#'
#' @method print circos_genome
#' @export
print.circos_genome <- function(x, ...) {
  cat("Genome: ", x$genome_build, " (", x$source, ")\n", sep = "")

  cat(
    "Chromosomes: ",
    paste(x$chromosomes, collapse = ", "),
    " (",
    length(x$chromosomes),
    ")\n",
    sep = ""
  )

  cat("Plot type: ", paste(x$plot_type, collapse = ", "), "\n", sep = "")

  invisible(x)
}


#' Standardize track column names
#'
#' @param track_data Data frame-like object.
#' @param track_type Character scalar. Track type.
#' @param track_name Character scalar. Track name used in warnings/errors.
#'
#' @return Data frame with standardized names when possible.
#'
#' @keywords internal
#' @noRd
.standardize_circos_track_columns <- function(track_data, track_type,
                                              track_name = "track") {
  track_data <- as.data.frame(track_data)
  track_data <- .circos_make_unique_names(track_data)

  aliases <- list(
    chr = c("chr", "chrom", "chromosome", "seqnames", "seqname", "#chrom"),
    start = c("start", "chromStart", "chrom_start", "pos", "position"),
    end = c("end", "chromEnd", "chrom_end", "stop"),
    value = c(
      "value", "score", "expr", "expression", "beta", "cnv",
      "log2", "log2ratio", "segmean", "segmentmean", "DNA_VAF",
      "vaf", "af"
    ),
    name = c(
      "name", "id", "gene", "symbol", "label", "probe",
      "probe_id", "probeid", "Amino_Acid_Change"
    )
  )

  for (target in names(aliases)) {
    if (!target %in% names(track_data)) {
      idx <- .circos_find_alias(names(track_data), aliases[[target]])

      if (!is.na(idx)) {
        names(track_data)[idx] <- target
      }
    }
  }

  if (!all(c("chr", "start", "end") %in% names(track_data)) &&
      ncol(track_data) >= 3) {
    names(track_data)[1:3] <- c("chr", "start", "end")
  }

  if (!"name" %in% names(track_data) &&
      ncol(track_data) >= 4 &&
      identical(track_type, "mutation")) {
    names(track_data)[4] <- "name"
  }

  if (!"value" %in% names(track_data) && ncol(track_data) >= 5) {
    candidate <- .circos_as_numeric(track_data[[5]])

    if (any(is.finite(candidate), na.rm = TRUE)) {
      names(track_data)[5] <- "value"
    }
  }

  if (!"value" %in% names(track_data) &&
      ncol(track_data) >= 4 &&
      !identical(track_type, "mutation")) {
    candidate <- .circos_as_numeric(track_data[[4]])

    if (any(is.finite(candidate), na.rm = TRUE)) {
      names(track_data)[4] <- "value"
    }
  }

  track_data
}


#' Normalize one circos track
#'
#' Converts file paths, `GRanges`, and data frames into a standardized track
#' representation used internally by [nice_circos()].
#'
#' @param track Data frame, `GRanges`, file path, or config list.
#' @param track_name Character scalar. Track name.
#' @param track_type Optional character scalar. Track type.
#'
#' @return A list with elements `data`, `type`, `color`, and `height`.
#'
#' @keywords internal
#' @noRd
.normalize_circos_track <- function(track, track_name, track_type = NULL) {
  track_data <- track
  config_type <- NULL
  config_color <- NULL
  config_height <- NULL
  config_gradient <- NULL

  is_config_list <- is.list(track) &&
    !is.data.frame(track) &&
    !inherits(track, "GRanges") &&
    any(c("data", "path", "file") %in% names(track))

  if (is_config_list) {
    if ("data" %in% names(track)) {
      track_data <- track$data
    } else if ("path" %in% names(track)) {
      track_data <- track$path
    } else {
      track_data <- track$file
    }

    config_type   <- track$type
    config_color  <- track$color
    config_height <- track$height
    config_gradient <- track$gradient
  }

  if (is.null(track_type)) {
    track_type <- config_type
  }

  if (is.null(track_type) || is.na(track_type) || !nzchar(track_type)) {
    lower_name <- tolower(track_name)

    track_type <- if (grepl("mut", lower_name)) {
      "mutation"
    } else if (grepl("cnv|copy", lower_name)) {
      "line"
    } else if (grepl("methyl|meth|beta", lower_name)) {
      "points"
    } else {
      "histogram"
    }
  }

  track_type <- tolower(as.character(track_type))

  allowed_types <- c("histogram", "scatter", "points", "line", "mutation")

  if (!track_type %in% allowed_types) {
    stop(
      "`track_types` values must be one of: ",
      paste(allowed_types, collapse = ", "),
      call. = FALSE
    )
  }

  if (identical(track_type, "scatter")) {
    track_type <- "points"
  }

  if (!is.null(config_height)) {
    config_height <- .circos_as_numeric(config_height)

    if (length(config_height) != 1 ||
        !is.finite(config_height) ||
        config_height <= 0) {
      stop(
        "`height` for track `", track_name,
        "` must be a positive numeric scalar.",
        call. = FALSE
      )
    }
  }

  if (!is.null(config_color)) {
    if (length(config_color) < 1 ||
        is.na(config_color[1]) ||
        !nzchar(as.character(config_color[1]))) {
      stop(
        "`color` for track `", track_name,
        "` must be a non-empty color value.",
        call. = FALSE
      )
    }

    if (length(config_color) > 1) {
      warning(
        "Track `", track_name,
        "` has multiple colors; using the first one.",
        call. = FALSE
      )
      config_color <- config_color[1]
    }

    config_color <- as.character(config_color)
  }

  if (is.character(track_data) && length(track_data) == 1) {
    if (!file.exists(track_data)) {
      stop("Track file not found: ", track_data, call. = FALSE)
    }

    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
      stop(
        "Package \"rtracklayer\" must be installed to import track files.",
        call. = FALSE
      )
    }

    ext <- tolower(tools::file_ext(track_data))

    if (!ext %in% c("bed", "bw", "bigwig")) {
      stop(
        "Unsupported file format: .", ext,
        "\nSupported: .bed, .bw, .bigwig",
        call. = FALSE
      )
    }

    track_data <- tryCatch(
      as.data.frame(rtracklayer::import(track_data)),
      error = function(e) {
        stop(
          "Could not import track file `", track_data,
          "`. Original error: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  } else if (inherits(track_data, "GRanges")) {
    track_data <- as.data.frame(track_data)
  } else {
    track_data <- as.data.frame(track_data)
  }

  track_data <- .standardize_circos_track_columns(
    track_data,
    track_type = track_type,
    track_name = track_name
  )

  required_cols <- c("chr", "start", "end")

  if (!all(required_cols %in% names(track_data))) {
    stop(
      "Track `", track_name,
      "` must have columns: chr, start, end. ",
      "Optional columns: value, name.",
      call. = FALSE
    )
  }

  had_value <- "value" %in% names(track_data)

  if (!had_value) {
    track_data$value <- 1

    if (!identical(track_type, "mutation")) {
      warning(
        "Track `", track_name,
        "` has no `value` or `score` column; using value = 1.",
        call. = FALSE
      )
    }
  }

  if (!"name" %in% names(track_data)) {
    track_data$name <- NA_character_
  }

  track_data <- track_data[, c("chr", "start", "end", "value", "name"),
                           drop = FALSE]

  track_data$chr   <- as.character(track_data$chr)
  track_data$start <- .circos_as_numeric(track_data$start)
  track_data$end   <- .circos_as_numeric(track_data$end)
  track_data$value <- .circos_as_numeric(track_data$value)
  track_data$name  <- as.character(track_data$name)

  bad_coords <- is.na(track_data$chr) |
    !is.finite(track_data$start) |
    !is.finite(track_data$end) |
    track_data$start > track_data$end

  if (any(bad_coords)) {
    stop(
      "Track `", track_name,
      "` contains invalid chr/start/end values.",
      call. = FALSE
    )
  }

  if (identical(track_type, "mutation")) {
    track_data$value[!is.finite(track_data$value)] <- 1
  } else {
    bad_value <- !is.finite(track_data$value)

    if (any(bad_value)) {
      warning(
        "Track `", track_name,
        "` contains ", sum(bad_value),
        " rows with non-finite values; dropping them.",
        call. = FALSE
      )

      track_data <- track_data[!bad_value, , drop = FALSE]
    }
  }

  if (nrow(track_data) > 0) {
    track_data <- track_data[
      order(track_data$chr, track_data$start, track_data$end),
      ,
      drop = FALSE
    ]
  }
  if (is.null(config_gradient)) {
    config_gradient <- TRUE
  }

  config_gradient <- isTRUE(config_gradient)

  list(
    data   = track_data,
    type   = track_type,
    color  = config_color,
    height = config_height,
    gradient = config_gradient
  )
}


#' Normalize metadata rings for circos plots
#'
#' @param metadata Named vector, one-row data frame, or genomic data frame with
#'   `chr`, `start`, `end`, and `group` columns.
#'
#' @return Normalized metadata data frame, or `NULL`.
#'
#' @keywords internal
#' @noRd
.normalize_circos_metadata <- function(metadata) {
  if (is.null(metadata)) {
    return(NULL)
  }

  if (is.atomic(metadata) && !is.null(names(metadata))) {
    metadata <- data.frame(
      variable = names(metadata),
      group = as.character(metadata),
      stringsAsFactors = FALSE
    )

    return(metadata)
  }

  metadata <- as.data.frame(metadata)
  metadata <- .circos_make_unique_names(metadata)

  aliases <- list(
    chr = c("chr", "chrom", "chromosome", "seqnames", "seqname", "#chrom"),
    start = c("start", "chromStart", "chrom_start", "pos", "position"),
    end = c("end", "chromEnd", "chrom_end", "stop"),
    group = c("group", "class", "category", "value", "label"),
    variable = c("variable", "track", "feature", "annotation")
  )

  for (target in names(aliases)) {
    if (!target %in% names(metadata)) {
      idx <- .circos_find_alias(names(metadata), aliases[[target]])

      if (!is.na(idx)) {
        names(metadata)[idx] <- target
      }
    }
  }

  if (all(c("chr", "start", "end", "group") %in% names(metadata))) {
    if (!"variable" %in% names(metadata)) {
      metadata$variable <- "metadata"
    }

    metadata$chr <- as.character(metadata$chr)
    metadata$start <- .circos_as_numeric(metadata$start)
    metadata$end <- .circos_as_numeric(metadata$end)
    metadata$group <- as.character(metadata$group)
    metadata$variable <- as.character(metadata$variable)

    bad <- is.na(metadata$chr) |
      !is.finite(metadata$start) |
      !is.finite(metadata$end) |
      metadata$start > metadata$end |
      is.na(metadata$group)

    if (any(bad)) {
      stop(
        "`metadata` contains invalid chr/start/end/group values.",
        call. = FALSE
      )
    }

    return(metadata)
  }

  if (nrow(metadata) != 1) {
    stop(
      "`metadata` must be a named vector, one-row data frame, or a data ",
      "frame with chr/start/end/group columns.",
      call. = FALSE
    )
  }

  data.frame(
    variable = names(metadata),
    group = as.character(metadata[1, , drop = TRUE]),
    stringsAsFactors = FALSE
  )
}


#' Infer chromosome end positions from normalized tracks
#'
#' @param normalized_tracks List of normalized tracks.
#'
#' @return Named numeric vector of chromosome end positions, or `NULL`.
#'
#' @keywords internal
#' @noRd
.circos_track_extents <- function(normalized_tracks) {
  if (length(normalized_tracks) == 0) {
    return(NULL)
  }

  dfs <- lapply(normalized_tracks, function(x) {
    x$data[, c("chr", "start", "end"), drop = FALSE]
  })

  dfs <- dfs[vapply(dfs, nrow, integer(1)) > 0]

  if (length(dfs) == 0) {
    return(NULL)
  }

  all_df <- do.call(rbind, dfs)

  x <- tapply(all_df$end, all_df$chr, max, na.rm = TRUE)
  stats::setNames(as.numeric(x), names(x))
}


#' Get chromosome ends for metadata rings
#'
#' @param chromosomes Character vector of chromosomes.
#' @param cytoband Optional cytoband data frame.
#' @param normalized_tracks List of normalized tracks.
#'
#' @return Named numeric vector of chromosome ends.
#'
#' @keywords internal
#' @noRd
.circos_chromosome_ends <- function(chromosomes, cytoband = NULL,
                                    normalized_tracks = list()) {
  if (!is.null(cytoband)) {
    cytoband <- .validate_circos_cytoband(cytoband)
    ends <- .circos_chr_lengths_from_cytoband(cytoband)
  } else {
    ends <- .circos_track_extents(normalized_tracks)
  }

  if (is.null(ends)) {
    stop(
      "Could not infer chromosome lengths for metadata rings.",
      call. = FALSE
    )
  }

  missing <- setdiff(chromosomes, names(ends))

  if (length(missing) > 0) {
    stop(
      "Could not infer chromosome lengths for: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  ends[chromosomes]
}


#' Compute y-axis limits for a circos track
#'
#' @param values Numeric vector of track values.
#' @param track_type Character scalar. Track type.
#'
#' @return Numeric vector of length two.
#'
#' @keywords internal
#' @noRd
.circos_get_ylim <- function(values, track_type) {
  values <- values[is.finite(values)]

  if (identical(track_type, "mutation")) {
    return(c(0, 1))
  }

  if (length(values) == 0) {
    return(c(0, 1))
  }

  rng <- if (identical(track_type, "histogram")) {
    range(c(0, values), na.rm = TRUE)
  } else {
    range(values, na.rm = TRUE)
  }

  if (!all(is.finite(rng))) {
    return(c(0, 1))
  }

  if (identical(rng[1], rng[2])) {
    pad <- if (identical(rng[1], 0)) 0.5 else abs(rng[1]) * 0.05

    if (!is.finite(pad) || pad == 0) {
      pad <- 0.5
    }

    rng <- rng + c(-pad, pad)
  }

  rng
}


#' Filter normalized tracks to selected chromosomes
#'
#' @param normalized_tracks List of normalized tracks.
#' @param chromosomes Character vector of selected chromosomes.
#'
#' @return List of normalized tracks.
#'
#' @keywords internal
#' @noRd
.circos_filter_tracks_to_chromosomes <- function(normalized_tracks,
                                                 chromosomes) {
  if (is.null(chromosomes) || length(chromosomes) == 0) {
    return(normalized_tracks)
  }

  out <- lapply(seq_along(normalized_tracks), function(i) {
    nm <- names(normalized_tracks)[i]
    track <- normalized_tracks[[i]]
    n_before <- nrow(track$data)

    track$data <- track$data[track$data$chr %in% chromosomes, , drop = FALSE]

    if (n_before > 0 && nrow(track$data) == 0) {
      warning(
        "Track `", nm,
        "` has no rows after applying `chromosome_index`.",
        call. = FALSE
      )
    }

    if (nrow(track$data) > 0) {
      track$data$chr <- factor(track$data$chr, levels = chromosomes)

      track$data <- track$data[
        order(track$data$chr, track$data$start, track$data$end),
        ,
        drop = FALSE
      ]

      track$data$chr <- as.character(track$data$chr)
    }

    track
  })

  names(out) <- names(normalized_tracks)

  out
}


#' Validate genomic link data
#'
#' @param link_data Data frame with paired genomic coordinates.
#' @param chromosomes Optional character vector of selected chromosomes.
#'
#' @return Validated link data frame, or `NULL`.
#'
#' @keywords internal
#' @noRd
.circos_validate_link_data <- function(link_data, chromosomes = NULL) {
  if (is.null(link_data)) {
    return(NULL)
  }

  link_data <- as.data.frame(link_data)
  required_link_cols <- c("chr1", "start1", "end1", "chr2", "start2", "end2")

  if (!all(required_link_cols %in% names(link_data))) {
    stop(
      "`link_data` must have columns: ",
      paste(required_link_cols, collapse = ", "),
      call. = FALSE
    )
  }

  link_data <- link_data[, required_link_cols, drop = FALSE]

  link_data$chr1 <- as.character(link_data$chr1)
  link_data$chr2 <- as.character(link_data$chr2)

  for (x in c("start1", "end1", "start2", "end2")) {
    link_data[[x]] <- .circos_as_numeric(link_data[[x]])
  }

  bad <- is.na(link_data$chr1) |
    is.na(link_data$chr2) |
    !is.finite(link_data$start1) |
    !is.finite(link_data$end1) |
    !is.finite(link_data$start2) |
    !is.finite(link_data$end2) |
    link_data$start1 > link_data$end1 |
    link_data$start2 > link_data$end2

  if (any(bad)) {
    stop("`link_data` contains invalid coordinates.", call. = FALSE)
  }

  if (!is.null(chromosomes)) {
    keep <- link_data$chr1 %in% chromosomes & link_data$chr2 %in% chromosomes

    if (any(!keep)) {
      warning(
        "Dropping ", sum(!keep),
        " link_data rows outside selected chromosomes.",
        call. = FALSE
      )
    }

    link_data <- link_data[keep, , drop = FALSE]
  }

  link_data
}


#' Plot genome-wide circos tracks for multi-omics data
#'
#' Builds a circular genome plot with optional data tracks, metadata rings, and
#' genomic links using `circlize`.
#'
#' Tracks must already use the selected genome build. For non-human builds
#' without UCSC cytobands, use [resolve_genome()] support for `ce10`, `ce11`,
#' or `dm6`, or provide `custom_cytoband`.
#'
#' Input tracks can be data frames, `GRanges`, file paths, or configuration
#' lists. File paths can point to `.bed`, `.bw`, or `.bigwig` files. Imported
#' `score` columns are automatically mapped to `value`.
#'
#' @param genome_build Character string. Genome build used to initialize the
#'   plot, such as `"hg38"`, `"hg19"`, `"mm10"`, `"ce10"`, `"ce11"`, or
#'   `"dm6"`.
#' @param data_tracks Named list of data frames, `GRanges`, file paths, or
#'   track config lists. Data frames must include `chr`, `start`, `end`, and
#'   optionally `value` and `name`.
#' @param track_types Optional character vector with one type per track.
#'   Supported values are `"histogram"`, `"line"`, `"points"`, `"scatter"`,
#'   and `"mutation"`. If `NULL`, types are inferred from track names.
#' @param metadata Optional named vector, one-row data frame, or genomic data
#'   frame with `chr`, `start`, `end`, and `group` columns for metadata rings.
#' @param meta_colors Optional named vector of colors for metadata groups.
#'   Names should correspond to values in `metadata$group`.
#' @param track_colors Optional named vector or list of colors keyed by track
#'   name.
#' @param ideogram Logical. If `TRUE`, initialize with chromosome ideogram,
#'   axis, and labels according to `plot_type`.
#' @param link_data Optional data frame with columns `chr1`, `start1`, `end1`,
#'   `chr2`, `start2`, and `end2`.
#' @param export_pdf Optional PDF file path. If `NULL`, draws on the active
#'   graphics device.
#' @param custom_cytoband Optional cytoband data frame with columns `chr`,
#'   `start`, `end`, `name`, and `gieStain`. Unnamed five-column UCSC-style
#'   cytobands are also accepted.
#' @param chromosome_index Optional character vector of chromosomes to keep and
#'   order.
#' @param plot_type Optional character vector passed as `plotType` to
#'   `circlize::circos.initializeWithIdeogram()`. Overrides the default from
#'   [resolve_genome()].
#' @param show_legend Logical. If `TRUE`, draw a simple legend using
#'   `ComplexHeatmap`.
#' @param show_labels Logical. If `TRUE`, label mutation tracks using the
#'   `name` column when present.
#' @param track_height Numeric scalar. Default height for data tracks.
#' @param metadata_height Numeric scalar. Height for metadata tracks.
#'
#' @return Invisibly returns a list with:
#' \describe{
#'   \item{`tracks`}{Normalized data tracks used for plotting.}
#'   \item{`metadata`}{Normalized metadata, or `NULL`.}
#'   \item{`track_settings`}{Track type, color, height, y-limits, and row count.}
#'   \item{`link_data`}{Validated link data, or `NULL`.}
#'   \item{`genome_build`}{Genome build used.}
#'   \item{`genome`}{`circos_genome` object returned by [resolve_genome()],
#'   or `NULL` if `custom_cytoband` was used.}
#'   \item{`chromosomes`}{Chromosomes shown in the plot.}
#'   \item{`plot_type`}{Final plot type passed to `circlize`.}
#'   \item{`cytoband`}{Cytoband used to initialize the plot.}
#'   \item{`pdf_opened`}{Logical indicating whether a PDF device was opened.}
#' }
#'
#' @section Track input:
#' Each data track may be:
#' \itemize{
#'   \item a data frame with `chr`, `start`, `end`, and optionally `value`,
#'   `name`;
#'   \item a `GRanges` object;
#'   \item a path to `.bed`, `.bw`, or `.bigwig`;
#'   \item a config list such as `list(path = "x.bw", type = "line",
#'   color = "steelblue", height = 0.1)`.
#' }
#'
#' @section Notes:
#' Coordinates are assumed to match `genome_build`. This function does not
#' perform liftOver.
#'
#' For sparse regional assays, such as methylation probes or mutation calls
#' restricted to a single gene window, a whole-genome circos plot may make the
#' signal difficult to see. In that case, use `chromosome_index` to focus on the
#' relevant chromosome.
#'
#' Mutation tracks are treated as discrete points. If no `value` column is
#' present, value `1` is used.
#'
#' @section Package requirements:
#' This function requires `circlize`. Importing `.bed`, `.bw`, or `.bigwig`
#' files requires `rtracklayer`. Drawing legends requires `ComplexHeatmap`.
#'
#' @references
#' Gu Z, Gu L, Eils R, Schlesner M, Brors B. circlize implements and enhances
#' circular visualization in R. Bioinformatics. 2014;30(19):2811-2812.
#'
#' @seealso [resolve_genome()], [is_circos_genome()]
#'
#' @examples
#' \dontrun{
#' tracks <- list(
#'   RNA = list(
#'     path = "rna_expr_hg38.bw",
#'     type = "histogram",
#'     height = 0.12
#'   ),
#'   CNV = list(
#'     path = "cnv_hg38.bw",
#'     type = "line",
#'     height = 0.10
#'   ),
#'   Methylation = list(
#'     path = "methyl_brca1_hg38.bw",
#'     type = "points",
#'     height = 0.10
#'   ),
#'   Mutation = list(
#'     path = "mutations_brca1_hg38.bed",
#'     type = "mutation",
#'     height = 0.08
#'   )
#' )
#'
#' nice_circos(
#'   genome_build = "hg38",
#'   data_tracks = tracks,
#'   chromosome_index = "chr17",
#'   plot_type = c("axis", "labels"),
#'   show_labels = TRUE,
#'   show_legend = TRUE,
#'   export_pdf = "TCGA_A1_A0SH_chr17_circos.pdf"
#' )
#' }
#'
#' \dontrun{
#' cnv <- data.frame(
#'   chr = c("chrI", "chrII"),
#'   start = c(1, 1),
#'   end = c(1e6, 1e6),
#'   value = c(0.3, -0.2)
#' )
#'
#' nice_circos(
#'   genome_build = "ce10",
#'   data_tracks = list(CNV = cnv),
#'   track_types = "line"
#' )
#' }
#'
#' @importFrom grDevices adjustcolor dev.off pdf
#' @importFrom grid gpar
#' @importFrom stats setNames
#' @importFrom tools file_ext
#'
#' @export
nice_circos <- function(
    genome_build  = "hg38",
    data_tracks   = list(),
    track_types   = NULL,
    metadata      = NULL,
    meta_colors   = NULL,
    track_colors  = NULL,
    ideogram      = TRUE,
    link_data     = NULL,
    export_pdf    = NULL,
    custom_cytoband = NULL,
    chromosome_index = NULL,
    plot_type = NULL,
    show_legend = FALSE,
    show_labels = FALSE,
    track_height = 0.10,
    metadata_height = 0.08
) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop(
      "Package \"circlize\" must be installed to use this function.",
      call. = FALSE
    )
  }

  track_height <- .circos_as_numeric(track_height)
  metadata_height <- .circos_as_numeric(metadata_height)

  if (length(track_height) != 1 ||
      !is.finite(track_height) ||
      track_height <= 0) {
    stop("`track_height` must be a positive numeric scalar.", call. = FALSE)
  }

  if (length(metadata_height) != 1 ||
      !is.finite(metadata_height) ||
      metadata_height <= 0) {
    stop("`metadata_height` must be a positive numeric scalar.",
         call. = FALSE)
  }

  if (!is.null(custom_cytoband)) {
    custom_cytoband <- .validate_circos_cytoband(
      custom_cytoband,
      arg = "custom_cytoband"
    )
  }

  if (!is.list(data_tracks) || is.data.frame(data_tracks)) {
    stop("`data_tracks` must be a list.", call. = FALSE)
  }

  n_tracks <- length(data_tracks)

  if (n_tracks > 0) {
    if (is.null(names(data_tracks))) {
      names(data_tracks) <- paste0("Track", seq_len(n_tracks))
    } else {
      unnamed <- is.na(names(data_tracks)) | !nzchar(names(data_tracks))
      names(data_tracks)[unnamed] <- paste0("Track", which(unnamed))
    }

    if (anyDuplicated(names(data_tracks))) {
      stop("`data_tracks` names must be unique.", call. = FALSE)
    }
  }

  if (!is.null(track_types)) {
    if (length(track_types) != n_tracks) {
      stop("`track_types` must have one value per data track.", call. = FALSE)
    }

    if (!is.null(names(track_types)) && any(nzchar(names(track_types)))) {
      missing_types <- setdiff(names(data_tracks), names(track_types))

      if (length(missing_types) > 0) {
        stop(
          "`track_types` is named but missing entries for: ",
          paste(missing_types, collapse = ", "),
          call. = FALSE
        )
      }

      track_types <- unname(track_types[names(data_tracks)])
    }
  }

  normalized_tracks <- vector("list", n_tracks)
  names(normalized_tracks) <- names(data_tracks)

  for (i in seq_along(data_tracks)) {
    type_i <- if (is.null(track_types)) NULL else track_types[[i]]

    normalized_tracks[[i]] <- .normalize_circos_track(
      data_tracks[[i]],
      track_name = names(data_tracks)[i],
      track_type = type_i
    )
  }

  normalized_metadata <- .normalize_circos_metadata(metadata)

  genome_info <- NULL
  chromosomes <- chromosome_index
  resolved_plot_type <- plot_type
  cytoband <- custom_cytoband

  if (isTRUE(ideogram)) {
    if (!is.null(custom_cytoband)) {
      filtered <- .filter_circos_cytoband(
        custom_cytoband,
        chromosome_index,
        genome_build = "custom_cytoband"
      )

      cytoband <- filtered$cytoband
      chromosomes <- filtered$chromosomes

      if (is.null(resolved_plot_type)) {
        resolved_plot_type <- c("ideogram", "axis", "labels")
      }
    } else {
      genome_info <- resolve_genome(genome_build, chromosome_index)

      if (!is_circos_genome(genome_info)) {
        stop(
          "`resolve_genome()` must return a circos_genome object.",
          call. = FALSE
        )
      }

      cytoband <- genome_info$cytoband
      chromosomes <- genome_info$chromosomes

      if (is.null(resolved_plot_type)) {
        resolved_plot_type <- genome_info$plot_type
      }
    }
  }

  if (is.null(chromosomes)) {
    chromosomes <- unique(unlist(lapply(normalized_tracks, function(x) {
      x$data$chr
    }), use.names = FALSE))
  }

  chromosomes <- as.character(chromosomes)

  if (length(chromosomes) == 0 && !isTRUE(ideogram)) {
    stop(
      "`data_tracks` with at least one non-empty track is required when ",
      "`ideogram = FALSE`.",
      call. = FALSE
    )
  }

  if (!is.null(chromosomes) && length(chromosomes) > 0) {
    normalized_tracks <- .circos_filter_tracks_to_chromosomes(
      normalized_tracks,
      chromosomes
    )
  }

  link_data <- .circos_validate_link_data(link_data, chromosomes)

  if (!is.null(normalized_metadata) &&
      all(c("chr", "start", "end", "group") %in% names(normalized_metadata)) &&
      length(chromosomes) > 0) {
    n_before <- nrow(normalized_metadata)

    normalized_metadata <- normalized_metadata[
      normalized_metadata$chr %in% chromosomes,
      ,
      drop = FALSE
    ]

    if (n_before > 0 && nrow(normalized_metadata) == 0) {
      warning(
        "`metadata` has no rows after applying `chromosome_index`.",
        call. = FALSE
      )
    }
  }

  pdf_opened <- FALSE

  if (!is.null(export_pdf)) {
    if (!is.character(export_pdf) || length(export_pdf) != 1 ||
        is.na(export_pdf) || !nzchar(export_pdf)) {
      stop("`export_pdf` must be a non-empty file path.", call. = FALSE)
    }

    export_dir <- dirname(export_pdf)

    if (!dir.exists(export_dir)) {
      stop(
        "Directory for `export_pdf` does not exist: ",
        export_dir,
        call. = FALSE
      )
    }

    grDevices::pdf(export_pdf, width = 8, height = 8)
    pdf_opened <- TRUE
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)

  if (isTRUE(ideogram)) {
    circlize::circos.initializeWithIdeogram(
      cytoband = cytoband,
      chromosome.index = chromosomes,
      plotType = resolved_plot_type
    )
  } else {
    init_dfs <- lapply(normalized_tracks, function(x) {
      x$data[, c("chr", "start", "end"), drop = FALSE]
    })

    init_dfs <- init_dfs[vapply(init_dfs, nrow, integer(1)) > 0]

    if (length(init_dfs) == 0) {
      stop(
        "`data_tracks` with at least one non-empty track is required when ",
        "`ideogram = FALSE`.",
        call. = FALSE
      )
    }

    init_df <- do.call(rbind, init_dfs)
    init_df$chr <- factor(init_df$chr, levels = chromosomes)

    init_df <- init_df[
      order(init_df$chr, init_df$start, init_df$end),
      ,
      drop = FALSE
    ]

    circlize::circos.genomicInitialize(init_df)
  }

  if (!is.null(normalized_metadata) && nrow(normalized_metadata) > 0) {
    if (all(c("chr", "start", "end", "group") %in% names(normalized_metadata))) {
      metadata_tracks <- split(
        normalized_metadata,
        normalized_metadata$variable
      )
    } else {
      chromosome_ends <- .circos_chromosome_ends(
        chromosomes = chromosomes,
        cytoband = cytoband,
        normalized_tracks = normalized_tracks
      )

      metadata_tracks <- lapply(seq_len(nrow(normalized_metadata)), function(i) {
        data.frame(
          chr = chromosomes,
          start = 0,
          end = as.numeric(chromosome_ends[chromosomes]),
          group = normalized_metadata$group[i],
          stringsAsFactors = FALSE
        )
      })

      names(metadata_tracks) <- normalized_metadata$variable
    }

    for (meta_name in names(metadata_tracks)) {
      meta_df <- metadata_tracks[[meta_name]]
      meta_df <- meta_df[meta_df$chr %in% chromosomes, , drop = FALSE]

      if (nrow(meta_df) == 0) {
        next
      }

      groups <- unique(as.character(meta_df$group))

      if (is.null(meta_colors)) {
        meta_palette <- stats::setNames(
          rep(.circos_default_colors, length.out = length(groups)),
          groups
        )
      } else {
        meta_palette <- meta_colors
      }

      circlize::circos.genomicTrack(
        meta_df[, c("chr", "start", "end", "group"), drop = FALSE],
        ylim = c(0, 1),
        panel.fun = function(region, value, ...) {
          group_values <- as.character(value[[1]])
          cols <- meta_palette[group_values]
          cols[is.na(cols)] <- "grey80"

          circlize::circos.genomicRect(
            region,
            ybottom = 0,
            ytop = 1,
            col = cols,
            border = NA,
            ...
          )
        },
        track.height = metadata_height
      )
    }
  }

  track_settings <- list()

  for (i in seq_along(normalized_tracks)) {
    nm <- names(normalized_tracks)[i]
    track <- normalized_tracks[[i]]
    bed <- track$data

    if (nrow(bed) == 0) {
      next
    }

    color <- if (!is.null(track$color)) {
      track$color
    } else if (!is.null(track_colors) && nm %in% names(track_colors)) {
      as.character(track_colors[[nm]])[1]
    } else {
      .circos_default_colors[
        ((i - 1) %% length(.circos_default_colors)) + 1
      ]
    }

    height <- if (!is.null(track$height)) {
      track$height
    } else {
      track_height
    }

    values <- bed$value
    ylim <- .circos_get_ylim(values, track$type)

    gradient <- if (!is.null(track$gradient)) {
      isTRUE(track$gradient)
      } else {
        TRUE
        }

    track_settings[[nm]] <- list(
      type = track$type,
      color = color,
      height = height,
      ylim = ylim,
      n = nrow(bed),
      gradient = gradient
    )

    local({
      bed_i <- bed
      type_i <- track$type
      color_i <- color
      ylim_i <- ylim
      values_i <- values
      show_labels_i <- show_labels
      height_i <- height
      gradient_i <- gradient

      circlize::circos.genomicTrack(
        bed_i,
        ylim = ylim_i,
        numeric.column = "value",
        panel.fun = function(region, value, ...) {
          v <- if ("value" %in% names(value)) {
            value[["value"]]
          } else {
            value[[1]]
          }

          if (identical(type_i, "histogram")) {
            finite_v <- values_i[is.finite(values_i)]

            if (isTRUE(gradient_i) && length(unique(finite_v)) > 1) {
              col_ramp <- circlize::colorRamp2(
                range(finite_v, na.rm = TRUE),
                c("#EEEDFE", color_i)
              )

              rect_cols <- col_ramp(v)
            } else {
              rect_cols <- rep(color_i, length(v))
            }

            circlize::circos.genomicRect(
              region,
              ybottom = pmin(0, v),
              ytop = pmax(0, v),
              col = rect_cols,
              border = NA,
              ...
            )
          } else if (identical(type_i, "line")) {
            circlize::circos.genomicLines(
              region,
              value,
              numeric.column = "value",
              col = color_i,
              lwd = 1,
              ...
            )
          } else if (identical(type_i, "points")) {
            circlize::circos.genomicPoints(
              region,
              value,
              numeric.column = "value",
              col = color_i,
              pch = 16,
              cex = 0.4,
              ...
            )
          } else if (identical(type_i, "mutation")) {
            circlize::circos.genomicPoints(
              region,
              value,
              numeric.column = "value",
              col = color_i,
              pch = 16,
              cex = 0.7,
              ...
            )

            if (isTRUE(show_labels_i) && "name" %in% names(value)) {
              labels <- as.character(value[["name"]])
              keep <- !is.na(labels) & nzchar(labels)

              if (any(keep)) {
                circlize::circos.genomicText(
                  region[keep, , drop = FALSE],
                  value[keep, , drop = FALSE],
                  y = 0.75,
                  labels = labels[keep],
                  cex = 0.45,
                  facing = "inside",
                  niceFacing = TRUE,
                  col = color_i,
                  ...
                )
              }
            }
          }
        },
        track.height = height_i
      )
    })
  }

  if (!is.null(link_data) && nrow(link_data) > 0) {
    for (i in seq_len(nrow(link_data))) {
      circlize::circos.genomicLink(
        link_data[i, c("chr1", "start1", "end1"), drop = FALSE],
        link_data[i, c("chr2", "start2", "end2"), drop = FALSE],
        col = grDevices::adjustcolor("#7F77DD", alpha.f = 0.35)
      )
    }
  }

  if (isTRUE(show_legend)) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop(
        "Package \"ComplexHeatmap\" must be installed when show_legend = TRUE.",
        call. = FALSE
      )
    }

    legend_labels <- names(track_settings)

    if (length(legend_labels) > 0) {
      legend_cols <- vapply(
        track_settings,
        function(x) as.character(x$color)[1],
        character(1)
      )

      legend <- ComplexHeatmap::Legend(
        labels = legend_labels,
        legend_gp = grid::gpar(fill = legend_cols),
        title = "Tracks"
      )

      ComplexHeatmap::draw(legend)
    }
  }

  invisible(list(
    tracks = normalized_tracks,
    metadata = normalized_metadata,
    track_settings = track_settings,
    link_data = link_data,
    genome_build = genome_build,
    genome = genome_info,
    chromosomes = chromosomes,
    plot_type = resolved_plot_type,
    cytoband = cytoband,
    pdf_opened = pdf_opened
  ))
}
