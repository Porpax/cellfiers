# encoding: UTF-8
# Internal helpers -------------------------------------------------------

#' Find a column in a data frame
#'
#' Look for a column by exact name first, then by a byte-level normalised
#' name. Normalisation lowercases, keeps only ASCII letters/digits, and maps
#' the micro sign used by QuPath coordinates (U+00B5, "Centroid X µm") onto
#' "u", so that names read with readr (spaces preserved) or read.csv (spaces
#' as dots) are matched regardless of R string encoding marks.
#'
#' @param data A data frame.
#' @param candidates Character vector of candidate column names.
#' @return The name of the first matching column, or NULL if none matches.
#' @keywords internal
.find_col <- function(data, candidates) {
  if (is.null(candidates)) return(NULL)
  for (c in candidates) {
    if (!is.na(c) && nzchar(c) && c %in% names(data)) return(c)
  }
  keys <- vapply(names(data), .norm_key_str, character(1))
  for (c in candidates) {
    if (is.na(c) || !nzchar(c)) next
    ck <- .norm_key_str(c)
    hit <- which(keys == ck)
    if (length(hit)) return(names(data)[hit[1]])
  }
  NULL
}

#' Byte-level normalisation key for a string
#'
#' Converts a string to a normalised ASCII key (lowercase letters and
#' digits), mapping the micro sign (U+00B5) to "u" so that "Centroid X µm"
#' matches "Centroid X um". Works on raw bytes so it is independent of R
#' string encoding marks.
#'
#' @param x A character string.
#' @return An integer vector of ASCII codes.
#' @keywords internal
.norm_key <- function(x) {
  ints <- as.integer(charToRaw(x))
  out <- integer(0)
  i <- 1L
  n <- length(ints)
  while (i <= n) {
    b <- ints[i]
    if (i < n && b == 0xC2 && ints[i + 1L] == 0xB5) {
      out <- c(out, 0x75L)
      i <- i + 2L
    } else if (b == 0xB5) {
      out <- c(out, 0x75L)
      i <- i + 1L
    } else {
      out <- c(out, b)
      i <- i + 1L
    }
  }
  lc <- out >= 0x41L & out <= 0x5AL
  out[lc] <- out[lc] + 32L
  keep <- (out >= 0x61L & out <= 0x7AL) | (out >= 0x30L & out <= 0x39L)
  out[keep]
}

#' @rdname dot-norm_key
#' @keywords internal
.norm_key_str <- function(x) paste(.norm_key(x), collapse = ",")

#' Normalise marker status values
#'
#' Convert marker status columns holding "+/-", "1/0", TRUE/FALSE, or any
#' common spelling of positive/negative into a canonical "pos"/"neg"/NA
#' character vector.
#'
#' @param x A vector of marker status values.
#' @return A character vector with values "pos", "neg" or NA.
#' @keywords internal
.normalize_marker <- function(x) {
  s <- as.character(x)
  out <- rep(NA_character_, length(s))
  out[!is.na(s) & s %in% c("+", "1", "true", "t", "positive", "pos")] <- "pos"
  out[!is.na(s) & s %in% c("-", "0", "false", "f", "negative", "neg")] <- "neg"
  out
}

#' Add a cell id column
#'
#' Adds a unique cell_id column of the form "<sample>_<index>", where
#' <index> is the row number within each sample.
#'
#' @param data A data frame.
#' @param sample_col Name of the sample column.
#' @param id_col Name of the new id column.
#' @return The data frame with an extra id column.
#' @keywords internal
.add_cell_id <- function(data, sample_col, id_col) {
  if (id_col %in% names(data)) return(data)
  if (!sample_col %in% names(data)) {
    data[[id_col]] <- paste0("cell_", seq_len(nrow(data)))
    return(data)
  }
  data[[id_col]] <- paste0(data[[sample_col]], "_", stats::ave(seq_len(nrow(data)), data[[sample_col]], FUN = seq_along))
  data
}

#' Extract a sample id from a file name
#'
#' @param file A file path.
#' @param pattern A regular expression with exactly one capture group. The
#'   captured group is used as the sample id.
#' @return The sample id (character).
#' @keywords internal
.extract_sample_id <- function(file, pattern) {
  m <- stringr::str_match(basename(file), pattern)
  if (is.na(m[1, 1]) || ncol(m) < 2 || is.na(m[1, 2])) {
    stop("Could not extract a sample id from file: ", basename(file),
         "\n  Provide sample_id explicitly or a sample_pattern with one capture group.",
         call. = FALSE)
  }
  m[1, 2]
}

#' Round a value up to a "nice" step
#'
#' Returns the smallest "nice" number (1, 2, 2.5, 5 or 10 times a power of
#' ten) that is >= x. Used to compute radar axis labels.
#'
#' @param x A positive number.
#' @return A nice step value.
#' @keywords internal
.nice_step <- function(x) {
  if (length(x) == 0L || is.na(x) || x <= 0) return(1)
  exp10 <- floor(log10(x))
  base <- 10^exp10
  for (m in c(1, 2, 2.5, 5, 10)) {
    if (m * base >= x) return(m * base)
  }
  10 * base
}