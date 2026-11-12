
#' Default categorical colour palette
#'
#' Returns a colour-blind friendly palette of the requested length, reusing
#' and shuffling colours when more categories than palette entries are needed.
#'
#' @param n Number of colours requested.
#' @return A character vector of colours.
#' @keywords internal
.default_palette <- function(n) {
  base <- c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
  )
  if (n <= length(base)) {
    out <- base[seq_len(n)]
  } else {
    out <- rep(base, ceiling(n / length(base)))[seq_len(n)]
  }
  set.seed(2026)
  sample(out)
}
#' Turbo-style colour ramp
#'
#' Returns a turbo-like colour ramp (30 colours) using `hcl.colors()` when
#' available, with a close approximation for older R versions that lack the
#' "turbo" palette.
#'
#' @return A character vector of 30 colours.
#' @keywords internal
.turbo_palette <- function() {
  if (requireNamespace("viridis", quietly = TRUE)) {
    viridis::viridis(30, option = "turbo")
  } else if ("turbo" %in% grDevices::hcl.pals()) {
    grDevices::hcl.colors(30, palette = "turbo")
  } else {
    grDevices::colorRampPalette(c(
      "#30123B", "#3B1B71", "#4651A0", "#2E86C1", "#57B8A7",
      "#A8DA5B", "#F9E84C", "#F7A81B", "#E6370C", "#7A0403"
    ))(30)
  }
}