#' cellfiers publication theme
#'
#' A clean theme in the style of the original workflow figures: classic axes,
#' a thin black panel border, black axis text, bold axis/strip titles and a
#' centred bold plot title. Use it on any ggplot, e.g.
#' `my_plot + theme_mif()`.
#'
#' @param base_size Base font size in points.
#'
#' @return A [ggplot2::theme()] object.
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) + geom_point() + theme_mif()
theme_mif <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        fill = NA, colour = "black", linewidth = 1
      ),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.line = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
}