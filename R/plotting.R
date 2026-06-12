#' ggplot2 theme for benchmark and feature-selection figures
#'
#' Minimal theme with black facet strips, inspired by the atlas feature
#' selection benchmark plotting utilities.
#'
#' @param ... Additional arguments passed to [ggplot2::theme_minimal()].
#'
#' @return A **ggplot2** `theme` object.
#'
#' @references
#' \url{https://github.com/theislab/atlas-feature-selection-benchmark/blob/b89fc0f66747062e6e1b4b35bd392b27ad035295/analysis/R/plotting.R}
#' @export
theme_features <- function(...) {
  ggplot2::theme_minimal(...) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(fill = NA),
      strip.text = ggplot2::element_text(colour = "white"),
      strip.background = ggplot2::element_rect(fill = "black")
    )
}

#' Publication-sized variant of [theme_features()]
#'
#' @param ... Additional arguments passed to [theme_features()].
#'
#' @return A **ggplot2** `theme` object.
#' @export
theme_features_pub <- function(...) {
  theme_features(base_size = 7, ...) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 7),
      plot.margin = ggplot2::margin(0.05, 0.05, 0.05, 0.05, "cm"),
      panel.spacing = grid::unit(0.06, "cm"),
      strip.text = ggplot2::element_text(size = 6),
      axis.title = ggplot2::element_text(size = 6),
      axis.text = ggplot2::element_text(size = 5),
      axis.ticks = ggplot2::element_line(linewidth = 0.25),
      legend.title = ggplot2::element_text(size = 6),
      legend.text = ggplot2::element_text(size = 5),
      legend.key.size = grid::unit(0.3, "cm"),
      legend.ticks = ggplot2::element_line(linewidth = 0.25)
    )
}
