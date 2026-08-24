# Generate hexagonal stickers for GastroDeconv2Fate
#
# Two versions, matching the DeCovarT logo pair:
#   * epurated  — Waddington landscape, bifurcation, two fates
#   * detailed  — landscape plus deconvolution, drivers, forecast / KO
#
# Workflow follows hexSticker + magick as in
# https://nelson-gon.github.io/12/06/2020/hex-sticker-creation-r/
# Painted subplots live in inst/logo/; this script hex-crops them and
# also builds a fully code-driven ggplot2 Waddington (epurated_ggplot).
#
# Run from the package root, bypassing renv:
#   RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript --vanilla \
#     scripts/00_generate_hex_logo.R
#
# Project colours (see scripts/05_02_predict_morphotype_with_bulk_omics.R):
#   TLS / organised morphotype   #3182bd
#   neural-biased morphotype     #e6550d
# Sticker accents (DeCovarT-like navy + teal border):
#   fill #071422, border #2EC4B6

library(hexSticker)
library(magick)
library(ggplot2)

pkg_root <- if (requireNamespace("here", quietly = TRUE)) {
  here::here()
} else {
  normalizePath(".")
}

logo_dir <- file.path(pkg_root, "inst", "logo")
fig_dir <- file.path(pkg_root, "man", "figures")
if (!dir.exists(logo_dir)) {
  dir.create(logo_dir, recursive = TRUE)
}
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

col_fill <- "#071422"
col_border <- "#2EC4B6"
col_title <- "#F4FBFF"
col_url <- "#8FB8C9"

wrap_hex <- function(
  subplot_path,
  outfile,
  package = "GastroDeconv2Fate",
  p_size = 10.5,
  p_y = 0.32,
  s_y = 1.08,
  s_width = 1.62
) {
  sticker(
    subplot = subplot_path,
    package = package,
    p_size = p_size,
    p_y = p_y,
    p_color = col_title,
    p_family = "sans",
    s_x = 1,
    s_y = s_y,
    s_width = s_width,
    s_height = s_width,
    h_fill = col_fill,
    h_color = col_border,
    h_size = 1.35,
    url = "gastruloid fate  ·  TLS vs neural-biased",
    u_size = 3.4,
    u_color = col_url,
    u_y = 0.08,
    spotlight = FALSE,
    white_around_sticker = FALSE,
    filename = outfile,
    dpi = 400
  )
}

# --- ggplot2 Waddington (fully reproducible epurated subplot) ---------------
waddington_z <- function(time, fate) {
  split <- pmax(0, (time - 0.28) / 0.72)
  (1 - split) * 0.62 * fate^2 +
    split * 0.48 * (fate^2 - 1)^2 -
    0.07 * split * fate +
    0.22 * (1 - time)
}

ridge_lines <- lapply(seq(-1.55, 1.55, length.out = 36), function(fate) {
  time <- seq(0.02, 0.98, length.out = 90)
  z <- waddington_z(time, fate)
  data.frame(
    fate = fate,
    u = time + 0.46 * fate,
    v = z - 0.40 * time
  )
})
ridges <- do.call(rbind, ridge_lines)

marble <- data.frame(
  u = 0.34,
  v = waddington_z(0.34, 0) - 0.40 * 0.34 + 0.04
)
fates <- data.frame(
  u = c(1.00 - 0.46, 1.00 + 0.46),
  v = c(
    waddington_z(1, -1) - 0.40,
    waddington_z(1, 1) - 0.40
  ),
  fill = c("#2EC4B6", "#E85D04")
)

ggplot_waddington <- ggplot2::ggplot() +
  ggplot2::geom_path(
    data = ridges,
    ggplot2::aes(u, v, group = fate, colour = fate),
    linewidth = 0.35,
    alpha = 0.9
  ) +
  ggplot2::scale_colour_gradient2(
    low = "#2EC4B6",
    mid = "#F4FBFF",
    high = "#E85D04",
    midpoint = 0
  ) +
  ggplot2::geom_point(
    data = marble,
    ggplot2::aes(u, v),
    colour = "#F4FBFF",
    fill = "#7FDFF0",
    shape = 21,
    size = 5.5,
    stroke = 0.8
  ) +
  ggplot2::geom_point(
    data = fates,
    ggplot2::aes(u, v, fill = fill),
    colour = "#F4FBFF",
    shape = 21,
    size = 3.2,
    stroke = 0.4
  ) +
  ggplot2::scale_fill_identity() +
  ggplot2::coord_equal(xlim = c(-0.55, 1.55), ylim = c(-0.55, 0.72)) +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "none",
    plot.background = ggplot2::element_rect(fill = col_fill, colour = NA),
    panel.background = ggplot2::element_rect(fill = col_fill, colour = NA)
  )

ggplot_png <- file.path(logo_dir, "subplot_waddington_ggplot.png")
ggplot2::ggsave(
  ggplot_png,
  ggplot_waddington,
  width = 5,
  height = 5,
  dpi = 320,
  bg = col_fill
)

epurated_in <- file.path(logo_dir, "subplot_waddington_epurated.png")
detailed_in <- file.path(logo_dir, "subplot_waddington_detailed.png")
epurated_out <- file.path(logo_dir, "gastrodeconv2fate_hex_epurated.png")
detailed_out <- file.path(logo_dir, "gastrodeconv2fate_hex_detailed.png")
ggplot_out <- file.path(logo_dir, "gastrodeconv2fate_hex_epurated_ggplot.png")

wrap_hex(epurated_in, epurated_out, p_size = 11, p_y = 0.30, s_y = 1.10)
wrap_hex(detailed_in, detailed_out, p_size = 10, p_y = 0.28, s_y = 1.06)
wrap_hex(ggplot_png, ggplot_out, p_size = 11, p_y = 0.32, s_y = 1.05, s_width = 1.35)

# pkgdown / README uses the painted epurated hex as the package logo
file.copy(epurated_out, file.path(fig_dir, "logo.png"), overwrite = TRUE)

message(
  "Wrote:\n  ",
  epurated_out,
  "\n  ",
  detailed_out,
  "\n  ",
  ggplot_out
)
