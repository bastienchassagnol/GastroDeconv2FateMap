# ==========================================================================
# 0. Libraries and Filename settings ----
# ==========================================================================

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/single-cell"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

suppinger_single_cell <- readRDS(
  file = "./data/intermediate/suppinger_single_cell_2026-05-06.rds"
)
dim(suppinger_single_cell)

# ==========================================================================
# 1. Load single-cell object ----
# ==========================================================================
# Retrieve feature metadata from the SingleCellExperiment object
feature_metadata <- SummarizedExperiment::colData(suppinger_single_cell) |>
  as.data.frame() |>
  dplyr::select(
    "batch",
    "timepoints",
    "celltypeannotation",
    "celltype_colour"
  )

# ==========================================================================
# 2. Cell-type distribution by batch and time point ----
# ==========================================================================

# 2.1 Cell type level: Set levels for cell type and time point ----

time_point_levels <- levels(feature_metadata$timepoints)
cell_type_colours <- stats::setNames(
  as.character(unique(feature_metadata$celltype_colour)),
  unique(feature_metadata$celltypeannotation)
)

# 2.2 Calculate cellular ratios per time point and batch ----

single_cell_ratio_data <- feature_metadata |>
  dplyr::mutate(
    batch = factor(.data$batch),
    celltypeannotation = factor(
      .data$celltypeannotation,
      levels = names(cell_type_colours)
    )
  ) |>
  dplyr::group_by(.data$batch, .data$timepoints) |>
  dplyr::mutate(total_cells = dplyr::n()) |>
  dplyr::group_by(
    .data$batch,
    .data$timepoints,
    .data$celltypeannotation,
    .add = TRUE
  ) |>
  dplyr::summarise(
    cellular_ratio = dplyr::n() / dplyr::first(.data$total_cells),
    .groups = "drop"
  ) |>
  # required for ggplot2::geom_area() to work, as no null input is allowed
  tidyr::complete(
    batch,
    timepoints,
    celltypeannotation,
    fill = list(cellular_ratio = 0)
  ) |>
  dplyr::arrange(.data$batch, .data$timepoints, .data$celltypeannotation) |>
  dplyr::mutate(
    time_elapsed_hours = as.numeric(sub("h$", "", .data$timepoints))
  )

tinytable::tt(
  single_cell_ratio_data |>
    tail(),
  caption = "Single-cell cytometry distribution data"
)

# 2.3 Plot the stacked single-cell cytometry distribution ----
single_cell_distribution_plot <- ggplot2::ggplot(
  single_cell_ratio_data,
  ggplot2::aes(
    x = .data$time_elapsed_hours,
    y = .data$cellular_ratio,
    fill = .data$celltypeannotation,
    group = .data$celltypeannotation
  )
) +
  ggplot2::geom_area(
    alpha = 0.6,
    colour = "white",
    linewidth = 0.5,
    position = "stack"
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$batch),
    scales = "fixed",
    switch = "y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = unique(single_cell_ratio_data$time_elapsed_hours),
    labels = time_point_levels,
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.25),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::coord_cartesian(
    xlim = c(
      0,
      max(unique(single_cell_ratio_data$time_elapsed_hours), na.rm = TRUE)
    ),
    ylim = c(0, 1)
  ) +
  ggplot2::scale_fill_manual(
    values = cell_type_colours,
    na.value = "grey65"
  ) +
  ggplot2::theme_bw(base_size = 12, base_family = "") +
  ggplot2::labs(
    x = "Elapsed time (hours from experiment start)",
    y = "Cellular ratio",
    fill = "Cell type",
    title = glue::glue("{study} single-cell cytometry distribution"),
    subtitle = paste(
      "Stacked proportions over time;",
      "normalised within each batch and time point."
    ),
    caption = paste(
      "Areas stack to 100% at each time point.",
      "Facets correspond to two distinct cell lines."
    )
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 13, face = "bold"),
    legend.text = ggplot2::element_text(size = 8),
    legend.key.size = grid::unit(0.5, "cm"),
    legend.key.height = grid::unit(0.45, "cm"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing.y = grid::unit(3.25, "lines"),
    strip.placement = "outside",
    strip.background.y = ggplot2::element_rect(
      fill = "grey88",
      colour = "grey35",
      linewidth = 0.4
    ),
    strip.text.y.left = ggplot2::element_text(
      size = 14,
      face = "bold",
      angle = 0,
      colour = "grey10"
    ),
    plot.margin = ggplot2::margin(t = 6, r = 10, b = 6, l = 10),
    plot.title = ggplot2::element_text(size = 15, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 11),
    axis.text.y = ggplot2::element_text(size = 11),
    plot.caption = ggplot2::element_text(size = 9, hjust = 0)
  )

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    glue::glue("{study}_single_cell_cytometry_stacked_plot_{today}.pdf")
  ),
  plot = single_cell_distribution_plot,
  width = 14,
  height = 12,
  units = "in",
  dpi = 500,
  device = grDevices::cairo_pdf
)


# 2.4 Plot the alluvial cell composition changes ----

single_cell_alluvial_plot <- single_cell_ratio_data |>
  dplyr::mutate(
    timepoints = factor(.data$timepoints, levels = time_point_levels),
    celltypeannotation = factor(
      .data$celltypeannotation,
      levels = names(cell_type_colours)
    )
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = .data$timepoints,
      stratum = .data$celltypeannotation,
      alluvium = .data$celltypeannotation,
      y = .data$cellular_ratio,
      fill = .data$celltypeannotation
    )
  ) +
  ggalluvial::geom_alluvium(
    alpha = 0.55,
    curve_type = "sine",
    linewidth = 0.2
  ) +
  ggalluvial::geom_stratum(
    alpha = 0.88,
    colour = "white",
    linewidth = 0.25
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$batch),
    scales = "fixed",
    switch = "y"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  ggplot2::scale_fill_manual(
    values = cell_type_colours,
    na.value = "grey65"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::labs(
    x = "Time point",
    y = "Cellular ratio",
    fill = "Cell type",
    title = glue::glue("{study} single-cell cytometry Sankey plot"),
    subtitle = paste(
      "Cell composition changes over time;",
      "ratios are normalised within each batch and time point."
    ),
    caption = paste(
      "Each vertical axis sums to 100%.",
      "Facets correspond to two distinct cell lines."
    )
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 13, face = "bold"),
    legend.text = ggplot2::element_text(size = 8),
    legend.key.size = grid::unit(0.5, "cm"),
    legend.key.height = grid::unit(0.45, "cm"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing.y = grid::unit(3.25, "lines"),
    strip.placement = "outside",
    strip.background.y = ggplot2::element_rect(
      fill = "grey88",
      colour = "grey35",
      linewidth = 0.4
    ),
    strip.text.y.left = ggplot2::element_text(
      size = 14,
      face = "bold",
      angle = 0,
      colour = "grey10"
    ),
    plot.margin = ggplot2::margin(t = 6, r = 10, b = 6, l = 10),
    plot.title = ggplot2::element_text(size = 15, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 11),
    axis.text.y = ggplot2::element_text(size = 11),
    plot.caption = ggplot2::element_text(size = 9, hjust = 0)
  )

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    glue::glue("{study}_single_cell_cytometry_sankey_{today}.pdf")
  ),
  plot = single_cell_alluvial_plot,
  width = 14,
  height = 12,
  units = "in",
  dpi = 500,
  device = grDevices::cairo_pdf
)
