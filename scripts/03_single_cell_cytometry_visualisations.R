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
    timepoints = factor(.data$timepoints, levels = time_point_levels),
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
    cell_count = dplyr::n(),
    cellular_ratio = dplyr::n() / dplyr::first(.data$total_cells),
    total_cells = dplyr::first(.data$total_cells),
    .groups = "drop"
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

saveRDS(
  single_cell_ratio_data,
  file = file.path(
    output_dir,
    glue::glue("{study}_single_cell_ratio_data_{today}.rds")
  )
)


# 2.3 Plot the stacked single-cell cytometry distribution ----
single_cell_distribution_plot <- ggplot2::ggplot(
  single_cell_ratio_data,
  ggplot2::aes(
    x = .data$timepoints,
    y = .data$cellular_ratio,
    fill = .data$celltypeannotation
  )
) +
  ggplot2::geom_col(
    alpha = 0.6,
    colour = "white",
    linewidth = 0.5,
    position = "stack"
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$batch),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) +
  ggplot2::scale_x_discrete(
    drop = TRUE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, by = 0.25),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
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
      "Stacked proportions over observed time points only;",
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
    glue::glue(
      "{study}_single_cell_celltype_abundances_stacked_plot_{today}.pdf"
    )
  ),
  plot = single_cell_distribution_plot,
  width = 14,
  height = 12,
  units = "in",
  dpi = 500,
  device = grDevices::cairo_pdf
)


# 2.4 Plot absolute cell abundances over time ----

single_cell_total_data <- single_cell_ratio_data |>
  dplyr::distinct(
    .data$batch,
    .data$timepoints,
    .data$time_elapsed_hours,
    .data$total_cells
  )

time_axis_labels <- single_cell_ratio_data |>
  dplyr::distinct(.data$time_elapsed_hours, .data$timepoints) |>
  dplyr::arrange(.data$time_elapsed_hours)

single_cell_abundance_line_plot <- ggplot2::ggplot() +
  ggplot2::geom_line(
    data = single_cell_ratio_data,
    ggplot2::aes(
      x = .data$time_elapsed_hours,
      y = .data$cell_count,
      colour = .data$celltypeannotation,
      group = .data$celltypeannotation
    ),
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    data = single_cell_ratio_data,
    ggplot2::aes(
      x = .data$time_elapsed_hours,
      y = .data$cell_count,
      colour = .data$celltypeannotation
    ),
    size = 2
  ) +
  ggplot2::geom_line(
    data = single_cell_total_data,
    ggplot2::aes(
      x = .data$time_elapsed_hours,
      y = .data$total_cells
    ),
    colour = "black",
    linewidth = 1.6,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = single_cell_total_data,
    ggplot2::aes(
      x = .data$time_elapsed_hours,
      y = .data$total_cells
    ),
    colour = "black",
    size = 2.5,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(.data$batch),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) +
  ggplot2::scale_x_continuous(
    breaks = time_axis_labels$time_elapsed_hours,
    labels = time_axis_labels$timepoints
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::comma_format(),
    expand = ggplot2::expansion(mult = c(0.02, 0.06))
  ) +
  ggplot2::scale_colour_manual(
    values = cell_type_colours,
    na.value = "grey65"
  ) +
  ggplot2::theme_bw(base_size = 12, base_family = "") +
  ggplot2::labs(
    x = "Elapsed time (hours from experiment start)",
    y = "Cell count",
    colour = "Cell type",
    title = glue::glue("{study} single-cell cytometry abundances"),
    subtitle = paste(
      "Absolute cell counts per type over observed time points;",
      "thick black line shows total cells per batch and time point."
    ),
    caption = paste(
      "Lines connect observed time points only.",
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
    glue::glue(
      "{study}_single_cell_celltype_abundances_lines_{today}.pdf"
    )
  ),
  plot = single_cell_abundance_line_plot,
  width = 14,
  height = 12,
  units = "in",
  dpi = 500,
  device = grDevices::cairo_pdf
)


# 2.5 Plot the alluvial cell composition changes ----

single_cell_alluvial_plot <- single_cell_ratio_data |>
  dplyr::mutate(
    timepoints = factor(.data$timepoints, levels = time_point_levels),
    celltypeannotation = factor(
      .data$celltypeannotation,
      levels = names(cell_type_colours)
    )
  ) |>
  dplyr::group_by(.data$batch, .data$timepoints) |>
  dplyr::filter(sum(.data$cellular_ratio, na.rm = TRUE) > 0) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    timepoints = factor(
      .data$timepoints,
      levels = time_point_levels[time_point_levels %in% .data$timepoints]
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
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) +
  ggplot2::scale_x_discrete(drop = TRUE) +
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
    glue::glue("{study}_single_cell_celltype_abundances_sankey_{today}.pdf")
  ),
  plot = single_cell_alluvial_plot,
  width = 14,
  height = 12,
  units = "in",
  dpi = 500,
  device = grDevices::cairo_pdf
)


# 2.6 Plot the mosaic (Marimekko) timeline ----

single_cell_mosaic_data <- single_cell_ratio_data |>
  dplyr::filter(.data$cell_count > 0) |>
  dplyr::mutate(
    timepoints = factor(.data$timepoints, levels = time_point_levels),
    celltypeannotation = factor(
      .data$celltypeannotation,
      levels = names(cell_type_colours)
    ),
    mosaic_label = paste0(
      scales::comma(.data$cell_count),
      "\n",
      scales::percent(.data$cellular_ratio, accuracy = 0.1)
    )
  )

mosaic_batch_theme <- ggplot2::theme(
  legend.position = "bottom",
  legend.title = ggplot2::element_text(size = 13, face = "bold"),
  legend.text = ggplot2::element_text(size = 8),
  legend.key.size = grid::unit(0.5, "cm"),
  legend.key.height = grid::unit(0.45, "cm"),
  panel.grid.major.x = ggplot2::element_blank(),
  panel.grid.minor = ggplot2::element_blank(),
  plot.margin = ggplot2::margin(t = 6, r = 10, b = 6, l = 10),
  plot.title = ggplot2::element_text(size = 14, face = "bold"),
  axis.text.y = ggplot2::element_text(size = 11)
)

single_cell_mosaic_plot_list <- single_cell_mosaic_data |>
  dplyr::group_split(.data$batch) |>
  purrr::map(function(batch_data) {
    batch_label <- as.character(batch_data$batch[1])

    batch_timepoints <- batch_data |>
      dplyr::distinct(.data$timepoints, .data$time_elapsed_hours) |>
      dplyr::arrange(.data$time_elapsed_hours) |>
      dplyr::pull(.data$timepoints) |>
      as.character()

    batch_data <- batch_data |>
      dplyr::filter(.data$timepoints %in% batch_timepoints) |>
      dplyr::mutate(
        timepoints = factor(.data$timepoints, levels = batch_timepoints)
      )

    mosaic_base <- ggplot2::ggplot(batch_data) +
      ggmosaic::geom_mosaic(
        ggplot2::aes(
          x = ggmosaic::product(timepoints),
          fill = celltypeannotation,
          weight = cell_count
        ),
        alpha = 0.6,
        colour = "white",
        linewidth = 0.5
      ) +
      ggmosaic::scale_x_productlist() +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        breaks = seq(0, 1, by = 0.25),
        expand = ggplot2::expansion(mult = c(0, 0.02))
      )

    mosaic_label_data <- ggplot2::layer_data(mosaic_base, 1) |>
      dplyr::filter(.data$.wt > 0) |>
      dplyr::group_by(.data$x__timepoints) |>
      dplyr::mutate(
        cellular_ratio = .data$.wt / sum(.data$.wt),
        mosaic_label = paste0(
          scales::comma(.data$.wt),
          "\n",
          scales::percent(.data$cellular_ratio, accuracy = 0.1)
        ),
        label_x = (.data$xmin + .data$xmax) / 2,
        label_y = (.data$ymin + .data$ymax) / 2
      ) |>
      dplyr::ungroup()

    mosaic_base +
      ggplot2::geom_text(
        data = mosaic_label_data,
        ggplot2::aes(
          x = .data$label_x,
          y = .data$label_y,
          label = .data$mosaic_label
        ),
        inherit.aes = FALSE,
        size = 2.2,
        colour = "grey10",
        lineheight = 0.9
      ) +
      ggplot2::scale_fill_manual(
        values = cell_type_colours,
        na.value = "grey65",
        drop = TRUE
      ) +
      ggmosaic::theme_mosaic() +
      ggplot2::theme_bw(base_size = 12, base_family = "") +
      ggplot2::labs(
        title = batch_label,
        x = "Elapsed time (hours from experiment start)",
        y = "Cellular ratio",
        fill = "Cell type"
      ) +
      mosaic_batch_theme
  })

single_cell_mosaic_plot <- patchwork::wrap_plots(
  single_cell_mosaic_plot_list,
  ncol = 1,
  guides = "collect"
) +
  patchwork::plot_annotation(
    title = glue::glue("{study} single-cell cytometry mosaic timeline"),
    subtitle = paste(
      "Bar width reflects absolute cell counts;",
      "stacked height reflects relative proportions within each time point."
    ),
    caption = paste(
      "Marimekko layout: wider bars indicate more cells sampled.",
      "Tile labels show absolute count and within-time-point ratio.",
      "Panels correspond to two distinct cell lines."
    )
  ) &
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 15, face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 11),
    plot.caption = ggplot2::element_text(size = 9, hjust = 0),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    glue::glue("{study}_single_cell_celltype_abundances_mosaic_{today}.pdf")
  ),
  plot = single_cell_mosaic_plot,
  width = 14,
  height = 16,
  units = "in",
  dpi = 500,
  device = grDevices::cairo_pdf
)
