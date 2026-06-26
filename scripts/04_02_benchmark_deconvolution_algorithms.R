# ==========================================================================
# 0. Libraries and Filename settings ----
# ==========================================================================

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/deconvolution_benchmark"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

source("./R/deconvolution_scorings.R")
source("./R/plotting.R")

library(ggplot2)

# ==========================================================================
# 1. Load omic objects----
# ==========================================================================

# 1.1 Load single-cell object ----

single_cell_ratio_data <- readRDS(file.path(
  "outputs/single-cell",
  "suppinger_single_cell_ratio_data_2026-06-24.rds"
))
tinytable::tt(
  single_cell_ratio_data,
  caption = "Single-cell cytometry distribution data"
)

# 1.2 Load deconvolution results ----
deconv_results <- readRDS(file.path(
  "outputs/deconvolution",
  "suppinger_deconvolution_results_2026-05-07.rds"
))
tinytable::tt(
  deconv_results,
  caption = "Deconvolution estimates"
)

# 1.3 join single-cell and deconvolution results

shared_timepoints <- intersect(
  unique(single_cell_ratio_data$timepoints),
  unique(deconv_results$time_point)
)
shared_cell_types <- intersect(
  unique(single_cell_ratio_data$celltypeannotation),
  unique(deconv_results$cell_type)
) |>
  sort()

single_cell_ratios_wider <- single_cell_ratio_data |>
  dplyr::filter(
    timepoints %in%
      shared_timepoints &
      celltypeannotation %in% shared_cell_types
  ) |>
  tidyr::pivot_wider(
    id_cols = c("batch", "timepoints"),
    names_from = "celltypeannotation",
    values_from = "cellular_ratio",
    values_fill = 0
  ) |>
  dplyr::group_by(timepoints) |>
  dplyr::summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::select(timepoints, all_of(shared_cell_types)) |>
  dplyr::mutate(
    timepoints = factor(
      as.character(timepoints),
      levels = c("48h", "72h", "96h"),
      ordered = TRUE
    )
  ) |>
  dplyr::arrange(timepoints)

tinytable::tt(
  single_cell_ratios_wider  |> 
  dplyr::filter(timepoints == "96h"),
  caption = "Single-cell cytometry distribution data at 96h"
)


deconv_results_wider <- deconv_results |>
  dplyr::filter(
    time_point %in%
      shared_timepoints &
      deconvolution_algorithm != "Scaden" &
      cell_type %in% shared_cell_types
  ) |>
  tidyr::pivot_wider(
    id_cols = c("deconvolution_algorithm", "time_point", "sample_id"),
    names_from = "cell_type",
    values_from = "cell_proportion",
    values_fill = 0
  ) |>
  dplyr::select(
    deconvolution_algorithm,
    time_point,
    sample_id,
    all_of(shared_cell_types)
  ) |>
  dplyr::mutate(
    sample_id = stringr::str_extract(sample_id, "rep\\d+$"),
    time_point = factor(
      as.character(time_point),
      levels = c("48h", "72h", "96h"),
      ordered = TRUE
    )
  ) |>
  dplyr::arrange(time_point, deconvolution_algorithm, sample_id)


  tinytable::tt(
    deconv_results_wider  |> 
    dplyr::filter(time_point == "96h" & deconvolution_algorithm == "BayesPrism"),
    caption = "Deconvolution estimates for BayesPrism at 96h"
  )


# ==========================================================================
# 2. Score each deconvolution replicate (purrr) ----
# ==========================================================================

# 2.1: Provide a shared vector of cell types ----

.composition_vector <- function(row, cell_types) {
  vals <- unlist(row[cell_types], use.names = FALSE)
  stats::setNames(vals, cell_types)
}

benchmark_metrics <- purrr::map_dfr(
  seq_len(nrow(deconv_results_wider)),
  \(i) {
    row <- deconv_results_wider[i, ]
    obs_row <- single_cell_ratios_wider[
      single_cell_ratios_wider$timepoints == row$time_point,
      ,
      drop = FALSE
    ]
    p_obs <- .composition_vector(obs_row, shared_cell_types)
    p_estimated <- .composition_vector(row, shared_cell_types)

    dplyr::bind_cols(
      row |>
        dplyr::select(
          deconvolution_algorithm,
          time_point,
          sample_id
        ),
      tibble::tibble(
        pearson = eval_Pearson(p_obs, p_estimated),
        jsd = eval_JSD(p_obs, p_estimated),
        aitchison = eval_Aitchison(p_obs, p_estimated),
        sdid = eval_SDID(p_obs, p_estimated),
        rmse = eval_RMSE(p_obs, p_estimated),
        mae = eval_MAE(p_obs, p_estimated)
      )
    )
  }
)

benchmark_hrrmse <- deconv_results_wider |>
  dplyr::group_by(deconvolution_algorithm, time_point) |>
  dplyr::group_modify(\(df, group_keys) {
    obs_row <- single_cell_ratios_wider[
      single_cell_ratios_wider$timepoints == group_keys$time_point,
      ,
      drop = FALSE
    ]
    p_obs <- .composition_vector(obs_row, shared_cell_types)
    p_estimated <- as.matrix(df[, shared_cell_types, drop = FALSE])
    rownames(p_estimated) <- df$sample_id
    tibble::tibble(hrrmse = eval_hrRMSE(p_obs, p_estimated))
  }) |>
  dplyr::ungroup()

benchmark_metrics <- benchmark_metrics |>
  dplyr::left_join(
    benchmark_hrrmse,
    by = c("deconvolution_algorithm", "time_point")
  )

saveRDS(
  benchmark_metrics,
  file.path(output_dir, glue::glue("{study}_benchmark_metrics_{today}.rds"))
)


# ==========================================================================
# 3. Scatter plot: observed vs estimated (rows = sample, cols = time point) ----
# ==========================================================================

joint_composition <- deconv_results_wider |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(shared_cell_types),
    names_to = "cell_type",
    values_to = "estimated"
  ) |>
  dplyr::rename(timepoints = time_point) |>
  dplyr::mutate(
    timepoints = factor(
      as.character(timepoints),
      levels = c("48h", "72h", "96h"),
      ordered = TRUE
    ),
    sample_id = factor(
      sample_id,
      levels = sort(unique(sample_id)),
      ordered = TRUE
    )
  ) |>
  dplyr::left_join(
    single_cell_ratios_wider |>
      dplyr::mutate(
        timepoints = factor(
          as.character(timepoints),
          levels = c("48h", "72h", "96h"),
          ordered = TRUE
        )
      ) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(shared_cell_types),
        names_to = "cell_type",
        values_to = "observed"
      ),
    by = c("timepoints", "cell_type")
  )  |> 
  dplyr::filter(!(observed == 0 & estimated == 0))

algorithm_levels <- sort(unique(joint_composition$deconvolution_algorithm))
algorithm_shapes <- stats::setNames(
  c(16L, 17L, 15L),
  algorithm_levels
)

scatter_pearson_labels <- benchmark_metrics |>
  dplyr::mutate(
    timepoints = factor(
      as.character(time_point),
      levels = c("48h", "72h", "96h"),
      ordered = TRUE
    ),
    sample_id = factor(
      sample_id,
      levels = levels(joint_composition$sample_id),
      ordered = TRUE
    )
  ) |>
  dplyr::arrange(deconvolution_algorithm) |>
  dplyr::group_by(sample_id, timepoints) |>
  dplyr::summarise(
    label = paste(
      glue::glue(
        "{deconvolution_algorithm}: r = ",
        "{formatC(pearson, format = 'f', digits = 3)}"
      ),
      collapse = "\n"
    ),
    .groups = "drop"
  )

plot_scatter_panel <- function(show_legend = FALSE) {
  p <- ggplot2::ggplot(
    joint_composition,
    ggplot2::aes(
      x = observed,
      y = estimated,
      colour = cell_type,
      shape = deconvolution_algorithm,
    )
  ) +
    ggplot2::geom_point(alpha = 0.85, size = 2.2) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = scatter_pearson_labels,
      ggplot2::aes(
        x = 0.05,
        y = 0.95,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.6,
      lineheight = 0.9
    ) +
    ggplot2::scale_shape_manual(values = algorithm_shapes) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.08, 0.05))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.08, 0.02))
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(sample_id),
      cols = ggplot2::vars(timepoints),
      labeller = ggplot2::label_value
    ) +
    ggplot2::labs(
      x = "Observed proportion",
      y = "Estimated proportion",
      colour = "Cell type",
      shape = "Algorithm"
    ) +
    theme_features(base_size = 11) +
    ggplot2::theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.box = "horizontal",
      axis.title = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(0.8, "lines")
    )

  if (show_legend) {
    p <- p +
      ggplot2::guides(
        colour = ggplot2::guide_legend(nrow = 2, title.position = "top"),
        shape = ggplot2::guide_legend(nrow = 1, title.position = "top")
      )
  }

  p
}

scatter_panels <- plot_scatter_panel(show_legend = FALSE)
scatter_legend <- cowplot::get_legend(plot_scatter_panel(show_legend = TRUE))

scatter_figure <- cowplot::plot_grid(
  scatter_panels,
  scatter_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

ggplot2::ggsave(
  file.path(
    output_dir,
    glue::glue("{study}_observed_vs_estimated_scatter_{today}.pdf")
  ),
  scatter_figure,
  width = 16,
  height = 14,
  dpi = 300,
  limitsize = FALSE
)

# ==========================================================================
# 4. Radar plots (mean metrics across replicates) ----
# ==========================================================================

n_cell_types <- length(shared_cell_types)
mae_upper <- 2 / n_cell_types
rmse_upper <- sqrt(2 / n_cell_types)

radar_long <- benchmark_metrics |>
  dplyr::group_by(time_point, deconvolution_algorithm) |>
  dplyr::summarise(
    dplyr::across(
      c(pearson, jsd, aitchison, sdid, rmse, mae),
      ~ mean(.x, na.rm = TRUE)
    ),
    hrrmse = dplyr::first(hrrmse),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Aitchison = 1 - aitchison / max(aitchison, na.rm = TRUE),
    SDID = 1 - sdid,
    Pearson = pmax(0, pearson),
    JSD = 1 - jsd,
    RMSE = pmax(0, 1 - rmse / rmse_upper),
    MAE = pmax(0, 1 - mae / mae_upper),
    hrRMSE = pmax(0, 1 - hrrmse)
  ) |>
  dplyr::select(
    time_point,
    deconvolution_algorithm,
    Pearson,
    JSD,
    Aitchison,
    SDID,
    RMSE,
    MAE,
    hrRMSE
  ) |>
  tidyr::pivot_longer(
    cols = c(Pearson, JSD, Aitchison, SDID, RMSE, MAE, hrRMSE),
    names_to = "metric",
    values_to = "score"
  ) |>
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c(
        "Pearson",
        "JSD",
        "Aitchison",
        "SDID",
        "RMSE",
        "MAE",
        "hrRMSE"
      )
    )
  )

# generate a function to plot a radar plot for each time point
algorithm_levels <- sort(unique(radar_long$deconvolution_algorithm))
algorithm_colours <- stats::setNames(
  scales::hue_pal()(length(algorithm_levels)),
  algorithm_levels
)

plot_radar_one_timepoint <- function(tp) {
  tp_data <- radar_long |>
    dplyr::filter(time_point == tp)

  tp_data_wide <- tp_data |>
    dplyr::select(deconvolution_algorithm, metric, score) |>
    tidyr::pivot_wider(
      names_from = metric,
      values_from = score
    ) |>
    dplyr::rename(group = deconvolution_algorithm)

  ggradar::ggradar(
    tp_data_wide,
    values.radar = c("0", "0.5", "1"),
    grid.min = 0,
    grid.max = 1,
    group.line.width = 0.7,
    group.point.size = 2.2,
    axis.label.size = 3.5,
    grid.label.size = 3,
    legend.text.size = 8,
    group.colours = unname(algorithm_colours[tp_data_wide$group])
  ) +
    ggplot2::labs(
      title = tp,
      x = NULL,
      y = NULL,
      colour = "Algorithm",
      fill = "Algorithm"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major = ggplot2::element_line(colour = "grey85"),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}
radar_plots <- purrr::map(
  unique(radar_long$time_point),
  plot_radar_one_timepoint
)
names(radar_plots) <- unique(radar_long$time_point)

radar_legend_plot <- plot_radar_one_timepoint(unique(radar_long$time_point)[[
  1L
]]) +
  ggplot2::theme(legend.position = "bottom")
radar_legend <- cowplot::get_legend(radar_legend_plot)
radar_panels <- patchwork::wrap_plots(radar_plots, nrow = 1)

radar_figure <- cowplot::plot_grid(
  radar_panels,
  radar_legend,
  ncol = 1,
  rel_heights = c(1, 0.12)
)

ggplot2::ggsave(
  file.path(
    output_dir,
    glue::glue("{study}_metric_radar_{today}.pdf")
  ),
  radar_figure,
  width = 14,
  height = 5.5
)

# ==========================================================================
# 5. Metric correlation heatmap (ellipse correlogram) ----
# ==========================================================================

metric_cor <- benchmark_metrics |>
  dplyr::select(pearson, jsd, aitchison, sdid, rmse, mae, hrrmse) |>
  stats::cor(use = "pairwise.complete.obs")

palette_cor <- grDevices::colorRampPalette(
  RColorBrewer::brewer.pal(5, "Spectral")
)(100)
cor_to_palette <- function(r) {
  palette_cor[pmax(1L, pmin(100L, round(r * 50 + 50)))]
}

grDevices::pdf(
  file.path(
    output_dir,
    glue::glue("{study}_metric_correlogram_{today}.pdf")
  ),
  width = 7.5,
  height = 7
)

graphics::layout(
  matrix(c(1L, 2L), nrow = 1L),
  widths = c(0.12, 0.88)
)

# Correlation colour scale (left; narrow panel needs tight margins)
par(mar = c(2, 0, 2, 0.8))
cor_scale_y <- seq(-1, 1, length.out = length(palette_cor))
graphics::image(
  x = 1,
  y = cor_scale_y,
  z = matrix(cor_scale_y, nrow = 1),
  col = palette_cor,
  axes = FALSE,
  xlab = "",
  ylab = ""
)
graphics::axis(
  4,
  at = seq(-1, 1, by = 0.5),
  labels = format(seq(-1, 1, by = 0.5), nsmall = 1),
  las = 1,
  cex.axis = 0.65,
  mgp = c(0.5, 0.3, 0)
)
graphics::mtext("r", side = 4, line = 1.5, cex = 0.8)

# Correlogram (right)
par(mar = c(1, 1, 1, 1))
ellipse::plotcorr(
  metric_cor,
  col = cor_to_palette(metric_cor),
  mar = c(1, 1, 1, 1)
)
grDevices::dev.off()
