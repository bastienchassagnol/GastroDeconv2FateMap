# ==========================================================================
# 0. Libraries and Filename settings ----
# ==========================================================================
study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/deconvolution"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ==========================================================================
# 1. Load the cytometry data ----
# ==========================================================================

deconv_results <- readRDS(file.path(
  output_dir,
  glue::glue("{study}_deconvolution_results_{today}.rds")
))

dim(deconv_results)
tinytable::tt(
  deconv_results,
  caption = "Deconvolution estimates"
)

# Plot the deconvolution estimates
algorithm_colours <- c(
  "BayesPrism" = "#e76f51",
  "DWLS" = "#2a9d8f",
  "SCDC" = "#457b9d",
  "Scaden" = "#6a994e"
)

# ==========================================================================
# 2. Deconvolution violin clouds ----
# ==========================================================================
time_points <- unique(deconv_results$time_point)
time_points <- time_points[order(as.numeric(sub("h$", "", time_points)))]

plot_list <- purrr::map(
  time_points,
  \(tp) {
    plot_data <- deconv_results |>
      dplyr::filter(time_point == tp) |>
      dplyr::mutate(
        cell_type = factor(cell_type, levels = unique(cell_type))
      )

    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = cell_type,
        y = cell_proportion,
        color = deconvolution_algorithm,
        fill = deconvolution_algorithm,
        group = interaction(cell_type, deconvolution_algorithm)
      )
    ) +
      ggplot2::geom_vline(
        xintercept = seq(
          1.5,
          nlevels(plot_data$cell_type) - 0.5,
          by = 1
        ),
        colour = "grey92",
        linetype = "dashed",
        linewidth = 0.25,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_boxplot(
        position = ggplot2::position_dodge(width = 0.55),
        width = 0.5,
        alpha = 0.45,
        outlier.shape = NA,
        linewidth = 0.35
      ) +
      ggplot2::geom_point(
        position = ggplot2::position_jitterdodge(
          jitter.width = 0.12,
          dodge.width = 0.55
        ),
        size = 3.6,
        alpha = 0.8
      ) +
      ggplot2::stat_summary(
        ggplot2::aes(group = interaction(cell_type, deconvolution_algorithm)),
        fun = mean,
        geom = "point",
        shape = 23,
        size = 3,
        fill = "white",
        color = "black",
        position = ggplot2::position_dodge(width = 0.55),
        show.legend = FALSE
      ) +
      ggplot2::geom_hline(
        yintercept = 1,
        linewidth = 0.8,
        linetype = "dashed",
        color = "grey20"
      ) +
      ggplot2::scale_color_manual(
        values = algorithm_colours,
        breaks = intersect(
          names(algorithm_colours),
          unique(plot_data$deconvolution_algorithm)
        )
      ) +
      ggplot2::scale_fill_manual(
        values = algorithm_colours,
        breaks = intersect(
          names(algorithm_colours),
          unique(plot_data$deconvolution_algorithm)
        )
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.25)
      ) +
      ggplot2::coord_flip() +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::labs(
        x = "Cell type",
        y = "Cell proportion",
        color = "Deconvolution algorithm",
        fill = "Deconvolution algorithm",
        title = glue::glue("Deconvolution estimates - {tp}"),
        caption = "White diamonds show means."
      ) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 13, face = "bold"),
        legend.text = ggplot2::element_text(size = 11),
        legend.key.size = grid::unit(0.8, "cm"),
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 12),
        plot.caption = ggplot2::element_text(size = 9, hjust = 0)
      )
  }
)
deconv_estimates_plot <- gridExtra::marrangeGrob(
  grobs = plot_list,
  nrow = 1,
  ncol = 1,
  top = NULL
)
# Save the plot
ggplot2::ggsave(
  file.path(
    output_dir,
    glue::glue("{study}_deconvolution_estimates_{today}.pdf")
  ),
  plot = deconv_estimates_plot,
  width = 18,
  height = 9
)


# ============================================================================
# 3. Alluvial plots (one page per algorithm, single PDF) ----
# ============================================================================

deconv_algorithms <- unique(deconv_results$deconvolution_algorithm)
deconv_algorithms <- c(
  intersect(names(algorithm_colours), deconv_algorithms),
  sort(setdiff(deconv_algorithms, names(algorithm_colours)))
)

alluvial_plot_list <- purrr::map(
  deconv_algorithms,
  \(algo) {
    algo_agg <- deconv_results |>
      dplyr::filter(.data$deconvolution_algorithm == algo) |>
      dplyr::group_by(.data$time_point, .data$cell_type) |>
      dplyr::summarise(
        cell_proportion = mean(.data$cell_proportion, na.rm = TRUE),
        .groups = "drop"
      )

    cell_levels <- sort(unique(as.character(algo_agg$cell_type)))
    plot_agg <- dplyr::mutate(
      algo_agg,
      time_point = factor(.data$time_point, levels = time_points),
      cell_type = factor(.data$cell_type, levels = cell_levels)
    )

    ggplot2::ggplot(
      plot_agg,
      ggplot2::aes(
        x = .data$time_point,
        stratum = .data$cell_type,
        alluvium = .data$cell_type,
        y = .data$cell_proportion,
        fill = .data$cell_type
      )
    ) +
      ggalluvial::geom_alluvium(
        alpha = 0.55,
        curve_type = "sine",
        linewidth = 0.2
      ) +
      ggalluvial::geom_stratum(alpha = 0.88, linewidth = 0.3) +
      ggplot2::scale_y_continuous(
        labels = scales::percent,
        limits = c(0, 1),
        expand = c(0, 0)
      ) +
      ggplot2::scale_fill_viridis_d(option = "C", end = 0.95) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        panel.grid.major.x = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        title = glue::glue("{algo} deconvolution - {study}"),
        subtitle = "Mean cell proportion per time point.",
        x = "Time point",
        y = "Cell proportion",
        fill = "Cell type"
      )
  }
)

deconv_alluvial_plot <- gridExtra::marrangeGrob(
  grobs = alluvial_plot_list,
  nrow = 1,
  ncol = 1,
  top = NULL
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    glue::glue("{study}_deconvolution_alluvials_{today}.pdf")
  ),
  plot = deconv_alluvial_plot,
  width = 12,
  height = 8
)
