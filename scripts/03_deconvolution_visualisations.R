# ==========================================================================
# 0. Libraries and Filename settings ----
# ==========================================================================
study <- "suppinger"
today <- "2026-07-11"
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
# 2. Deconvolution stacked bar plots ----
# ==========================================================================

time_points <- unique(deconv_results$time_point)
time_points <- time_points[order(as.numeric(sub("h$", "", time_points)))]

cell_type_colours <- c(
  "Anterior primitive streak/Def. endoderm" = "#faa38a",
  "Caudal epiblast" = "#56d312",
  "Caudal epiblast/primitive streak" = "#ffc36a",
  "Caudal mesoderm" = "#01ef92",
  "Cd63+ ectoderm-like artefact" = "#e9b000",
  "Ectopic pluripotency" = "#01d9bd",
  "Epiblast" = "#bdfe0b",
  "Epiblast/primitive streak" = "#70cb94",
  "Exiting naïve pluripotency" = "#efff4e",
  "Gut" = "#8affc4",
  "Hemogenic endothelium" = "#cfba1d",
  "Naïve pluripotency" = "#35d365",
  "Neuromesodermal progenitors" = "#edff9b",
  "Paraxial mesoderm" = "#53d240",
  "Pre-somitic mesoderm" = "#b8dfa2",
  "Primitive streak" = "#9ec72a",
  "Somite" = "#78cc6e",
  "Somite differentiation front" = "#baff73",
  "Zscan4+ Artefact" = "#a5c54a"
)

deconv_algorithms <- unique(deconv_results$deconvolution_algorithm)
deconv_algorithms <- c(
  intersect(names(algorithm_colours), deconv_algorithms),
  sort(setdiff(deconv_algorithms, names(algorithm_colours)))
)

treatment_labels <- c(
  controls = "Controls",
  early_treatment = "Early treatment"
)
treatment_shapes <- c(
  controls = 16L,
  early_treatment = 17L
)

deconv_results_annotated <- deconv_results |>
  dplyr::mutate(
    treatment = dplyr::case_when(
      grepl("_CTL_", .data$sample_id) ~ "controls",
      grepl("_Early_", .data$sample_id) ~ "early_treatment",
      TRUE ~ NA_character_
    ),
    time_point = factor(.data$time_point, levels = time_points),
    deconvolution_algorithm = factor(
      .data$deconvolution_algorithm,
      levels = deconv_algorithms
    ),
    treatment = factor(.data$treatment, levels = names(treatment_labels))
  )

facet_labeller <- ggplot2::labeller(
  treatment = treatment_labels,
  deconvolution_algorithm = ggplot2::label_value
)

cell_type_levels <- deconv_results_annotated |>
  dplyr::group_by(.data$cell_type) |>
  dplyr::summarise(
    overall = mean(.data$cell_proportion, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$overall) |>
  dplyr::pull(.data$cell_type)

# ==========================================================================
# 3. Deconvolution boxplot clouds ----
# ==========================================================================

plot_list <- purrr::map(
  time_points,
  \(tp) {
    plot_data <- deconv_results_annotated |>
      dplyr::filter(.data$time_point == tp, !is.na(.data$treatment))

    cell_type_levels_tp <- intersect(
      cell_type_levels,
      unique(plot_data$cell_type)
    )

    plot_data <- plot_data |>
      dplyr::mutate(
        cell_type = factor(.data$cell_type, levels = cell_type_levels_tp)
      )

    dodge_width <- 0.82

    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$cell_type,
        y = .data$cell_proportion,
        color = .data$deconvolution_algorithm,
        fill = .data$deconvolution_algorithm,
        shape = .data$treatment,
        group = interaction(
          .data$cell_type,
          .data$deconvolution_algorithm,
          .data$treatment
        )
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
        position = ggplot2::position_dodge(width = dodge_width),
        width = 0.35,
        alpha = 0.45,
        outlier.shape = NA,
        linewidth = 0.35
      ) +
      ggplot2::geom_point(
        position = ggplot2::position_jitterdodge(
          jitter.width = 0.1,
          dodge.width = dodge_width
        ),
        size = 2.8,
        alpha = 0.85
      ) +
      ggplot2::stat_summary(
        ggplot2::aes(
          group = interaction(
            .data$cell_type,
            .data$deconvolution_algorithm,
            .data$treatment
          ),
          shape = .data$treatment
        ),
        fun = mean,
        geom = "point",
        size = 2.8,
        fill = "white",
        color = "black",
        position = ggplot2::position_dodge(width = dodge_width),
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
      ggplot2::scale_shape_manual(
        values = treatment_shapes,
        labels = treatment_labels,
        name = "Treatment"
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
        caption = "Points show biological replicates; outlined shapes show means."
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
# 4. Alluvial plots (one page per algorithm, single PDF) ----
# ============================================================================

alluvial_plot_list <- purrr::map(
  deconv_algorithms,
  \(algo) {
    algo_agg <- deconv_results_annotated |>
      dplyr::filter(
        .data$deconvolution_algorithm == algo,
        !is.na(.data$treatment)
      ) |>
      dplyr::group_by(.data$treatment, .data$time_point, .data$cell_type) |>
      dplyr::summarise(
        cell_proportion = mean(.data$cell_proportion, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::group_by(.data$treatment, .data$time_point) |>
      dplyr::mutate(
        cell_proportion = .data$cell_proportion / sum(.data$cell_proportion)
      ) |>
      dplyr::ungroup()

    plot_agg <- dplyr::mutate(
      algo_agg,
      time_point = factor(.data$time_point, levels = time_points),
      cell_type = factor(.data$cell_type, levels = cell_type_levels),
      treatment = factor(.data$treatment, levels = names(treatment_labels))
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
      ggplot2::facet_grid(
        treatment ~ .,
        labeller = facet_labeller
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent,
        limits = c(0, 1),
        expand = c(0, 0)
      ) +
      ggplot2::scale_fill_manual(
        values = cell_type_colours,
        name = "Cell type"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        panel.grid.major.x = ggplot2::element_blank(),
        strip.text.y = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::labs(
        title = glue::glue("{algo} deconvolution - {study}"),
        subtitle = paste(
          "Mean cell proportion per time point across three biological",
          "replicates, normalised to 100% within each treatment and time point."
        ),
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
  height = 10
)
