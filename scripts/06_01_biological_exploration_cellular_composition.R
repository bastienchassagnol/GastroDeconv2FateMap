# ==========================================================================
# 0. Libraries and filename settings ----
# ==========================================================================

library(dplyr)
library(flextable)
library(forcats)
library(ggdist)
library(ggplot2)
library(glue)
library(patchwork)
library(readr)
library(scales)
library(sccomp)
library(tidyr)
library(tibble)

study <- "luque"
intermediate_dir <- "data/intermediate"


# ==========================================================================
# 1. Load merged object, realign barcodes, subset 120h ----
# ==========================================================================

# h5ad_path <- file.path(intermediate_dir, "luque_120h_sccoda.h5ad")
# 1.1 Extract 120h, and raw counts assay----
# luque_merged <- readRDS(
#   file = "./data/intermediate/luque_single_cell_merged_2026-06-07.rds"
# )

# # Merged object colnames are numeric indices; RNA assay uses cell barcodes.
# counts <- SeuratObject::GetAssayData(
#   object = luque_merged,
#   assay = "RNA",
#   layer = "counts"
# )
# meta <- luque_merged@meta.data
# rownames(meta) <- colnames(counts)

# luque_merged <- CreateSeuratObject(
#   counts = counts,
#   meta.data = meta
# )

# # Subset on the last time point only
# luque_120h <- subset(luque_merged, subset = timepoint == "120h")

# cluster_mapping <- readr::read_csv(
#   "./data/intermediate/mapping_seurat_villaronga.csv",
#   show_col_types = FALSE
# ) |>
#   dplyr::mutate(
#     seurat_clusters = as.character(.data$seurat_clusters)
#   )

# meta_with_annotation <- luque_120h@meta.data |>
#   tibble::rownames_to_column("cell_barcode") |>
#   dplyr::mutate(
#     seurat_clusters = as.character(.data$seurat_clusters)
#   ) |>
#   dplyr::select(-dplyr::any_of("luque_cluster_annotation")) |>
#   dplyr::inner_join(
#     cluster_mapping,
#     by = c("timepoint", "seurat_clusters")
#   ) |>
#   dplyr::filter(
#     !is.na(.data$Sample.barcode),
#     !is.na(.data$luque_cluster_annotation),
#     !is.na(.data$Morphotype),
#     .data$Morphotype %in% c("TLS", "neural_bias")
#   )

# matched_cells <- meta_with_annotation$cell_barcode
# luque_merged <- subset(luque_merged, cells = matched_cells)

# # 1.2 export phenotype data  only----

# luque_120h@meta.data <- luque_120h@meta.data[,
#   c(
#     "Sample.barcode",
#     "luque_cluster_annotation",
#     "Morphotype"
#   ),
#   drop = FALSE
# ]

# saveRDS(
#   luque_120h,
#   file = file.path(intermediate_dir, "luque_120h_seurat.rds")
# )

# # 1.3 export cell-level AnnData for Python scCODA (via AnnDataR) ----
# anndataR::write_h5ad(
#   object = luque_merged,
#   path = h5ad_path,
#   mode = "w"
# )

# message(
#   "Exported ",
#   ncol(luque_merged),
#   " cells to ",
#   h5ad_path,
#   " and sample counts to ",
#   counts_path
# )

# ==========================================================================
# 3. Cellular composition — counts, summaries, visualisations ----
# ==========================================================================

plot_dir <- "outputs/compositional-analysis"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

luque_120h <- readRDS(file.path(intermediate_dir, "luque_120h_seurat.rds"))

morphotype_colours <- stats::setNames(
  c("#377EB8", "#E41A1C"),
  c("TLS", "neural_bias")
)
cell_type_levels <- c(
  "pluripotent",
  "neuromesodermal progenitors",
  "somitic",
  "neural",
  "unknown"
)
cell_type_colours <- stats::setNames(
  c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"),
  cell_type_levels
)
reference_cell_type <- "pluripotent"
eps <- 1e-8

# Per-organoid cell-type counts (shared by all composition plots and models).
composition_counts <- luque_120h@meta.data |>
  dplyr::count(
    .data$Sample.barcode,
    .data$Morphotype,
    .data$luque_cluster_annotation,
    name = "n_cells"
  ) |>
  dplyr::filter(.data$Morphotype %in% c("TLS", "neural_bias")) |>
  dplyr::mutate(
    Morphotype = stats::relevel(factor(.data$Morphotype), ref = "TLS"),
    Sample.barcode = factor(.data$Sample.barcode),
    luque_cluster_annotation = factor(
      .data$luque_cluster_annotation,
      levels = cell_type_levels
    )
  )

# 3a. Stacked bar — mean organoid composition per morphotype (sums to 100%) ----
# Per organoid (Sample.barcode), counts are summed by cell type. For each
# morphotype, we take the unweighted mean of those counts across organoids
# (mean_n). Relative abundance (pct) is then mean_n / sum(mean_n) within each
# morphotype, so both bars sum to 100%. Cell types are reordered with
# fct_reorder (least abundant first in factor levels = top of bar). ggplot2
# stacks in reverse level order (last level at the bottom), so ymin/ymax are
# computed with cumsum after arranging by descending factor level.
stack_summary <- composition_counts |>
  dplyr::group_by(.data$Morphotype, .data$luque_cluster_annotation) |>
  dplyr::summarise(mean_n = mean(.data$n_cells), .groups = "drop") |>
  dplyr::mutate(
    luque_cluster_annotation = forcats::fct_reorder(
      .data$luque_cluster_annotation,
      .data$mean_n,
      .fun = sum,
      .desc = FALSE
    )
  ) |>
  dplyr::group_by(.data$Morphotype) |>
  dplyr::arrange(desc(as.integer(.data$luque_cluster_annotation)), .by_group = TRUE) |>
  dplyr::mutate(
    pct = .data$mean_n / sum(.data$mean_n) * 100,
    ymax = cumsum(.data$pct),
    ymin = .data$ymax - .data$pct,
    y_label = (.data$ymin + .data$ymax) / 2,
    label = glue::glue(
      "{round(.data$mean_n, 0)} ({sprintf('%.1f', .data$pct)}%)"
    )
  ) |>
  dplyr::ungroup()

stacked_bar_plot <- ggplot2::ggplot(
  stack_summary,
  ggplot2::aes(
    x = .data$Morphotype,
    y = .data$pct,
    fill = .data$luque_cluster_annotation
  )
) +
  ggplot2::geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
  ggplot2::geom_label(
    ggplot2::aes(y = .data$y_label, label = .data$label),
    position = ggplot2::position_identity(),
    show.legend = FALSE,
    size = 2.8,
    colour = "grey10",
    fill = "white",
    linewidth = 0.2,
    label.padding = ggplot2::unit(0.15, "lines")
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0, 0.08))
  ) +
  ggplot2::scale_fill_manual(values = cell_type_colours, name = "Cell type") +
  ggplot2::labs(
    x = "Morphotype",
    y = "Relative abundance (%)"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.title = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  filename = file.path(plot_dir, "stacked_bar_morphotype_composition.pdf"),
  plot = stacked_bar_plot,
  width = 7,
  height = 6,
  device = cairo_pdf
)

# 3b. Raincloud — per-organoid cell-type frequencies by morphotype ----
raincloud_dodge <- 0.35
morphotype_x_offset <- stats::setNames(
  c(-raincloud_dodge / 2, raincloud_dodge / 2),
  c("TLS", "neural_bias")
)

set.seed(42)
raincloud_data <- composition_counts |>
  dplyr::group_by(.data$Sample.barcode, .data$Morphotype) |>
  dplyr::mutate(frequency = .data$n_cells / sum(.data$n_cells)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    x_id = as.numeric(.data$luque_cluster_annotation),
    x_center = .data$x_id +
      morphotype_x_offset[as.character(.data$Morphotype)],
    x_point = .data$x_center - 0.1 +
      stats::runif(dplyr::n(), min = -0.02, max = 0.02)
  )

raincloud_plot <- ggplot2::ggplot(
  raincloud_data,
  ggplot2::aes(
    x = .data$x_center,
    y = .data$frequency,
    fill = .data$Morphotype
  )
) +
  ggdist::stat_halfeye(
    orientation = "vertical",
    side = "right",
    density = "unbounded",
    limits = c(0, 1),
    trim = FALSE,
    expand = TRUE,
    slab_alpha = 0.45,
    slab_colour = NA,
    interval_alpha = 0,
    point_alpha = 0,
    scale = 0.65,
    normalize = "groups"
  ) +
  ggplot2::geom_boxplot(
    ggplot2::aes(
      group = interaction(
        .data$luque_cluster_annotation,
        .data$Morphotype
      )
    ),
    width = 0.12,
    outlier.shape = NA,
    alpha = 0.75,
    colour = "grey20",
    linewidth = 0.35
  ) +
  ggplot2::geom_point(
    ggplot2::aes(x = .data$x_point, y = .data$frequency, colour = Morphotype),
    size = 1.1,
    alpha = 0.75,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_fill_manual(values = morphotype_colours) +
  ggplot2::scale_colour_manual(values = morphotype_colours) +
  ggplot2::scale_x_continuous(
    breaks = seq_along(cell_type_levels),
    labels = cell_type_levels,
    expand = ggplot2::expansion(mult = c(0.08, 0.08))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    x = "Cell type",
    y = "Relative abundance",
    fill = "Morphotype",
    colour = "Morphotype",
    caption = paste(
      "Luque GSE250136, 120h: per-organoid cell-type frequencies.",
      "Raincloud — scatter (left), boxplot (centre), half-eye (right)."
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
    legend.position = "top",
    plot.caption = ggplot2::element_text(
      hjust = 0,
      size = 8,
      colour = "grey30",
      margin = ggplot2::margin(t = 8)
    )
  )

ggplot2::ggsave(
  filename = file.path(plot_dir, "raincloud_morphotype_composition.pdf"),
  plot = raincloud_plot,
  width = 11,
  height = 6,
  device = cairo_pdf
)

# 3c. Empirical log2 fold-change (sample-mean; pluripotent reference) ----
condition_means <- composition_counts |>
  tidyr::pivot_wider(
    names_from = luque_cluster_annotation,
    values_from = n_cells,
    values_fill = 0
  ) |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(cell_type_levels),
      ~ .x / rowSums(dplyr::pick(dplyr::all_of(cell_type_levels)), na.rm = TRUE)
    )
  ) |>
  dplyr::group_by(.data$Morphotype) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(cell_type_levels), mean),
    .groups = "drop"
  )

nb_means <- condition_means |>
  dplyr::filter(.data$Morphotype == "neural_bias") |>
  dplyr::select(dplyr::all_of(cell_type_levels)) |>
  as.list()
tls_means <- condition_means |>
  dplyr::filter(.data$Morphotype == "TLS") |>
  dplyr::select(dplyr::all_of(cell_type_levels)) |>
  as.list()

empirical_log2_fc <- tibble::tibble(
  `Cell Type` = cell_type_levels,
  log2_neural_bias_vs_TLS = log2(
    ((unlist(nb_means[cell_type_levels]) + eps) /
      (nb_means[[reference_cell_type]] + eps)) /
      ((unlist(tls_means[cell_type_levels]) + eps) /
        (tls_means[[reference_cell_type]] + eps))
  )
)

# ==========================================================================
# 4. scCODA log2 fold-change vs TLS (forest plot + effects table) ----
# ==========================================================================

sccoda_colours <- cell_type_colours[
  c("neural", "neuromesodermal progenitors", "somitic", "unknown")
]

sccoda_credible <- readr::read_csv(
  file.path(
    "outputs",
    "biological-exploration",
    "pertpy",
    "luque_GSE250136_120h_sccoda_credible_effects.csv"
  ),
  show_col_types = FALSE
)

sccoda_diag <- readr::read_csv(
  file.path(
    "outputs",
    "biological-exploration",
    "pertpy",
    "luque_GSE250136_120h_sccoda_mcmc_diagnostics.csv"
  ),
  show_col_types = FALSE
) |>
  dplyr::filter(grepl("^beta\\[", .data$parameter)) |>
  dplyr::mutate(
    `Cell Type` = sub(".*, ([^]]+)\\]$", "\\1", .data$parameter)
  ) |>
  dplyr::select("Cell Type", "ess_bulk", "ess_tail")

effects_plot <- readr::read_csv(
  file.path(
    "outputs",
    "biological-exploration",
    "pertpy",
    "luque_GSE250136_120h_sccoda_effects_extended.csv"
  ),
  show_col_types = FALSE
) |>
  dplyr::filter(.data$`Cell Type` != "pluripotent") |>
  dplyr::left_join(
    sccoda_credible |>
      dplyr::select("Cell Type", "credible"),
    by = "Cell Type"
  ) |>
  dplyr::left_join(sccoda_diag, by = "Cell Type") |>
  dplyr::mutate(
    cell_type = factor(
      .data$`Cell Type`,
      levels = c("unknown", "somitic", "neuromesodermal progenitors", "neural")
    ),
    significant = .data$relative_log2_low > 0 | .data$relative_log2_high < 0,
    xmin = pmin(.data$relative_log2_low, .data$relative_log2_high),
    xmax = pmax(.data$relative_log2_low, .data$relative_log2_high),
    star = ifelse(.data$significant, "*", ""),
    star_x = ifelse(
      .data$relative_log2_high < 0,
      .data$xmin - 0.06,
      .data$xmax + 0.06
    ),
    y_idx = as.numeric(.data$cell_type)
  ) |>
  dplyr::left_join(empirical_log2_fc, by = "Cell Type")

effects_table <- effects_plot |>
  dplyr::transmute(
    `Cell Type` = as.character(.data$cell_type),
    `Log2 FC` = .data$log2_neural_bias_vs_TLS,
    `ETI Low` = .data$relative_log2_low,
    `ETI High` = .data$relative_log2_high,
    `Exp Sample` = .data$`Expected Sample`,
    `Incl Prob` = .data$`Inclusion probability`,
    `Credible` = ifelse(.data$credible, "Yes", "No"),
    `ESS Bulk` = .data$ess_bulk,
    `ESS Tail` = .data$ess_tail
  )

effects_flex <- effects_table |>
  flextable::flextable() |>
  flextable::colformat_double(
    j = c("Log2 FC", "ETI Low", "ETI High"),
    digits = 3
  ) |>
  flextable::colformat_double(j = "Exp Sample", digits = 1) |>
  flextable::colformat_double(j = "Incl Prob", digits = 3) |>
  flextable::colformat_double(j = c("ESS Bulk", "ESS Tail"), digits = 0) |>
  flextable::align(align = "center", part = "all") |>
  flextable::align(j = 1, align = "left", part = "body") |>
  flextable::fontsize(size = 8, part = "all") |>
  flextable::padding(padding = 2, part = "all") |>
  flextable::autofit() |>
  flextable::theme_vanilla() |>
  flextable::footnote(
    i = 1,
    j = 2,
    ref_symbols = "a",
    value = flextable::as_paragraph(
      "Log2 FC: empirical sample-mean log2 ratio-of-ratios vs pluripotent."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1,
    j = 3,
    ref_symbols = "b",
    value = flextable::as_paragraph(
      "ETI Low: lower bound of the 95% equal-tailed credible interval (log2)."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1,
    j = 4,
    ref_symbols = "c",
    value = flextable::as_paragraph(
      "ETI High: upper bound of the 95% equal-tailed credible interval (log2)."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1,
    j = 5,
    ref_symbols = "d",
    value = flextable::as_paragraph(
      "Exp Sample: posterior expected cell count."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1,
    j = 6,
    ref_symbols = "e",
    value = flextable::as_paragraph(
      "Incl Prob: spike-and-slab posterior inclusion probability."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1,
    j = 7,
    ref_symbols = "f",
    value = flextable::as_paragraph(
      "Credible: effect selected at FDR 5%."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1,
    j = 8,
    ref_symbols = "g",
    value = flextable::as_paragraph(
      "ESS Bulk: effective sample size for the main posterior mass; higher is better."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1,
    j = 9,
    ref_symbols = "h",
    value = flextable::as_paragraph(
      "ESS Tail: effective sample size for posterior tails; higher is better."
    ),
    part = "header"
  ) |>
  flextable::fontsize(size = 7, part = "footer")

x_lim <- max(abs(c(effects_plot$xmin, effects_plot$xmax)), na.rm = TRUE) + 0.2
y_top <- nrow(effects_plot) + 0.85

sccoda_l2fc_plot <- ggplot2::ggplot(effects_plot) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.4
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = .data$xmin,
      xmax = .data$xmax,
      y = .data$y_idx,
      colour = .data$cell_type
    ),
    height = 0.3,
    linewidth = 0.85
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = .data$log2_neural_bias_vs_TLS,
      y = .data$y_idx,
      colour = .data$cell_type
    ),
    size = 3.6
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = .data$star_x, y = .data$y_idx, label = .data$star),
    size = 5,
    colour = "grey20"
  ) +
  ggplot2::annotate(
    "text",
    x = -x_lim * 0.92,
    y = y_top,
    label = "TLS enriched",
    hjust = 0,
    size = 3.5
  ) +
  ggplot2::annotate(
    "text",
    x = x_lim * 0.92,
    y = y_top,
    label = "neural_bias enriched",
    hjust = 1,
    size = 3.5
  ) +
  ggplot2::annotate(
    "segment",
    x = -0.05,
    xend = -0.35,
    y = y_top,
    yend = y_top,
    arrow = grid::arrow(length = grid::unit(0.15, "cm"), ends = "last")
  ) +
  ggplot2::annotate(
    "segment",
    x = 0.05,
    xend = 0.35,
    y = y_top,
    yend = y_top,
    arrow = grid::arrow(length = grid::unit(0.15, "cm"), ends = "last")
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = -Inf,
      y = .data$y_idx,
      label = .data$cell_type,
      colour = .data$cell_type
    ),
    hjust = 1.08,
    size = 3.5,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = sccoda_colours) +
  ggplot2::scale_x_continuous(
    limits = c(-x_lim, x_lim),
    breaks = scales::pretty_breaks(n = 5),
    expand = ggplot2::expansion(mult = c(0.03, 0.03))
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq_len(nrow(effects_plot)),
    labels = NULL,
    limits = c(0.5, y_top + 0.45)
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    x = expression(log[2] ~ "(ratio of ratios vs pluripotent)"),
    y = NULL,
    caption = paste(
      "Luque GSE250136, 120h scCODA (pertpy): neural_bias vs TLS.",
      "Effect = log2 ratio-of-ratios relative to pluripotent (reference cell type):",
      "morphotype-associated change for each cell type beyond that of the reference,",
      "i.e. a compositional difference-in-differences (Eitzinger et al., scCODA).",
      "Whiskers = 95% ETI on log2 scale; dots = empirical sample-mean estimate;",
      "* = ETI excludes 0; credible = spike-and-slab selection at FDR 5%."
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(
      hjust = 0,
      size = 8,
      colour = "grey30",
      margin = ggplot2::margin(t = 10)
    ),
    plot.margin = ggplot2::margin(t = 28, r = 15, b = 10, l = 200)
  )

sccoda_table_panel <- ggplot2::ggplot() +
  ggplot2::theme_void() +
  ggplot2::coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  patchwork::inset_element(
    p = flextable::gen_grob(effects_flex),
    left = 0,
    right = 1,
    bottom = 0.1,
    top = 0.9,
    align_to = "full"
  )

sccoda_figure <- sccoda_l2fc_plot |
  sccoda_table_panel +
  patchwork::plot_layout(widths = c(1.4, 1.1))

ggplot2::ggsave(
  filename = file.path(plot_dir, "sccoda_log2fc_morphotype.pdf"),
  plot = sccoda_figure,
  width = 22,
  height = 6.5,
  device = cairo_pdf
)

# ==========================================================================
# 5. sccomp — differential composition (neural_bias vs TLS) ----
# ==========================================================================

sccomp_result <- composition_counts |>
  sccomp_estimate(
    formula_composition = ~Morphotype,
    formula_variability = ~1,
    sample = "Sample.barcode",
    cell_group = "luque_cluster_annotation",
    abundance = "n_cells",
    cores = 4,
    verbose = FALSE
  ) |>
  sccomp_test()

readr::write_csv(
  sccomp_result,
  file.path(plot_dir, "luque_120h_sccomp.csv")
)

ggplot2::ggsave(
  filename = file.path(plot_dir, "sccomp_boxplot_morphotype.pdf"),
  plot = sccomp_boxplot(sccomp_result, factor = "Morphotype"),
  width = 10,
  height = 6
)

# 5b. sccomp forest plot + effects table (aligned with scCODA layout) ----
# Empirical Log2 FC: log2(mean proportion in neural_bias / mean proportion in TLS)
# per cell type (no compositional reference cell type; all five types shown).
empirical_log2_sccomp <- tibble::tibble(
  `Cell Type` = cell_type_levels,
  log2_neural_bias_vs_TLS = log2(
    (unlist(nb_means[cell_type_levels]) + eps) /
      (unlist(tls_means[cell_type_levels]) + eps)
  )
)

# Model ETI on log2 scale: sccomp estimates morphotype contrasts on the logit
# scale (c_effect, c_lower, c_upper). We re-express these as log2 proportion
# fold-changes at the empirical TLS mean proportion for each cell type:
#   p_ref = TLS mean proportion (closure),
#   logit_ref = qlogis(p_ref),
#   p_contrast = plogis(logit_ref + c_*),
#   log2 FC = log2(p_contrast / p_ref).
tls_props <- tls_means[cell_type_levels]
p_ref_safe <- pmin(pmax(unlist(tls_props), eps), 1 - eps)
logit_ref <- stats::qlogis(p_ref_safe)

sccomp_mean_n <- composition_counts |>
  dplyr::group_by(.data$luque_cluster_annotation) |>
  dplyr::summarise(`Exp Sample` = mean(.data$n_cells), .groups = "drop") |>
  dplyr::rename(`Cell Type` = luque_cluster_annotation)

sccomp_effects <- sccomp_result |>
  dplyr::filter(.data$parameter == "Morphotypeneural_bias") |>
  dplyr::rename(`Cell Type` = luque_cluster_annotation) |>
  dplyr::left_join(sccomp_mean_n, by = "Cell Type") |>
  dplyr::mutate(
    p_low = stats::plogis(
      logit_ref[match(as.character(.data$`Cell Type`), cell_type_levels)] +
        .data$c_lower
    ),
    p_eff = stats::plogis(
      logit_ref[match(as.character(.data$`Cell Type`), cell_type_levels)] +
        .data$c_effect
    ),
    p_high = stats::plogis(
      logit_ref[match(as.character(.data$`Cell Type`), cell_type_levels)] +
        .data$c_upper
    ),
    relative_log2_low = log2(
      .data$p_low /
        p_ref_safe[match(as.character(.data$`Cell Type`), cell_type_levels)]
    ),
    relative_log2_effect = log2(
      .data$p_eff /
        p_ref_safe[match(as.character(.data$`Cell Type`), cell_type_levels)]
    ),
    relative_log2_high = log2(
      .data$p_high /
        p_ref_safe[match(as.character(.data$`Cell Type`), cell_type_levels)]
    ),
    `Incl Prob` = 1 - .data$c_pH0
  )

sccomp_plot_data <- sccomp_effects |>
  dplyr::mutate(
    cell_type = factor(
      .data$`Cell Type`,
      levels = cell_type_levels
    ),
    significant = .data$relative_log2_low > 0 | .data$relative_log2_high < 0,
    xmin = pmin(.data$relative_log2_low, .data$relative_log2_high),
    xmax = pmax(.data$relative_log2_low, .data$relative_log2_high),
    star = ifelse(.data$significant, "*", ""),
    star_x = ifelse(
      .data$relative_log2_high < 0,
      .data$xmin - 0.06,
      .data$xmax + 0.06
    ),
    y_idx = as.numeric(.data$cell_type),
    credible = .data$c_FDR < 0.05
  ) |>
  dplyr::left_join(empirical_log2_sccomp, by = "Cell Type")

sccomp_table <- sccomp_plot_data |>
  dplyr::transmute(
    `Cell Type` = as.character(.data$cell_type),
    `Log2 FC` = .data$log2_neural_bias_vs_TLS,
    `ETI Low` = .data$relative_log2_low,
    `ETI High` = .data$relative_log2_high,
    `Exp Sample` = .data$`Exp Sample`,
    `Incl Prob` = .data$`Incl Prob`,
    `Credible` = ifelse(.data$credible, "Yes", "No"),
    `ESS Bulk` = .data$c_ess_bulk,
    `ESS Tail` = .data$c_ess_tail
  )

sccomp_flex <- sccomp_table |>
  flextable::flextable() |>
  flextable::colformat_double(
    j = c("Log2 FC", "ETI Low", "ETI High"),
    digits = 3
  ) |>
  flextable::colformat_double(j = "Exp Sample", digits = 1) |>
  flextable::colformat_double(j = "Incl Prob", digits = 3) |>
  flextable::colformat_double(j = c("ESS Bulk", "ESS Tail"), digits = 0) |>
  flextable::align(align = "center", part = "all") |>
  flextable::align(j = 1, align = "left", part = "body") |>
  flextable::fontsize(size = 8, part = "all") |>
  flextable::padding(padding = 2, part = "all") |>
  flextable::autofit() |>
  flextable::theme_vanilla() |>
  flextable::footnote(
    i = 1, j = 2, ref_symbols = "a",
    value = flextable::as_paragraph(
      "Log2 FC: empirical sample-mean log2(neural_bias / TLS) proportion."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1, j = 3, ref_symbols = "b",
    value = flextable::as_paragraph(
      "ETI Low: lower bound of 95% credible interval on log2 scale."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1, j = 4, ref_symbols = "c",
    value = flextable::as_paragraph(
      "ETI High: upper bound of 95% credible interval on log2 scale."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1, j = 5, ref_symbols = "d",
    value = flextable::as_paragraph(
      "Exp Sample: mean organoid cell count per cell type."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1, j = 6, ref_symbols = "e",
    value = flextable::as_paragraph(
      "Incl Prob: 1 - P(H0); posterior mass above sccomp null threshold."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1, j = 7, ref_symbols = "f",
    value = flextable::as_paragraph(
      "Credible: composition effect with FDR < 5% (sccomp_test)."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1, j = 8, ref_symbols = "g",
    value = flextable::as_paragraph(
      "ESS Bulk: effective sample size for the main posterior mass; higher is better."
    ),
    part = "header"
  ) |>
  flextable::footnote(
    i = 1, j = 9, ref_symbols = "h",
    value = flextable::as_paragraph(
      "ESS Tail: effective sample size for posterior tails; higher is better."
    ),
    part = "header"
  ) |>
  flextable::fontsize(size = 7, part = "footer")

sccomp_x_lim <- max(
  abs(c(sccomp_plot_data$xmin, sccomp_plot_data$xmax)),
  na.rm = TRUE
) + 0.2
sccomp_y_top <- nrow(sccomp_plot_data) + 0.85

sccomp_l2fc_plot <- ggplot2::ggplot(sccomp_plot_data) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.4
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(
      xmin = .data$xmin,
      xmax = .data$xmax,
      y = .data$y_idx,
      colour = .data$cell_type
    ),
    height = 0.3,
    linewidth = 0.85
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = .data$log2_neural_bias_vs_TLS,
      y = .data$y_idx,
      colour = .data$cell_type
    ),
    size = 3.6
  ) +
  ggplot2::geom_text(
    ggplot2::aes(x = .data$star_x, y = .data$y_idx, label = .data$star),
    size = 5,
    colour = "grey20"
  ) +
  ggplot2::annotate(
    "text",
    x = -sccomp_x_lim * 0.92,
    y = sccomp_y_top,
    label = "TLS enriched",
    hjust = 0,
    size = 3.5
  ) +
  ggplot2::annotate(
    "text",
    x = sccomp_x_lim * 0.92,
    y = sccomp_y_top,
    label = "neural_bias enriched",
    hjust = 1,
    size = 3.5
  ) +
  ggplot2::annotate(
    "segment",
    x = -0.05,
    xend = -0.35,
    y = sccomp_y_top,
    yend = sccomp_y_top,
    arrow = grid::arrow(length = grid::unit(0.15, "cm"), ends = "last")
  ) +
  ggplot2::annotate(
    "segment",
    x = 0.05,
    xend = 0.35,
    y = sccomp_y_top,
    yend = sccomp_y_top,
    arrow = grid::arrow(length = grid::unit(0.15, "cm"), ends = "last")
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = -Inf,
      y = .data$y_idx,
      label = .data$cell_type,
      colour = .data$cell_type
    ),
    hjust = 1.08,
    size = 3.5,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = cell_type_colours) +
  ggplot2::scale_x_continuous(
    limits = c(-sccomp_x_lim, sccomp_x_lim),
    breaks = scales::pretty_breaks(n = 5),
    expand = ggplot2::expansion(mult = c(0.03, 0.03))
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq_len(nrow(sccomp_plot_data)),
    labels = NULL,
    limits = c(0.5, sccomp_y_top + 0.45)
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    x = expression(log[2] ~ "(neural_bias / TLS proportion)"),
    y = NULL,
    caption = paste(
      "Luque GSE250136, 120h sccomp: neural_bias vs TLS (TLS reference morphotype).",
      "Dots = empirical log2(mean proportion neural_bias / mean proportion TLS).",
      "Whiskers = 95% credible interval on log2 scale, obtained by mapping",
      "sccomp logit-scale contrasts (c_lower, c_effect, c_upper) through",
      "plogis(logit(p_TLS) + c_*) and log2(p_contrast / p_TLS); no reference",
      "cell type is used (all five cell types shown).",
      "* = credible interval excludes 0; Credible = FDR < 5%."
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(
      hjust = 0,
      size = 8,
      colour = "grey30",
      margin = ggplot2::margin(t = 10)
    ),
    plot.margin = ggplot2::margin(t = 28, r = 15, b = 10, l = 200)
  )

sccomp_table_panel <- ggplot2::ggplot() +
  ggplot2::theme_void() +
  ggplot2::coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  patchwork::inset_element(
    p = flextable::gen_grob(sccomp_flex),
    left = 0,
    right = 1,
    bottom = 0.08,
    top = 0.92,
    align_to = "full"
  )

sccomp_figure <- sccomp_l2fc_plot |
  sccomp_table_panel +
  patchwork::plot_layout(widths = c(1.4, 1.1))

ggplot2::ggsave(
  filename = file.path(plot_dir, "sccomp_log2fc_morphotype.pdf"),
  plot = sccomp_figure,
  width = 22,
  height = 7,
  device = cairo_pdf
)
