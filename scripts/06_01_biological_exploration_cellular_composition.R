# ==========================================================================
# 0. Libraries and filename settings ----
# ==========================================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(sccomp)
library(speckle)

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
# 3. Empirical log2 fold-change from sample compositions (Seurat meta) ----
# ==========================================================================
# Dots on the plot below use a simple sample-level estimate; bars come from
# scCODA (section 4). Both are on the same scale: log2 ratio-of-ratios with
# pluripotent as the compositional reference (matching pertpy / scCODA).

luque_120h <- readRDS(file.path(intermediate_dir, "luque_120h_seurat.rds"))

# Tiny pseudocount so log2 is defined when a cell type is absent in one arm.
eps <- 1e-8

# Biological reference chosen in the Python scCODA run (not data-driven).
reference_cell_type <- "pluripotent"

# Per-sample cell-type counts (one row per Sample.barcode).
sample_relative_longer <- luque_120h@meta.data |>
  dplyr::count(
    .data$Sample.barcode,
    .data$Morphotype,
    .data$luque_cluster_annotation,
    name = "n"
  ) 

sample_relative <- sample_relative_longer |>
  tidyr::pivot_wider(
    names_from = luque_cluster_annotation,
    values_from = n,
    values_fill = 0
  )
cell_types <- setdiff(
  colnames(sample_relative),
  c("Sample.barcode", "Morphotype")
)

# Compositional closure: each sample sums to 1 (relative abundance).
sample_relative <- sample_relative |>
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(cell_types),
      ~ .x / rowSums(dplyr::pick(dplyr::all_of(cell_types)), na.rm = TRUE)
    )
  )

# Mean composition per morphotype (unweighted across samples).
condition_means <- sample_relative |>
  dplyr::group_by(.data$Morphotype) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(cell_types), mean),
    .groups = "drop"
  )

nb_means <- condition_means |>
  dplyr::filter(.data$Morphotype == "neural_bias") |>
  dplyr::select(dplyr::all_of(cell_types)) |>
  as.list()
tls_means <- condition_means |>
  dplyr::filter(.data$Morphotype == "TLS") |>
  dplyr::select(dplyr::all_of(cell_types)) |>
  as.list()

# Ratio of ratios vs pluripotent:
#   log2( [p_nb(ct)/p_nb(pluripotent)] / [p_TLS(ct)/p_TLS(pluripotent)] )
# Positive => enrichment in neural_bias relative to TLS (on log2 scale).
empirical_log2_fc <- tibble::tibble(
  `Cell Type` = cell_types,
  log2_neural_bias_vs_TLS = log2(
    ((unlist(nb_means[cell_types]) + eps) /
      (nb_means[[reference_cell_type]] + eps)) /
      ((unlist(tls_means[cell_types]) + eps) /
        (tls_means[[reference_cell_type]] + eps))
  )
)

# ==========================================================================
# 4. ggplot2 — log2 fold-change vs TLS (scCODA effects) ----
# ==========================================================================
# Bars = Bayesian 95% equal-tailed intervals (ETI) from pertpy scCODA.
# Dots  = empirical sample-mean estimate from section 3.

plot_dir <- "outputs/compositional-analysis"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
cell_type_colours <- stats::setNames(
  c(
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A"
  ),
  c(
    "neural",
    "neuromesodermal progenitors",
    "somitic",
    "unknown"
  )
)

effects_plot <- readr::read_csv(
  file.path(
    "outputs",
    "biological-exploration",
    "pertpy",
    "luque_GSE250136_120h_sccoda_effects_extended.csv"
  ),
  show_col_types = FALSE
) |>
  # Reference cell type has zero effect by construction in scCODA.
  dplyr::filter(.data$`Cell Type` != "pluripotent") |>
  dplyr::mutate(
    cell_type = factor(
      .data$`Cell Type`,
      levels = c("unknown", "somitic", "neuromesodermal progenitors", "neural")
    ),
    # * marks ETI excluding 0 (directional change). Distinct from spike-and-
    # slab "credible" at FDR 5%, which was non-significant for all types here.
    significant = .data$relative_log2_low > 0 | .data$relative_log2_high < 0,
    # Bar endpoints: 95% ETI on log2 scale (natural-log ETI / log(2) in Python).
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

# Symmetric x-axis from the widest 95% ETI bar (not from empirical dots).
x_lim <- max(abs(c(effects_plot$xmin, effects_plot$xmax)), na.rm = TRUE) + 0.15

n_types <- nrow(effects_plot)
y_top <- n_types + 0.65

sccoda_l2fc_plot <- ggplot2::ggplot(effects_plot) +
  # Zero = no morphotype difference on the log2 ratio-of-ratios scale.
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.4
  ) +
  # Horizontal bars: scCODA 95% ETI (uncertainty from the Bayesian model).
  ggplot2::geom_rect(
    ggplot2::aes(
      xmin = .data$xmin,
      xmax = .data$xmax,
      ymin = .data$y_idx - 0.35,
      ymax = .data$y_idx + 0.35,
      fill = .data$cell_type
    ),
    alpha = 0.75,
    colour = NA
  ) +
  # Empirical point estimate from section 3 (sample-mean composition).
  ggplot2::geom_point(
    ggplot2::aes(
      x = .data$log2_neural_bias_vs_TLS,
      y = .data$y_idx,
      colour = .data$cell_type
    ),
    size = 2.8
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = .data$star_x,
      y = .data$y_idx,
      label = .data$star
    ),
    size = 5,
    hjust = 0.5,
    vjust = 0.5,
    colour = "grey20"
  ) +
  # Direction labels: left = TLS-enriched, right = neural_bias-enriched.
  ggplot2::annotate(
    "text",
    x = -x_lim,
    y = y_top,
    label = "Relative increase in TLS",
    hjust = 0,
    size = 3.8,
    fontface = "plain"
  ) +
  ggplot2::annotate(
    "text",
    x = x_lim,
    y = y_top,
    label = "Relative increase in neural_bias",
    hjust = 1,
    size = 3.8,
    fontface = "plain"
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
    data = effects_plot,
    ggplot2::aes(
      x = -Inf,
      y = .data$y_idx,
      label = .data$cell_type,
      colour = .data$cell_type
    ),
    hjust = 1.08,
    vjust = 0.5,
    size = 3.5,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(values = cell_type_colours) +
  ggplot2::scale_colour_manual(values = cell_type_colours) +
  ggplot2::scale_x_continuous(
    limits = c(-x_lim, x_lim),
    breaks = scales::pretty_breaks(n = 5),
    expand = ggplot2::expansion(mult = c(0.03, 0.03))
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq_len(n_types),
    labels = NULL,
    limits = c(0.5, y_top + 0.35)
  ) +
  # clip = "off" lets cell-type labels sit left of the panel (x = -Inf).
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    x = expression(log[2] ~ "(ratio of ratios vs pluripotent)"),
    y = NULL,
    caption = paste(
      "Luque GSE250136, 120h scCODA (pertpy): neural_bias vs TLS, pluripotent reference.",
      "Bars = 95% ETI on log2-scaled compositional effects;",
      "* = interval excludes 0 (no spike-and-slab selection at FDR 5%)."
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
    plot.margin = ggplot2::margin(t = 20, r = 15, b = 10, l = 200)
  )

ggplot2::ggsave(
  filename = file.path(plot_dir, "sccoda_log2fc_morphotype.pdf"),
  plot = sccoda_l2fc_plot,
  width = 10,
  height = 5
)

message(
  "Saved ",
  file.path(plot_dir, "sccoda_log2fc_morphotype.pdf")
)

# ==========================================================================
# 5. sccomp — differential composition (neural_bias vs TLS) ----
# ==========================================================================
# Per-organoid cell-type counts (Sample.barcode); TLS is the compositional
# reference morphotype. Effects are on the logit scale (c_effect).

counts_long <- sample_relative_longer |>
  dplyr::rename(n_cells = n) |>
  dplyr::filter(.data$Morphotype %in% c("TLS", "neural_bias")) |>
  dplyr::mutate(
    Morphotype = stats::relevel(factor(.data$Morphotype), ref = "TLS"),
    Sample.barcode = factor(.data$Sample.barcode),
    luque_cluster_annotation = factor(.data$luque_cluster_annotation)
  )

sccomp_result <- counts_long |>
  sccomp_estimate(
    formula_composition = ~ Morphotype,
    formula_variability = ~ 1,
    sample = "Sample.barcode",
    cell_group = "luque_cluster_annotation",
    abundance = "n_cells",
    cores = 4,
    verbose = FALSE
  ) |>
  sccomp_test()

sccomp_neural_bias_effects <- sccomp_result |>
  dplyr::filter(.data$parameter == "Morphotypeneural_bias") |>
  dplyr::select(
    luque_cluster_annotation,
    c_effect,
    c_lower,
    c_upper,
    c_pH0,
    c_FDR
  )

readr::write_csv(
  sccomp_result,
  file.path(plot_dir, "luque_120h_sccomp.csv")
)

sccomp_boxplot <- sccomp_result |>
  sccomp_boxplot(factor = "Morphotype")

ggplot2::ggsave(
  filename = file.path(plot_dir, "sccomp_boxplot_morphotype.pdf"),
  plot = sccomp_boxplot,
  width = 10,
  height = 6
)
