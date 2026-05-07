# ==========================================================================
# 0. Libraries and Filename settings ----
# ==========================================================================

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/bulk-omics"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==========================================================================
# 1. Bulk RNA-seq data QC ----
# ==========================================================================

bulk_se <- readRDS(
  "./data/intermediate/suppinger_bulk_summarized_experiment_2026-05-06.rds"
)

dim(bulk_se)
tinytable::tt(
  bulk_se@colData |>
    as.data.frame(),
  caption = "Sample phenotype data (bulk RNA-seq)"
)

# Log-normalise raw counts; transpose so samples are rows (required by PCA)
log_counts <- SummarizedExperiment::assay(bulk_se, "counts") |>
  log1p() |>
  t()

# Keep the top most variable genes (standard bulk RNA-seq PCA criterion).
# Removes zero-variance (constant/all-zero) genes that would break scale. = TRUE
N_HVG <- 500L
gene_vars <- apply(log_counts, 2, var)
hvg_idx <- order(gene_vars, decreasing = TRUE)[
  seq_len(min(N_HVG, sum(gene_vars > 0)))
]
log_counts_hvg <- log_counts[, hvg_idx]

# Compute PCA using the top 500 most variable genes
pca_res <- prcomp(log_counts_hvg, center = TRUE, scale. = TRUE)

# Bind PC1/PC2 scores with sample metadata
col_meta <- bulk_se@colData |>
  as.data.frame()

pca_data <- as.data.frame(pca_res$x[, 1:2]) |>
  dplyr::mutate(
    time_point_id = col_meta$time_point_id,
    treatment_status = col_meta$treatment_status,
    batch_id = col_meta$batch_id
  )

var_pct <- round(pca_res$sdev[1:2]^2 / sum(pca_res$sdev^2) * 100, 1)

# Aesthetic mappings -----------------------------------------------------------
treatment_colours <- c("control" = "#3182bd", "early_treatment" = "#e6550d")
# Filled square / triangle / diamond: maximally distinct at small sizes
time_shapes <- c("48h" = 15, "72h" = 17, "96h" = 18)

pca_bulk_plot <- ggplot2::ggplot(
  pca_data,
  ggplot2::aes(x = PC1, y = PC2)
) +
  # Ellipses for treatment status only
  ggplot2::stat_ellipse(
    ggplot2::aes(color = treatment_status, fill = treatment_status),
    level = 0.95,
    geom = "polygon",
    alpha = 0.10,
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  # Points: colour = treatment status, shape = time point
  ggplot2::geom_point(
    ggplot2::aes(color = treatment_status, shape = time_point_id),
    size = 3.5,
    stroke = 0.8
  ) +
  # Batch replicate as repelled text labels (third aesthetic dimension)
  ggrepel::geom_text_repel(
    ggplot2::aes(label = batch_id),
    size = 2.8,
    color = "grey35",
    max.overlaps = Inf,
    segment.size = 0.3
  ) +
  ggplot2::scale_color_manual(values = treatment_colours) +
  ggplot2::scale_fill_manual(values = treatment_colours, guide = "none") +
  ggplot2::scale_shape_manual(values = time_shapes) +
  ggplot2::labs(
    x = glue::glue("PC1 ({var_pct[1]}%)"),
    y = glue::glue("PC2 ({var_pct[2]}%)"),
    color = "Treatment",
    shape = "Time point",
    title = "PCA \u2014 Bulk RNA-seq samples"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "right")

pca_bulk_plot
ggplot2::ggsave(
  filename = file.path(output_dir, paste0(study, "_pca_bulk_plot_after_log_", today, ".pdf")),
  plot = pca_bulk_plot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)
