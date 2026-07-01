# ==========================================================================
# 0. Libraries and filename settings ----
# ==========================================================================

library(Seurat)
library(SummarizedExperiment)
library(tinytable)
library(ggplot2)
library(ggrepel)
library(mixOmics)

study <- "luque"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/bulk-omics"
source("./R/utils.R")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==========================================================================
# 1. Load pseudo-bulk data---
# ==========================================================================

naive_pseudo_bulk <- readRDS(
  file = "./data/intermediate/luque_naive_pseudo_bulk_2026-06-10.rds"
)


tinytable::tt(
  SummarizedExperiment::colData(naive_pseudo_bulk) |>
    as.data.frame(),
  caption = "Phenotype distribution in naive pseudo-bulk"
)
# Log-normalise raw counts; transpose so samples are rows (required by PCA)
log_counts <- SummarizedExperiment::assay(naive_pseudo_bulk, "counts") |>
  log1p() |>
  t()

col_meta <- SummarizedExperiment::colData(naive_pseudo_bulk) |>
  as.data.frame()

phenotype_colours <- stats::setNames(
  c("#3182bd", "#e6550d"),
  c("TLS", "neural_bias")
)


# ==========================================================================
# 2. PCA analysis on naive pseudo-bulk data----
# ==========================================================================

# Keep the top most variable genes (standard bulk RNA-seq PCA criterion).

N_HVG <- 500L
gene_vars <- apply(log_counts, 2, var)
hvg_idx <- order(gene_vars, decreasing = TRUE)[
  seq_len(min(N_HVG, sum(gene_vars > 0)))
]
log_counts_hvg <- log_counts[, hvg_idx]

# Compute PCA using the top 500 most variable genes
pca_res <- prcomp(log_counts_hvg, center = TRUE, scale. = TRUE)

# Add PC1/PC2 scores with sample metadata
pca_data <- as.data.frame(pca_res$x[, 1:2]) |>
  dplyr::mutate(
    phenotype = col_meta$Morphotype,
    barcode_id = col_meta$barcode_id
  )
var_pct <- round(pca_res$sdev[1:2]^2 / sum(pca_res$sdev^2) * 100, 1)

pca_bulk_plot <- ggplot2::ggplot(
  pca_data,
  ggplot2::aes(x = PC1, y = PC2, color = phenotype)
) +
  # Ellipses for treatment status only
  ggplot2::stat_ellipse(
    ggplot2::aes(color = phenotype, fill = phenotype),
    level = 0.95,
    geom = "polygon",
    alpha = 0.10,
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  # Points: colour = morphotype; barcode as repelled label (not shape:
  # 24 barcodes exceed ggplot2's default shape palette).
  ggplot2::geom_point(
    ggplot2::aes(color = phenotype),
    size = 3.5,
    stroke = 0.8
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = barcode_id),
    size = 2.8,
    color = "grey35",
    max.overlaps = Inf,
    segment.size = 0.3
  ) +
  ggplot2::scale_color_manual(values = phenotype_colours) +
  ggplot2::scale_fill_manual(values = phenotype_colours, guide = "none") +
  ggplot2::labs(
    x = glue::glue("PC1 ({var_pct[1]}%)"),
    y = glue::glue("PC2 ({var_pct[2]}%)"),
    color = "Morphotype",
    title = "PCA - naive pseudo-bulk samples"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "right")

pca_bulk_plot
ggplot2::ggsave(
  filename = file.path(output_dir, paste0(study, "_pca_naive_pseudo_bulk_", today, ".pdf")),
  plot = pca_bulk_plot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)


# ==========================================================================
# 3. Sparse PLS-DA on naive pseudo-bulk data ----
# ==========================================================================

# mixOmics::splsda: sparse partial least squares discriminant analysis
# (Le Cao et al., mixOmics vignette section 5.3). Y = Morphotype;
# keepX = 500 genes selected per latent dimension.

N_COMP <- 2L
N_KEEPX <- 500L
morphotype_y <- factor(col_meta$Morphotype)

splsda_fit <- mixOmics::splsda(
  X = log_counts,
  Y = morphotype_y,
  ncomp = N_COMP,
  keepX = rep(N_KEEPX, N_COMP)
)

# --- Sample scores: variates on components 1 and 2 ------------------------
splsda_scores <- as.data.frame(splsda_fit$variates$X[, 1:2, drop = FALSE])
colnames(splsda_scores) <- c("comp1", "comp2")
var_expl <- round(splsda_fit$prop_expl_var$X[1:2] * 100, 1)

splsda_data <- splsda_scores |>
  dplyr::mutate(
    phenotype = col_meta$Morphotype,
    barcode_id = col_meta$barcode_id
  )

splsda_sample_plot <- ggplot2::ggplot(
  splsda_data,
  ggplot2::aes(x = comp1, y = comp2, color = phenotype)
) +
  ggplot2::stat_ellipse(
    ggplot2::aes(color = phenotype, fill = phenotype),
    level = 0.95,
    geom = "polygon",
    alpha = 0.10,
    linewidth = 0.5,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = phenotype),
    size = 3.5,
    stroke = 0.8
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = barcode_id),
    size = 2.8,
    color = "grey35",
    max.overlaps = Inf,
    segment.size = 0.3
  ) +
  ggplot2::scale_color_manual(values = phenotype_colours) +
  ggplot2::scale_fill_manual(values = phenotype_colours, guide = "none") +
  ggplot2::labs(
    x = glue::glue("component 1 ({var_expl[1]}% expl. var.)"),
    y = glue::glue("component 2 ({var_expl[2]}% expl. var.)"),
    color = "Morphotype",
    title = "sPLS-DA - naive pseudo-bulk samples"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "right")

splsda_sample_plot
ggplot2::ggsave(
  filename = file.path(
    output_dir,
    paste0(study, "_splsda_samples_naive_pseudo_bulk_", today, ".pdf")
  ),
  plot = splsda_sample_plot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)

# --- Selected features: loadings from selectVar() per component -------------
extract_splsda_loadings <- function(fit, comp) {
  sel <- mixOmics::selectVar(fit, comp = comp)
  data.frame(
    gene = sel$name,
    loading = sel$value[[1L]],
    component = paste0("comp", comp),
    stringsAsFactors = FALSE
  )
}

splsda_loadings <- dplyr::bind_rows(
  extract_splsda_loadings(splsda_fit, comp = 1L),
  extract_splsda_loadings(splsda_fit, comp = 2L)
) |>
  dplyr::mutate(
    abs_loading = abs(.data$loading)
  )

# Top 25 genes per component by absolute loading (discriminant signature)
N_TOP_LOADINGS <- 25L
splsda_top_loadings <- splsda_loadings |>
  dplyr::group_by(.data$component) |>
  dplyr::slice_max(.data$abs_loading, n = N_TOP_LOADINGS, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    gene = forcats::fct_reorder(.data$gene, .data$abs_loading)
  )

splsda_loading_plot <- ggplot2::ggplot(
  splsda_top_loadings,
  ggplot2::aes(
    x = .data$loading,
    y = .data$gene,
    fill = .data$component
  )
) +
  ggplot2::geom_col() +
  ggplot2::facet_wrap(~ component, scales = "free_y") +
  ggplot2::labs(
    x = "sPLS-DA loading",
    y = NULL,
    fill = "Component",
    title = glue::glue(
      "Top {N_TOP_LOADINGS} selected genes per sPLS-DA component"
    ),
    subtitle = glue::glue(
      "Sparse selection: {N_KEEPX} genes per component (mixOmics::splsda)"
    )
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "none",
    strip.background = ggplot2::element_rect(fill = "grey92")
  )

splsda_loading_plot
ggplot2::ggsave(
  filename = file.path(
    output_dir,
    paste0(study, "_splsda_loadings_naive_pseudo_bulk_", today, ".pdf")
  ),
  plot = splsda_loading_plot,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300
)

tinytable::tt(
  splsda_top_loadings |>
    dplyr::select("component", "gene", "loading", "abs_loading") |>
    dplyr::arrange(.data$component, dplyr::desc(.data$abs_loading)),
  caption = glue::glue(
    "Top {N_TOP_LOADINGS} sPLS-DA selected genes per component ",
    "(keepX = {N_KEEPX})"
  )
)

