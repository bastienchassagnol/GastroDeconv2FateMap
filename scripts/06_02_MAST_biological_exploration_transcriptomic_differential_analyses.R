# ============================================================================
# 0. Filename settings ----
# ============================================================================

# bash command to run this script ----
# nohup Rscript --no-save --no-restore \
#   scripts/06_02_MAST_biological_exploration_transcriptomic_differential_analyses.R \
#   > "logs/GSE250136_$(date +%F)_mast_hurdle_glmer.log" 2>&1 &

library(MAST)
library(SeuratObject)
library(Matrix)
library(S4Vectors)
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(ggplot2)
library(EnhancedVolcano)
library(ComplexUpset)
library(cowplot)
library(gridExtra)
library(grid)

study <- "GSE250136"
technique <- "mast"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/biological-exploration/DEA-analyses"
source("./R/utils.R")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# set the number of cores to use
options(mc.cores = parallel::detectCores() %/% 2L)

# ============================================================================
# 1. Load merged object, subset to 120h, extract normalized expression ----
# ============================================================================

# Methodological note:
# MAST models log-normalised single-cell expression with a hurdle GLM
# (Finak et al., Genome Biology 2015). We use RNA/data (not raw counts),
# adjust for cellular detection rate (cngeneson), and include organoid-level
# random effects to avoid pseudo-replication.
obj <- readRDS("./data/intermediate/luque_single_cell_merged_2026-06-07.rds")
SeuratObject::DefaultAssay(obj) <- "RNA"

cells_120h <- obj[[]]$timepoint == "120h" &
  obj[[]]$Morphotype %in% c("TLS", "neural_bias") &
  !is.na(obj[[]]$luque_cluster_annotation)

meta_120h <- obj[[]][cells_120h, , drop = FALSE]
expr_data <- SeuratObject::LayerData(
  object = obj,
  assay = "RNA",
  layer = "data"
)[, cells_120h, drop = FALSE]

meta_120h <- meta_120h[
  match(colnames(expr_data), meta_120h$cell_id),
  ,
  drop = FALSE
]
rownames(meta_120h) <- colnames(expr_data)

meta_120h$celltype <- make.names(meta_120h$luque_cluster_annotation)
meta_120h$Morphotype <- stats::relevel(
  factor(meta_120h$Morphotype),
  ref = "TLS"
)
meta_120h$Sample.barcode <- factor(meta_120h$Sample.barcode)
meta_120h$wellKey <- rownames(meta_120h)

# ============================================================================
# 2. Gene filtering and MAST SingleCellAssay ----
# ============================================================================

min_detection_fraction <- 0.80
keep_genes <- Matrix::rowMeans(expr_data > 0) >= min_detection_fraction
expr_data <- expr_data[keep_genes, , drop = FALSE]
expr_mat <- as.matrix(expr_data)

message(
  "Retained ",
  nrow(expr_mat),
  " genes detected in >= ",
  100 * min_detection_fraction,
  "% of 120h cells (n = ",
  ncol(expr_mat),
  ")."
)

meta_120h$cngeneson <- scale(Matrix::colSums(expr_mat > 0))[, 1]

fdata <- S4Vectors::DataFrame(
  primerid = rownames(expr_mat),
  row.names = rownames(expr_mat)
)
sca <- MAST::FromMatrix(
  exprsArray = expr_mat,
  cData = meta_120h,
  fData = fdata,
  check_sanity = FALSE
)

# ============================================================================
# 3. MAST helpers ----
# ============================================================================

mast_pvalue_col <- function(df) {
  hits <- grep("Pr.*Chisq", colnames(df), value = TRUE)
  if (length(hits) == 0L) {
    stop("Could not find MAST hurdle p-value column.")
  }
  hits[[1L]]
}

compute_detection_freq <- function(celltype_level) {
  idx_ct <- meta_120h$celltype == celltype_level
  idx_tls <- idx_ct & meta_120h$Morphotype == "TLS"
  idx_nb <- idx_ct & meta_120h$Morphotype == "neural_bias"
  tibble::tibble(
    gene = rownames(expr_mat),
    pct_TLS = Matrix::rowMeans(expr_mat[, idx_tls, drop = FALSE] > 0),
    pct_neural_bias = Matrix::rowMeans(expr_mat[, idx_nb, drop = FALSE] > 0),
    detection_freq = Matrix::rowMeans(expr_mat[, idx_ct, drop = FALSE] > 0)
  )
}

morphotype_lrt_name <- function(celltype_level) {
  paste0("celltype", celltype_level, ":Morphotypeneural_bias")
}

morphotype_lfc_contrast <- function(celltype_level) {
  paste0("celltype", celltype_level, ":Morphotypeneural_bias")
}

extract_mast_hurdle <- function(zlm_fit, summary_fit, lfc_contrast) {
  res_dt <- summary_fit$datatable
  res_h <- res_dt[res_dt$component == "H", , drop = FALSE]
  p_col <- mast_pvalue_col(res_h)

  lfc_dt <- MAST::getLogFC(zlm_fit)
  lfc_h <- lfc_dt[
    as.character(lfc_dt$contrast) == lfc_contrast,
    c("primerid", "logFC"),
    drop = FALSE
  ]

  res_h |>
    as.data.frame() |>
    dplyr::rename(
      gene = primerid,
      pvalue = dplyr::all_of(p_col)
    ) |>
    dplyr::left_join(
      as.data.frame(lfc_h),
      by = c("gene" = "primerid")
    ) |>
    dplyr::mutate(log2FoldChange = logFC)
}

celltypes <- sort(unique(meta_120h$celltype))

# ============================================================================
# 4. MAST: cell-type interaction model ----
# ============================================================================

message("Fitting MAST hurdle mixed model for cell-type interaction...")
zlm_celltype <- MAST::zlm(
  formula = ~ 0 +
    celltype +
    celltype:Morphotype +
    cngeneson +
    (1 | Sample.barcode),
  sca = sca,
  method = "glmer",
  ebayes = FALSE,
  parallel = TRUE,
  silent = FALSE
)

extract_mast_celltype <- function(celltype_level) {
  message("Extracting results for cell type: ", celltype_level)
  lrt_term <- morphotype_lrt_name(celltype_level)
  lfc_term <- morphotype_lfc_contrast(celltype_level)
  summary_fit <- summary(zlm_celltype, doLRT = lrt_term)

  extract_mast_hurdle(zlm_celltype, summary_fit, lfc_term) |>
    dplyr::mutate(
      celltype = celltype_level,
      padj_celltype_BH = stats::p.adjust(pvalue, method = "BH")
    ) |>
    dplyr::left_join(compute_detection_freq(celltype_level), by = "gene")
}

de_celltype_long <- dplyr::bind_rows(lapply(
  celltypes,
  extract_mast_celltype
)) |>
  dplyr::mutate(
    padj_global_BH = stats::p.adjust(pvalue, method = "BH"),
    study = study,
    technique = technique,
    analysis_level = "celltype_single_cell"
  ) |>
  dplyr::relocate(study, technique, analysis_level, .before = gene)

de_celltype_wide <- de_celltype_long |>
  tidyr::pivot_wider(
    id_cols = c(study, technique, analysis_level, gene),
    names_from = celltype,
    values_from = c(
      log2FoldChange,
      pvalue,
      padj_celltype_BH,
      pct_TLS,
      pct_neural_bias,
      detection_freq
    ),
    names_glue = "{celltype}__{.value}"
  )

# ============================================================================
# 5. MAST: sample-level global model (~ Morphotype) ----
# ============================================================================

message("Fitting MAST hurdle mixed model for sample-level global contrast...")
zlm_sample <- MAST::zlm(
  formula = ~ Morphotype + cngeneson + (1 | Sample.barcode),
  sca = sca,
  method = "glmer",
  ebayes = FALSE,
  parallel = TRUE
)

summary_sample <- summary(zlm_sample, doLRT = "Morphotypeneural_bias")

de_sample <- extract_mast_hurdle(
  zlm_sample,
  summary_sample,
  "Morphotypeneural_bias"
) |>
  dplyr::mutate(
    padj_BH = stats::p.adjust(pvalue, method = "BH"),
    pct_TLS = Matrix::rowMeans(
      expr_mat[, meta_120h$Morphotype == "TLS", drop = FALSE] > 0
    ),
    pct_neural_bias = Matrix::rowMeans(
      expr_mat[, meta_120h$Morphotype == "neural_bias", drop = FALSE] > 0
    ),
    detection_freq = Matrix::rowMeans(expr_mat > 0),
    study = study,
    technique = technique,
    analysis_level = "sample_single_cell"
  ) |>
  dplyr::relocate(study, technique, analysis_level, .before = gene)

# ============================================================================
# 6. Save DE tables ----
# ============================================================================

readr::write_csv(
  de_celltype_wide,
  file.path(
    output_dir,
    paste0(
      study,
      "_",
      technique,
      "_120h_celltype_model_results_",
      today,
      ".csv"
    )
  )
)
readr::write_csv(
  de_sample,
  file.path(
    output_dir,
    paste0(study, "_", technique, "_120h_sample_model_results_", today, ".csv")
  )
)

# ============================================================================
# 7. Volcano plots ----
# ============================================================================

p_raw_cutoff <- 0.05

build_volcano <- function(
  df_de,
  plot_title,
  padj_col,
  show_legend = FALSE
) {
  df_de <- df_de |>
    dplyr::mutate(neg_log10_p = -log10(pvalue))

  x_limit <- max(abs(df_de$log2FoldChange), na.rm = TRUE)
  x_limit <- max(0.6, x_limit)

  y_limit <- max(df_de$neg_log10_p, na.rm = TRUE)
  if (!is.finite(y_limit)) {
    y_limit <- 2
  }
  y_limit <- max(2, y_limit)

  raw_line <- -log10(p_raw_cutoff)
  padj_sig <- df_de[[padj_col]] < p_raw_cutoff & !is.na(df_de[[padj_col]])
  if (any(padj_sig, na.rm = TRUE)) {
    padj_line <- -log10(max(df_de$pvalue[padj_sig], na.rm = TRUE))
  } else {
    padj_line <- raw_line
  }

  raw_label_y <- raw_line + 0.05 * y_limit
  padj_label_y <- padj_line + 0.15 * y_limit
  if (
    is.finite(raw_label_y) &&
      is.finite(padj_label_y) &&
      abs(padj_label_y - raw_label_y) < 0.08 * y_limit
  ) {
    padj_label_y <- raw_label_y + 0.12 * y_limit
  }

  EnhancedVolcano::EnhancedVolcano(
    df_de,
    lab = df_de$gene,
    x = "log2FoldChange",
    y = "pvalue",
    title = plot_title,
    subtitle = paste0(study, " | ", technique),
    xlab = bquote(Log[2] ~ "fold change (neural_bias vs TLS)"),
    ylab = bquote(-Log[10] ~ "hurdle p-value"),
    pCutoff = p_raw_cutoff,
    FCcutoff = 0.5,
    xlim = c(-x_limit, x_limit),
    ylim = c(0, y_limit),
    pointSize = 1.5,
    labSize = 2.5,
    drawConnectors = FALSE,
    legendPosition = if (show_legend) "bottom" else "none"
  ) +
    ggplot2::geom_hline(
      yintercept = raw_line,
      colour = "red",
      linewidth = 0.4,
      linetype = "solid"
    ) +
    ggplot2::geom_hline(
      yintercept = padj_line,
      colour = "blue",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    ggplot2::annotate(
      "text",
      x = -x_limit * 0.95,
      y = raw_label_y,
      label = "raw p = 0.05",
      colour = "red",
      hjust = 0,
      size = 3
    ) +
    ggplot2::annotate(
      "text",
      x = -x_limit * 0.95,
      y = padj_label_y,
      label = "BH padj = 0.05",
      colour = "blue",
      hjust = 0,
      size = 3
    )
}

volcano_plot_list <- lapply(celltypes, function(ct) {
  df_ct <- de_celltype_long[de_celltype_long$celltype == ct, , drop = FALSE]
  build_volcano(
    df_de = df_ct,
    plot_title = ct,
    padj_col = "padj_celltype_BH",
    show_legend = FALSE
  )
})
names(volcano_plot_list) <- celltypes

volcano_sample <- build_volcano(
  df_de = de_sample,
  plot_title = "sample-level global",
  padj_col = "padj_BH",
  show_legend = FALSE
)

volcano_legend_plot <- build_volcano(
  df_de = de_celltype_long[
    de_celltype_long$celltype == celltypes[[1L]],
    ,
    drop = FALSE
  ],
  plot_title = celltypes[[1L]],
  padj_col = "padj_celltype_BH",
  show_legend = TRUE
)
volcano_legend <- cowplot::get_legend(volcano_legend_plot)

n_celltypes <- length(volcano_plot_list)
n_row <- ceiling((n_celltypes + 1L) / 2L)
n_slot <- n_row * 2L
lower_left_slot <- (n_row - 1L) * 2L + 1L
volcano_slots <- vector("list", n_slot)
volcano_slots[[lower_left_slot]] <- volcano_sample
ct_idx <- 1L
for (slot_idx in seq_len(n_slot)) {
  if (is.null(volcano_slots[[slot_idx]])) {
    volcano_slots[[slot_idx]] <- volcano_plot_list[[ct_idx]]
    ct_idx <- ct_idx + 1L
  }
}

volcano_grob_list <- lapply(volcano_slots, ggplot2::ggplotGrob)
volcano_panel_grid <- do.call(
  gridExtra::arrangeGrob,
  c(volcano_grob_list, list(ncol = 2L))
)
volcano_figure <- gridExtra::arrangeGrob(
  volcano_panel_grid,
  volcano_legend,
  ncol = 1L,
  heights = grid::unit.c(grid::unit(1, "null"), grid::unit(1.2, "cm"))
)

volcano_file <- file.path(
  output_dir,
  paste0(study, "_", technique, "_120h_celltype_model_volcano_", today, ".pdf")
)
ggplot2::ggsave(
  filename = volcano_file,
  plot = volcano_figure,
  width = 14,
  height = 3.5 * n_row + 1,
  units = "in",
  dpi = 500
)

# ============================================================================
# 8. P-value distribution histograms ----
# ============================================================================

build_pvalue_histogram <- function(
  p_values,
  plot_title,
  p_label,
  fill_colour
) {
  df_hist <- tibble::tibble(p_value = p_values) |>
    dplyr::filter(!is.na(p_value))

  ggplot2::ggplot(df_hist, ggplot2::aes(x = p_value)) +
    ggplot2::geom_histogram(
      bins = 50,
      fill = fill_colour,
      colour = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_vline(
      xintercept = p_raw_cutoff,
      colour = "red",
      linewidth = 0.4
    ) +
    ggplot2::labs(
      title = plot_title,
      subtitle = p_label,
      x = "p-value",
      y = "Gene count"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 9),
      plot.subtitle = ggplot2::element_text(size = 8)
    )
}

build_celltype_hist_row <- function(ct) {
  df_ct <- de_celltype_long[de_celltype_long$celltype == ct, , drop = FALSE]
  cowplot::plot_grid(
    build_pvalue_histogram(
      p_values = df_ct$pvalue,
      plot_title = ct,
      p_label = "raw hurdle p-value",
      fill_colour = "#4C78A8"
    ),
    build_pvalue_histogram(
      p_values = df_ct$padj_celltype_BH,
      plot_title = ct,
      p_label = "BH adjusted p-value",
      fill_colour = "#F58518"
    ),
    ncol = 2L,
    labels = c("raw", "adjusted")
  )
}

hist_page_list <- list()
celltype_pairs <- split(
  celltypes,
  ceiling(seq_along(celltypes) / 2L)
)
for (pair_idx in seq_along(celltype_pairs)) {
  pair_rows <- lapply(celltype_pairs[[pair_idx]], build_celltype_hist_row)
  hist_page_list[[pair_idx]] <- cowplot::plot_grid(
    plotlist = pair_rows,
    ncol = 1L
  )
}

hist_page_global <- cowplot::plot_grid(
  build_pvalue_histogram(
    p_values = de_sample$pvalue,
    plot_title = "sample-level global",
    p_label = "raw hurdle p-value",
    fill_colour = "#4C78A8"
  ),
  build_pvalue_histogram(
    p_values = de_sample$padj_BH,
    plot_title = "sample-level global",
    p_label = "BH adjusted p-value",
    fill_colour = "#F58518"
  ),
  ncol = 2L,
  labels = c("raw", "adjusted")
)

hist_pages <- c(hist_page_list, list(hist_page_global))
hist_figure <- gridExtra::marrangeGrob(
  grobs = hist_pages,
  nrow = 1,
  ncol = 1,
  top = paste0(study, " | ", technique, " | p-value distributions")
)

pvalue_hist_file <- file.path(
  output_dir,
  paste0(study, "_", technique, "_120h_pvalue_histograms_", today, ".pdf")
)
ggplot2::ggsave(
  filename = pvalue_hist_file,
  plot = hist_figure,
  width = 12,
  height = 7,
  units = "in",
  dpi = 500
)

# ============================================================================
# 9. ComplexUpset plot ----
# ============================================================================

deg_hits <- de_celltype_long |>
  dplyr::filter(!is.na(pvalue), !is.na(log2FoldChange)) |>
  dplyr::filter(pvalue < 0.05, abs(log2FoldChange) > 0.5) |>
  dplyr::select(gene, celltype) |>
  dplyr::distinct()

upset_file <- file.path(
  output_dir,
  paste0(
    study,
    "_",
    technique,
    "_120h_celltype_model_deg_overlap_",
    today,
    ".pdf"
  )
)

if (nrow(deg_hits) > 0) {
  membership <- deg_hits |>
    dplyr::mutate(present = TRUE) |>
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = present,
      values_fill = FALSE
    )

  membership_cols <- setdiff(colnames(membership), "gene")

  upset_plot <- ComplexUpset::upset(
    membership,
    intersect = membership_cols,
    min_size = 1,
    n_intersections = 20,
    name = "Shared DEGs"
  ) +
    ggplot2::ggtitle(
      paste0(study, " | ", technique, " | DEG overlap (p < 0.05, |LFC| > 0.5)")
    )

  ggplot2::ggsave(
    filename = upset_file,
    plot = upset_plot,
    width = 18,
    height = 10,
    units = "in",
    dpi = 300
  )
} else {
  message(
    "No DEGs passed pvalue < 0.05 and |log2FoldChange| > 0.5; ",
    "skipping ComplexUpset plot."
  )
}

message("MAST hurdle mixed-model workflow completed.")
