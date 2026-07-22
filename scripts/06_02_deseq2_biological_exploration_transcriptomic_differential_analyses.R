# ============================================================================
# 0. Filename settings ----
# ============================================================================

study <- "GSE250136"
technique <- "pseudobulk_deseq2"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/biological-exploration/DEA-analyses"
source("./R/utils.R")
source("./R/simulate_pseudo_bulk_samples.R")
source("./R/filtering_unexpressed_genes.R")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# see also https://rdrr.io/bioc/muscat/man/pbDS.html, describing distinct pseudo-bulk strategies
# at sample, or cell-type elvel, using traditional bulk approaches

# ============================================================================
# 1. Load merged object and subset to 120h ----
# ============================================================================

# Methodological note:
# DESeq2 requires unnormalised raw counts (Love et al., DESeq2 vignette).
# Pseudo-bulking by sample and cell type avoids cell-level pseudoreplication
# (Seurat DE vignette; Nguyen et al., Nat Commun 2023).
obj <- readRDS("./data/intermediate/luque_single_cell_merged_2026-06-07.rds")
SeuratObject::DefaultAssay(obj) <- "RNA"

cells_120h <- obj[[]]$timepoint == "120h" &
  obj[[]]$Morphotype %in% c("TLS", "neural_bias") &
  !is.na(obj[[]]$luque_cluster_annotation)

# ============================================================================
# 2. Pseudo-bulk aggregation ----
# ============================================================================

pb_celltype_barcode <- aggregate_celltype_barcode_pseudo_bulk(
  seurat_obj = obj,
  celltype_col = "luque_cluster_annotation",
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  cell_mask = cells_120h
)

pb_sample_barcode <- aggregate_barcode_pseudo_bulk(
  seurat_obj = obj,
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  cell_mask = cells_120h
)

# ============================================================================
# 3. DESeq2: per-cell-type pseudo-bulk models (~ Morphotype) ----
# ============================================================================

# Filter genes separately within each cell type after pseudo-bulk aggregation.
# edgeR::filterByExpr() applies a design-aware CPM/count filter, preserving
# genes that are sufficiently expressed in biological replicates of at least
# one Morphotype instead of using a single global expression universe.
min_pseudobulk_cells <- 10L
deseq2_keep_genes_by_celltype <- filter_deseq2_unexpressed_genes(
  seurat_obj = obj,
  celltype_col = "luque_cluster_annotation",
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  assay = "RNA",
  cell_mask = cells_120h,
  min_cells_per_sample = min_pseudobulk_cells,
  min_count = 5L,
  min_total_count = 15L
)
names(deseq2_keep_genes_by_celltype) <- make.names(names(
  deseq2_keep_genes_by_celltype
))

col_data_celltype <- as.data.frame(SummarizedExperiment::colData(
  pb_celltype_barcode
))
col_data_celltype$celltype <- make.names(col_data_celltype$celltype)
count_data_celltype <- round(SummarizedExperiment::assay(
  pb_celltype_barcode,
  "counts"
))

selected_genes_celltype <- tibble::tibble(
  celltype = names(deseq2_keep_genes_by_celltype),
  n_selected_genes = lengths(deseq2_keep_genes_by_celltype)
)

# Run differential expression independently for each cell type. Each fit tests
# neural_bias versus TLS within that cell type, using only the genes retained
# by the cell-type-specific filter above.
run_deseq2_celltype <- function(ct) {
  keep_samples <- col_data_celltype$celltype == ct &
    col_data_celltype$n_cells >= min_pseudobulk_cells
  col_data_ct <- col_data_celltype[keep_samples, , drop = FALSE]
  genes_ct <- deseq2_keep_genes_by_celltype[[ct]]

  if (nrow(col_data_ct) == 0L || length(genes_ct) == 0L) {
    return(NULL)
  }

  col_data_ct$Morphotype <- stats::relevel(
    factor(col_data_ct$Morphotype),
    ref = "TLS"
  )
  if (length(unique(col_data_ct$Morphotype)) < 2L) {
    return(NULL)
  }

  count_data_ct <- count_data_celltype[
    genes_ct,
    rownames(col_data_ct),
    drop = FALSE
  ]
  dds_ct <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data_ct,
    colData = col_data_ct,
    design = ~Morphotype
  )
  dds_ct <- DESeq2::DESeq(
    dds_ct,
    test = "LRT",
    reduced = ~1,
    minmu = 1e-6,
    minReplicatesForReplace = Inf
  )

  DESeq2::results(
    dds_ct,
    contrast = c("Morphotype", "neural_bias", "TLS"),
    alpha = 0.05,
    independentFiltering = TRUE
  ) |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(
      celltype = ct,
      n_selected_genes = length(genes_ct),
      padj_celltype_BH = stats::p.adjust(pvalue, method = "BH")
    )
}

celltypes <- sort(unique(col_data_celltype$celltype))
de_celltype_long <- dplyr::bind_rows(lapply(
  celltypes,
  run_deseq2_celltype
)) |>
  dplyr::mutate(
    padj_global_BH = stats::p.adjust(pvalue, method = "BH"),
    study = study,
    technique = technique,
    analysis_level = "celltype_barcode"
  ) |>
  dplyr::relocate(study, technique, analysis_level, .before = gene)

de_celltype_wide <- de_celltype_long |>
  tidyr::pivot_wider(
    id_cols = c(study, technique, analysis_level, gene),
    names_from = celltype,
    values_from = c(log2FoldChange, pvalue, padj_celltype_BH),
    names_glue = "{celltype}__{.value}"
  )

# ============================================================================
# 4. DESeq2: sample-level global model (~ Morphotype) ----
# ============================================================================

col_data_sample <- as.data.frame(SummarizedExperiment::colData(
  pb_sample_barcode
))
col_data_sample$Morphotype <- stats::relevel(
  factor(col_data_sample$Morphotype),
  ref = "TLS"
)
count_data_sample <- round(SummarizedExperiment::assay(
  pb_sample_barcode,
  "counts"
))
sample_design <- stats::model.matrix(~Morphotype, data = col_data_sample)
keep_genes_sample <- edgeR::filterByExpr(
  y = count_data_sample,
  design = sample_design,
  min.count = 5L,
  min.total.count = 15L
)
count_data_sample <- count_data_sample[keep_genes_sample, , drop = FALSE]

dds_sample <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_data_sample,
  colData = col_data_sample,
  design = ~Morphotype
)
dds_sample <- DESeq2::DESeq(
  dds_sample,
  test = "LRT",
  reduced = ~1,
  minmu = 1e-6,
  minReplicatesForReplace = Inf
)

de_sample <- DESeq2::results(
  dds_sample,
  contrast = c("Morphotype", "neural_bias", "TLS"),
  alpha = 0.05
) |>
  as.data.frame() |>
  tibble::rownames_to_column("gene") |>
  dplyr::mutate(
    n_selected_genes = nrow(count_data_sample),
    padj_BH = stats::p.adjust(pvalue, method = "BH"),
    study = study,
    technique = technique,
    analysis_level = "sample_barcode"
  ) |>
  dplyr::relocate(study, technique, analysis_level, .before = gene)

# ============================================================================
# 5. Save DE tables ----
# ============================================================================

readr::write_csv(
  de_celltype_long,
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
# 6. Volcano plots (cell-type model; raw p-value, |log2FC| > 0.5) ----
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
  n_selected <- unique(df_de$n_selected_genes)
  n_selected <- n_selected[!is.na(n_selected)][[1L]]

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

  p_volcano <- EnhancedVolcano::EnhancedVolcano(
    df_de,
    lab = df_de$gene,
    x = "log2FoldChange",
    y = "pvalue",
    title = plot_title,
    subtitle = paste0(study, " | ", technique, " | tested genes: ", n_selected),
    xlab = bquote(Log[2] ~ "fold change (neural_bias vs TLS)"),
    ylab = bquote(-Log[10] ~ "p-value"),
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

  p_volcano
}

celltypes <- sort(unique(de_celltype_long$celltype))
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
# 7. P-value distribution histograms ----
# ============================================================================

build_pvalue_histogram <- function(
  p_values,
  plot_title,
  p_label,
  fill_colour,
  n_selected_genes
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
      subtitle = paste0(p_label, " | tested genes: ", n_selected_genes),
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
  n_selected <- unique(df_ct$n_selected_genes)
  n_selected <- n_selected[!is.na(n_selected)][[1L]]
  cowplot::plot_grid(
    build_pvalue_histogram(
      p_values = df_ct$pvalue,
      plot_title = ct,
      p_label = "raw p-value",
      fill_colour = "#4C78A8",
      n_selected_genes = n_selected
    ),
    build_pvalue_histogram(
      p_values = df_ct$padj_celltype_BH,
      plot_title = ct,
      p_label = "BH adjusted p-value",
      fill_colour = "#F58518",
      n_selected_genes = n_selected
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
    p_label = "raw p-value",
    fill_colour = "#4C78A8",
    n_selected_genes = unique(de_sample$n_selected_genes)[[1L]]
  ),
  build_pvalue_histogram(
    p_values = de_sample$padj_BH,
    plot_title = "sample-level global",
    p_label = "BH adjusted p-value",
    fill_colour = "#F58518",
    n_selected_genes = unique(de_sample$n_selected_genes)[[1L]]
  ),
  ncol = 2L,
  labels = c("raw", "adjusted")
)

hist_pages <- c(
  hist_page_list,
  list(hist_page_global)
)
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
# 8. ComplexUpset plot (cell-type model DEG overlaps) ----
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

message("Pseudo-bulk DESeq2 workflow completed.")
