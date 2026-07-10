# ============================================================================
# Compare DEGs from MAST vs pseudobulk DESeq2 (Venn diagrams)
# ============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(ggVennDiagram)
library(patchwork)

study <- "GSE250136"
output_dir <- "outputs/biological-exploration/DEA-analyses"
results_date <- "2026-07-09"

pvalue_cutoff <- 0.05
lfc_cutoff <- 0.5

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# 1. Load DE result tables ----
# ============================================================================

degs_mast_celltype <- readr::read_csv(
  file.path(
    output_dir,
    paste0(
      study,
      "_mast_120h_celltype_model_results_",
      results_date,
      ".csv"
    )
  ),
  show_col_types = FALSE
)
degs_deseq2_celltype <- readr::read_csv(
  file.path(
    output_dir,
    paste0(
      study,
      "_pseudobulk_deseq2_120h_celltype_model_results_",
      results_date,
      ".csv"
    )
  ),
  show_col_types = FALSE
)
degs_mast_sample <- readr::read_csv(
  file.path(
    output_dir,
    paste0(
      study,
      "_mast_120h_sample_model_results_",
      results_date,
      ".csv"
    )
  ),
  show_col_types = FALSE
)
degs_deseq2_sample <- readr::read_csv(
  file.path(
    output_dir,
    paste0(
      study,
      "_pseudobulk_deseq2_120h_sample_model_results_",
      results_date,
      ".csv"
    )
  ),
  show_col_types = FALSE
)

# ============================================================================
# 2. Helpers ----
# ============================================================================

filter_degs <- function(df_de) {
  df_de |>
    dplyr::filter(
      !is.na(pvalue),
      !is.na(log2FoldChange),
      pvalue < pvalue_cutoff,
      abs(log2FoldChange) > lfc_cutoff
    ) |>
    dplyr::distinct(gene) |>
    dplyr::pull(gene)
}

build_venn_plot <- function(gene_sets, plot_title) {
  ggVennDiagram::ggVennDiagram(
    gene_sets,
    set_size = 4,
    label = "count"
  ) +
    ggplot2::scale_fill_gradient(low = "grey90", high = "red") +
    ggplot2::ggtitle(plot_title) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.22, 0.08))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
      legend.position = "none"
    )
}

pretty_celltype <- function(celltype) {
  gsub("\\.", " ", celltype, fixed = FALSE) |>
    tools::toTitleCase()
}

# ============================================================================
# 3. Per-cell-type Venn diagrams ----
# ============================================================================

celltypes <- sort(unique(c(
  degs_mast_celltype$celltype,
  degs_deseq2_celltype$celltype
)))

venn_plots_celltype <- lapply(celltypes, function(ct) {
  mast_genes <- degs_mast_celltype |>
    dplyr::filter(celltype == ct) |>
    filter_degs()

  deseq2_genes <- degs_deseq2_celltype |>
    dplyr::filter(celltype == ct) |>
    filter_degs()

  gene_sets <- list(
    MAST = mast_genes,
    "DESeq2 pseudobulk" = deseq2_genes
  )

  build_venn_plot(
    gene_sets,
    plot_title = paste0(
      pretty_celltype(ct),
      " (n = ",
      length(mast_genes),
      " vs ",
      length(deseq2_genes),
      ")"
    )
  )
})
names(venn_plots_celltype) <- celltypes

# ============================================================================
# 4. Sample-level (global) Venn diagram ----
# ============================================================================

mast_sample_genes <- filter_degs(degs_mast_sample)
deseq2_sample_genes <- filter_degs(degs_deseq2_sample)

venn_plot_sample <- build_venn_plot(
  list(
    MAST = mast_sample_genes,
    "DESeq2 pseudobulk" = deseq2_sample_genes
  ),
  plot_title = paste0(
    "Sample level (n = ",
    length(mast_sample_genes),
    " vs ",
    length(deseq2_sample_genes),
    ")"
  )
)

# ============================================================================
# 5. Combine and save single PDF ----
# ============================================================================

venn_figure <- patchwork::wrap_plots(
  c(venn_plots_celltype, list(Sample = venn_plot_sample)),
  ncol = 2
) +
  patchwork::plot_annotation(
    title = paste0(
      study,
      " | MAST vs DESeq2 DEG overlap (p < ",
      pvalue_cutoff,
      ", |log2FC| > ",
      lfc_cutoff,
      ")"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
    )
  )

venn_file <- file.path(
  output_dir,
  paste0(study, "_mast_vs_deseq2_deg_venn_120h.pdf")
)

ggplot2::ggsave(
  filename = venn_file,
  plot = venn_figure,
  width = 14,
  height = 18,
  units = "in",
  dpi = 300
)

message("Saved: ", venn_file)
