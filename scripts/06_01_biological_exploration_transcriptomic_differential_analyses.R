# ============================================================================
# 0. Libraries and filename settings ----
# ============================================================================

library(Seurat)
library(SummarizedExperiment)
library(DESeq2)
library(dplyr)
library(tibble)
library(purrr)
library(readr)
library(ggplot2)
library(EnhancedVolcano)
library(ComplexUpset)
library(patchwork)

study <- "GSE250136"
technique <- "pseudobulk_deseq2"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/biological-exploration/DEA-analyses"
source("./R/utils.R")
source("./R/simulate_pseudo_bulk_samples.R")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

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
# 3. DESeq2: cell-type-aware model (~ celltype * Morphotype) ----
# ============================================================================

col_data_celltype <- as.data.frame(SummarizedExperiment::colData(pb_celltype_barcode))
col_data_celltype$celltype <- factor(col_data_celltype$celltype)
col_data_celltype$Morphotype <- relevel(
  factor(col_data_celltype$Morphotype),
  ref = "TLS"
)

dds_celltype <- DESeqDataSetFromMatrix(
  countData = round(SummarizedExperiment::assay(pb_celltype_barcode, "counts")),
  colData = col_data_celltype,
  design = ~ celltype * Morphotype
)
dds_celltype <- DESeq(dds_celltype)

morphotype_coef <- "Morphotype_neural_bias_vs_TLS"
ref_celltype <- levels(col_data_celltype$celltype)[[1L]]

extract_celltype_contrast <- function(dds, celltype_level) {
  coef_names <- resultsNames(dds)
  contrast_vec <- numeric(length(coef_names))
  names(contrast_vec) <- coef_names
  contrast_vec[morphotype_coef] <- 1

  if (!identical(celltype_level, ref_celltype)) {
    interaction_coef <- grep(
      paste0("^celltype", make.names(celltype_level), "\\.Morphotype"),
      coef_names,
      value = TRUE
    )
    if (length(interaction_coef) == 0L) {
      stop("Missing interaction coefficient for cell type: ", celltype_level)
    }
    contrast_vec[interaction_coef[[1L]]] <- 1
  }

  results(dds, contrast = contrast_vec)
}

de_celltype_long <- bind_rows(lapply(levels(col_data_celltype$celltype), function(ct) {
  extract_celltype_contrast(dds_celltype, ct) |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    mutate(
      celltype = ct,
      padj_celltype_BH = p.adjust(pvalue, method = "BH")
    )
})) |>
  mutate(
    padj_global_BH = p.adjust(pvalue, method = "BH"),
    study = study,
    technique = technique,
    analysis_level = "celltype_barcode"
  ) |>
  relocate(study, technique, analysis_level, .before = gene)

fdr_celltype_wide <- de_celltype_long |>
  select(gene, celltype, padj_celltype_BH) |>
  tidyr::pivot_wider(
    names_from = celltype,
    values_from = padj_celltype_BH,
    names_prefix = "padj_"
  ) |>
  mutate(study = study, technique = technique, analysis_level = "celltype_barcode") |>
  relocate(study, technique, analysis_level, .before = gene)

# ============================================================================
# 4. DESeq2: sample-level global model (~ Morphotype) ----
# ============================================================================

col_data_sample <- as.data.frame(SummarizedExperiment::colData(pb_sample_barcode))
col_data_sample$Morphotype <- relevel(
  factor(col_data_sample$Morphotype),
  ref = "TLS"
)

dds_sample <- DESeqDataSetFromMatrix(
  countData = round(SummarizedExperiment::assay(pb_sample_barcode, "counts")),
  colData = col_data_sample,
  design = ~ Morphotype
)
dds_sample <- DESeq(dds_sample)

de_sample <- results(
  dds_sample,
  contrast = c("Morphotype", "neural_bias", "TLS"),
  alpha = 0.05
) |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  mutate(
    padj_BH = p.adjust(pvalue, method = "BH"),
    study = study,
    technique = technique,
    analysis_level = "sample_barcode"
  ) |>
  relocate(study, technique, analysis_level, .before = gene)

# ============================================================================
# 5. Save DE tables ----
# ============================================================================

readr::write_csv(
  de_celltype_long,
  file.path(
    output_dir,
    paste0(study, "_", technique, "_120h_celltype_model_results_", today, ".csv")
  )
)
readr::write_csv(
  fdr_celltype_wide,
  file.path(
    output_dir,
    paste0(study, "_", technique, "_120h_celltype_model_fdr_by_celltype_", today, ".csv")
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

volcano_df <- de_celltype_long |>
  mutate(neg_log10_p = -log10(pvalue))

x_limit <- max(abs(volcano_df$log2FoldChange), na.rm = TRUE)
x_limit <- max(0.6, x_limit)
y_limit <- max(volcano_df$neg_log10_p, na.rm = TRUE)
y_limit <- max(2, y_limit)

build_volcano <- function(df_ct) {
  EnhancedVolcano::EnhancedVolcano(
    df_ct,
    lab = df_ct$gene,
    x = "log2FoldChange",
    y = "pvalue",
    title = unique(df_ct$celltype),
    subtitle = paste0(study, " | ", technique, " | cell-type model"),
    xlab = bquote(Log[2] ~ "fold change (neural_bias vs TLS)"),
    ylab = bquote(-Log[10] ~ "p-value"),
    pCutoff = 0.05,
    FCcutoff = 0.5,
    xlim = c(-x_limit, x_limit),
    ylim = c(0, y_limit),
    pointSize = 1.5,
    labSize = 2.5,
    drawConnectors = FALSE,
    legendPosition = "right"
  )
}

volcano_plots <- split(de_celltype_long, de_celltype_long$celltype) |>
  purrr::map(build_volcano)

volcano_file <- file.path(
  output_dir,
  paste0(study, "_", technique, "_120h_celltype_model_volcano_", today, ".pdf")
)

gr_devices <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
gr_devices(file = volcano_file, width = 24, height = 14)
print(
  patchwork::wrap_plots(volcano_plots, ncol = 4) +
    patchwork::plot_annotation(
      title = paste0(study, " | ", technique, " | 120h cell-type DESeq2 model")
    )
)
invisible(grDevices::dev.off())

# ============================================================================
# 7. ComplexUpset plot (cell-type model DEG overlaps) ----
# ============================================================================

deg_hits <- de_celltype_long |>
  filter(!is.na(pvalue), !is.na(log2FoldChange)) |>
  filter(pvalue < 0.05, abs(log2FoldChange) > 0.5) |>
  select(gene, celltype) |>
  distinct()

upset_file <- file.path(
  output_dir,
  paste0(study, "_", technique, "_120h_celltype_model_deg_overlap_", today, ".pdf")
)

if (nrow(deg_hits) > 0) {
  membership <- deg_hits |>
    mutate(present = TRUE) |>
    tidyr::pivot_wider(
      names_from = celltype,
      values_from = present,
      values_fill = FALSE
    )

  membership_cols <- setdiff(colnames(membership), "gene")

  gr_devices(file = upset_file, width = 18, height = 10)
  print(
    ComplexUpset::upset(
      membership,
      intersect = membership_cols,
      min_size = 1,
      n_intersections = 20,
      name = "Shared DEGs"
    ) +
      ggplot2::ggtitle(
        paste0(study, " | ", technique, " | DEG overlap (p < 0.05, |LFC| > 0.5)")
      )
  )
  invisible(grDevices::dev.off())
} else {
  message(
    "No DEGs passed pvalue < 0.05 and |log2FoldChange| > 0.5; ",
    "skipping ComplexUpset plot."
  )
}

message("Pseudo-bulk DESeq2 workflow completed.")
