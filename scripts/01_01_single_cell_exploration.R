# ==========================================================================
# 0. Libraries and Filename settings ----
# ==========================================================================

library(Seurat)
library(tinytable)
library(patchwork)
library(ggplot2)
library(gridExtra)

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/single-cell"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

source("./R/utils.R")


# ==========================================================================
# 1. Load the Seurat object ----
# ==========================================================================

GSE229513_gastruloids_seurat <- readRDS(
  "./data/raw/GSE229513_gastruloidsobject.rds"
)

dim(GSE229513_gastruloids_seurat)
table(Seurat::Idents(GSE229513_gastruloids_seurat))


raw_counts <- SeuratObject::GetAssayData(
  object = GSE229513_gastruloids_seurat,
  assay = "RNA",
  layer = "counts"
)

summarise_seurat_assays_layers(GSE229513_gastruloids_seurat)


# Quick inspection of counts matrix (first genes x first cells)

# Retrieve feature (gene) information from RNA assay
feature_info <- GSE229513_gastruloids_seurat@assays$RNA@meta.features
feature_info$gene <- rownames(feature_info)

if (ncol(feature_info) == 1) {
  message("No extra feature metadata found; showing gene names only.")
}

print(utils::head(feature_info))
print(colnames(raw_counts))


phenotype_data <- GSE229513_gastruloids_seurat@meta.data
tinytable::tt(
  phenotype_data |>
    head(),
  caption = "Phenotype data for the single-cell data"
)

tinytable::tt(
  phenotype_data |>
    dplyr::group_by(batch, timepoints) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(batch, timepoints),
  caption = "Number of cells per batch and timepoint"
)

tinytable::tt(
  phenotype_data |>
    dplyr::group_by(celltypeannotation) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop"),
  caption = "Number of cells per cell type"
)

GSE229513_gastruloids_seurat@meta.data$batch <- forcats::fct_recode(
  factor(GSE229513_gastruloids_seurat@meta.data$batch),
  "B-S" = "batch1",
  "SBR" = "batch2"
)

cell_type_colours <- c(
  "Anterior primitive streak/Def. endoderm" = "#faa38a",
  "Caudal epiblast"                         = "#56d312",
  "Caudal epiblast/primitive streak"        = "#ffc36a",
  "Caudal mesoderm"                         = "#01ef92",
  "Cd63+ ectoderm-like artefact"            = "#e9b000",
  "Ectopic pluripotency"                    = "#01d9bd",
  "Epiblast"                                = "#bdfe0b",
  "Epiblast/primitive streak"               = "#70cb94",
  "Exiting naïve pluripotency"              = "#efff4e",
  "Gut"                                     = "#8affc4",
  "Hemogenic endothelium"                   = "#cfba1d",
  "Naïve pluripotency"                      = "#35d365",
  "Neuromesodermal progenitors"             = "#edff9b",
  "Paraxial mesoderm"                       = "#53d240",
  "Pre-somitic mesoderm"                    = "#b8dfa2",
  "Primitive streak"                        = "#9ec72a",
  "Somite"                                  = "#78cc6e",
  "Somite differentiation front"            = "#baff73",
  "Zscan4+ Artefact"                        = "#a5c54a"
)

GSE229513_gastruloids_seurat@meta.data$celltype_colour <- cell_type_colours[
  GSE229513_gastruloids_seurat@meta.data$celltypeannotation
]
batch_shapes <- c("B-S" = 6, "SBR" = 3)

# ==========================================================================
# 2. Export PCA + UMAP for each timepoint (one PDF page per ident)
# ==========================================================================
plot_prefix <- file.path(
  output_dir,
  paste(study, "projection_plots_", today, sep = "_")
)
time_points <- levels(Seurat::Idents(GSE229513_gastruloids_seurat))
plot_list <- vector(mode = "list", length = length(time_points))
names(plot_list) <- time_points

for (tp in time_points) {
  tp_object <- subset(
    x = GSE229513_gastruloids_seurat,
    idents = tp
  )

  Seurat::Idents(tp_object) <- "celltypeannotation"

  pca_plot <- Seurat::DimPlot(
    object = tp_object,
    reduction = "pca",
    group.by = "celltypeannotation",
    cols = cell_type_colours,
    pt.size = 1.2,
    shape.by = "batch",
    label.box = TRUE,
    alpha = 0.7,
    stroke.size = 0.5
  ) +
    ggplot2::scale_shape_manual(values = batch_shapes) +
    ggplot2::ggtitle("PCA")

  umap_plot <- Seurat::DimPlot(
    object = tp_object,
    reduction = "umap",
    group.by = "celltypeannotation",
    cols = cell_type_colours,
    pt.size = 1.2,
    shape.by = "batch",
    label.box = TRUE,
    alpha = 0.7,
    stroke.size = 0.5
  ) +
    ggplot2::scale_shape_manual(values = batch_shapes) +
    ggplot2::ggtitle("UMAP")

  combined_plot <- (pca_plot | umap_plot) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = paste("Time point:", tp)
    ) &
    ggplot2::theme(legend.position = "bottom")
  plot_list[[tp]] <- patchwork::patchworkGrob(combined_plot)
}

all_timepoints_pages <- gridExtra::marrangeGrob(
  grobs = plot_list,
  nrow = 1,
  ncol = 1,
  top = NULL
)


ggplot2::ggsave(
  filename = paste0(plot_prefix, ".pdf"),
  plot = all_timepoints_pages,
  width = 14,
  height = 7,
  units = "in",
  dpi = 500
)


# ==========================================================================
# 3. Identify cell-type marker genes (phenotype features)
# ==========================================================================

# Switch default assay to SCT for differential expression
SeuratObject::DefaultAssay(GSE229513_120h) <- "SCT"

all_markers <- Seurat::FindAllMarkers(
  object = GSE229513_120h,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# --- 3a. Top 10 markers per cell type ------------------------------------
top_markers <- all_markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = avg_log2FC, n = 20) |>
  dplyr::ungroup()

# --- 3b. Summary table with tinytable ------------------------------------

marker_summary <- top_markers |>
  dplyr::select(
    `Cell type` = cluster,
    Gene = gene,
    `Log2 FC` = avg_log2FC,
    `Adj. p` = p_val_adj,
    `% Expressed (in)` = pct.1,
    `% Expressed (out)` = pct.2
  )

marker_table <- tinytable::tt(
  marker_summary,
  caption = "Top 10 marker genes per cell type (120h gastruloids)"
) |>
  tinytable::format_tt(
    j = "Log2 FC",
    fn = function(mean_signature_matrix) {
      formatC(mean_signature_matrix, digits = 2, format = "f")
    }
  ) |>
  tinytable::format_tt(
    j = "Adj. p",
    fn = function(mean_signature_matrix) {
      formatC(mean_signature_matrix, digits = 2, format = "e")
    }
  ) |>
  tinytable::style_tt(
    i = which(marker_summary$`Adj. p` < 0.01),
    bold = TRUE
  )

print(marker_table)


# ==========================================================================
# 4. Build barcode distribution table from barcodes.tsv ----
# ==========================================================================

barcode_names <- readr::read_tsv(
  "./data/raw/GSE229513_barcodes.tsv.gz",
  col_names = "barcode_raw",
  show_col_types = FALSE
) |>
  tidyr::separate_wider_regex(
    cols = barcode_raw,
    patterns = c(
      batch_id = "[A-Za-z0-9.]+",
      "[ -]",
      time_id = "[^:]+",
      ":",
      barcode_id = ".+"
    )
  )

barcodes_distribution <- barcode_names |>
  dplyr::group_by(barcode_id) |>
  dplyr::summarise(n_barcodes = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(n_barcodes))

tinytable::tt(
  barcodes_distribution,
  caption = "Distribution of barcodes per barcode ID"
)

tinytable::tt(
  barcode_names |>
    dplyr::group_by(batch_id, time_id) |>
    dplyr::arrange(batch_id, time_id) |>
    dplyr::summarise(n_barcodes = dplyr::n(), .groups = "drop"),
  caption = "Number of unique barcodes per batch and time ID combination"
)


/mnt/DATA_11TB/projects/dtoo_project/gastruloids/02_programming/Mohamad_Al_charif_DTOO/GRNITE/Data/Gastruloids/

/mnt/DATA_1/p/d/GastroDeconv2FateMap/outputs/gene_regulatory_networks/suppinger_2026-05-13_log_counts_for_grnboost_time_points.tar.gz


cp -v \
  "/mnt/DATA_11TB/projects/dtoo_project/GastroDeconv2FateMap/outputs/gene_regulatory_networks/suppinger_2026-05-13_log_counts_for_grnboost_time_points.tar.gz" \
  "/mnt/DATA_11TB/projects/dtoo_project/gastruloids/02_programming/Mohamad_Al_charif_DTOO/GRNITE/Data/Gastruloids/"
