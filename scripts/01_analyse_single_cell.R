library(Seurat)
library(tinytable)
library(patchwork)
library(ggplot2)
library(gridExtra)

# ==========================================================================
# 0. Plotting and output settings ----
# ==========================================================================

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/single-cell"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ==========================================================================
# 1. Load the Seurat object ----
# ==========================================================================

GSE229513_gastruloids_seurat <- readRDS(
  "./data/raw-data/GSE229513_gastruloidsobject.rds"
)

dim(GSE229513_gastruloids_seurat)
table(Seurat::Idents(GSE229513_gastruloids_seurat))


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
    object    = tp_object,
    reduction = "pca",
    group.by  = "celltypeannotation",
    pt.size   = 0.6
  ) + ggplot2::ggtitle("PCA")

  umap_plot <- Seurat::DimPlot(
    object    = tp_object,
    reduction = "umap",
    group.by  = "celltypeannotation",
    pt.size   = 0.6
  ) + ggplot2::ggtitle("UMAP")

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
  dpi = 300
)



# ==========================================================================
# 5. Identify cell-type marker genes (phenotype features)
# ==========================================================================

# Switch default assay to SCT for differential expression
SeuratObject::DefaultAssay(GSE229513_120h) <- "SCT"

all_markers <- Seurat::FindAllMarkers(
  object          = GSE229513_120h,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)

# --- 5a. Top 10 markers per cell type ------------------------------------
top_markers <- all_markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = avg_log2FC, n = 20) |>
  dplyr::ungroup()


# ==========================================================================
# 6. Summary table with tinytable
# ==========================================================================

marker_summary <- top_markers |>
  dplyr::select(
    `Cell type`  = cluster,
    Gene         = gene,
    `Log2 FC`    = avg_log2FC,
    `Adj. p`     = p_val_adj,
    `% Expressed (in)` = pct.1,
    `% Expressed (out)` = pct.2
  )

marker_table <- tinytable::tt(
  marker_summary,
  caption = "Top 10 marker genes per cell type (120h gastruloids)"
) |>
  tinytable::format_tt(
    j = "Log2 FC",
    fn = function(mean_signature_matrix) formatC(mean_signature_matrix, digits = 2, format = "f")
  ) |>
  tinytable::format_tt(
    j = "Adj. p",
    fn = function(mean_signature_matrix) formatC(mean_signature_matrix, digits = 2, format = "e")
  ) |>
  tinytable::style_tt(
    i    = which(marker_summary$`Adj. p` < 0.01),
    bold = TRUE
  )

print(marker_table)
