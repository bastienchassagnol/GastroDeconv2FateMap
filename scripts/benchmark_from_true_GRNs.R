
library(Seurat)


# ==========================================================================
# 1. Load the Seurat object ----
# ==========================================================================

GSE229513_gastruloids_seurat <- readRDS(
  "./data/raw-data/GSE229513_gastruloidsobject.rds"
)

dim(GSE229513_gastruloids_seurat)
levels(mean_signature_matrix = GSE229513_gastruloids_seurat)
table(Seurat::Idents(GSE229513_gastruloids_seurat))


# ==========================================================================
# 2. Subset to the 120h timepoint
# ==========================================================================

GSE229513_120h <- Seurat::subset(
  mean_signature_matrix = GSE229513_gastruloids_seurat,
  idents = "120h"
)

saveRDS(
  object = GSE229513_120h,
  file   = "./data/raw-data/GSE229513_120h.rds.gz",
  compress = "gzip"
)


# ==========================================================================
# 3. Extract main layers from the Seurat object
# ==========================================================================

# --- 3a. Raw counts (RNA assay) -------------------------------------------
raw_counts <- SeuratObject::GetAssayData(
  object = GSE229513_120h,
  assay  = "RNA",
  layer   = "counts"
)
cat("Raw counts (RNA):    ", nrow(raw_counts), "genes mean_signature_matrix",
    ncol(raw_counts), "cells\n")

# --- 3b. Normalised data (RNA assay) --------------------------------------
normalised_data <- SeuratObject::GetAssayData(
  object = GSE229513_120h,
  assay  = "RNA",
  layer   = "data"
)

# --- 3c. SCT-corrected counts ---------------------------------------------
sct_counts <- SeuratObject::GetAssayData(
  object = GSE229513_120h,
  assay  = "SCT",
  layer   = "counts"
)

# --- 3d. SCT log-normalised data ------------------------------------------
sct_data <- SeuratObject::GetAssayData(
  object = GSE229513_120h,
  assay  = "SCT",
  layer   = "data"
)

# > dim(sct_data)
# [1] 2944 2819

# --- 3e. Scaled data (integrated assay, used for PCA/clustering) ----------
scaled_data <- SeuratObject::GetAssayData(
  object = GSE229513_120h,
  assay  = "integrated",
  layer   = "scale.data"
)

# --- 3f. Cell-level metadata ----------------------------------------------
cell_metadata <- GSE229513_120h@meta.data
str(cell_metadata)

# --- 3g. Pre-computed embeddings ------------------------------------------
pca_embeddings  <- SeuratObject::Embeddings(GSE229513_120h, reduction = "pca")
umap_embeddings <- SeuratObject::Embeddings(GSE229513_120h, reduction = "umap")

cat("PCA dims:  ", ncol(pca_embeddings), "\n")
cat("UMAP dims: ", ncol(umap_embeddings), "\n")

# --- 3h. Graph structures (SNN / NN) -------------------------------------
nn_graph  <- GSE229513_120h@graphs$integrated_nn
snn_graph <- GSE229513_120h@graphs$integrated_snn


# ==========================================================================
# 4. Visualise PCA and UMAP coloured by cell type
#    (uses pre-computed embeddings already stored in the object)
# ==========================================================================

# Set identity to the cell-type annotation column
Seurat::Idents(GSE229513_120h) <- "celltypeannotation"

pca_plot <- Seurat::DimPlot(
  object    = GSE229513_120h,
  reduction = "pca",
  group.by  = "celltypeannotation",
  pt.size   = 0.6
) + ggplot2::ggtitle("PCA — 120h gastruloids by cell type")

umap_plot <- Seurat::DimPlot(
  object    = GSE229513_120h,
  reduction = "umap",
  group.by  = "celltypeannotation",
  pt.size   = 0.6
) + ggplot2::ggtitle("UMAP — 120h gastruloids by cell type")

print(pca_plot)
print(umap_plot)


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

