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

# --- 3h. Graph structures (SNN / NN) -------------------------------------
nn_graph  <- GSE229513_120h@graphs$integrated_nn
snn_graph <- GSE229513_120h@graphs$integrated_snn
