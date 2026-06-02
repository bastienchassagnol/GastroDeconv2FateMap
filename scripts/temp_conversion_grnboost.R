study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/gene_regulatory_networks"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

source("./R/utils.R")

# Run with:
# nohup Rscript --no-save --no-restore scripts/temp_conversion_grnboost.R > logs/temp_conversion_grnboost.log 2>&1 &

# ==========================================================================
# 1. Load the Seurat object ----
# ==========================================================================

GSE229513_gastruloids_seurat <- readRDS(
  "./data/raw/GSE229513_gastruloidsobject.rds"
)

dim(GSE229513_gastruloids_seurat)

summarise_seurat_assays_layers(GSE229513_gastruloids_seurat)

# --------------------------------------------------------------------------
# 1.1 Repair object for Seurat v5 (NormalizeData / assay replacement) ----
# --------------------------------------------------------------------------
# Some objects (e.g. after integration) keep RNA counts for all genes while
# `meta.features` still rows the integrated feature set only.
Seurat::Idents(GSE229513_gastruloids_seurat) <-
  Seurat::Idents(GSE229513_gastruloids_seurat)[
    colnames(GSE229513_gastruloids_seurat)
  ]
rna_assay <- GSE229513_gastruloids_seurat[["RNA"]]
feature_names <- rownames(
  SeuratObject::LayerData(rna_assay, layer = "counts")
)
feature_df <- data.frame(matrix(nrow = length(feature_names), ncol = 0))
rownames(feature_df) <- feature_names
methods::slot(rna_assay, "meta.features") <- feature_df
GSE229513_gastruloids_seurat[["RNA"]] <- rna_assay

# ==========================================================================
# 2. CPM normalisation and log transformation ----
# ==========================================================================

# SCT transform is not recommended for GRNBoost approaches.

# 2.1 CPM normalisation, followed by Log1p transformation ----

GSE229513_gastruloids_seurat <- Seurat::NormalizeData(
  object = GSE229513_gastruloids_seurat,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# 2.2 (Optional) Gene selection for export ----
# Gold standard (Seurat): normalise first (§2.1), then rank variable genes on
# the log-normalised `RNA@data` matrix. Do not run FindVariableFeatures on raw
# counts.
#
# Choose one approach below (integration intersection OR HVG).

# 2.2.1 Integration feature set (SCT ∩ integrated ∩ RNA, ~2944 genes) ----
# sct_genes <- rownames(
#   SeuratObject::LayerData(GSE229513_gastruloids_seurat[["SCT"]], layer = "data")
# )
# integrated_genes <- rownames(
#   SeuratObject::LayerData(
#     GSE229513_gastruloids_seurat[["integrated"]],
#     layer = "data"
#   )
# )
# export_genes <- sort(intersect(
#   intersect(sct_genes, integrated_genes),
#   rownames(
#     SeuratObject::LayerData(
#       GSE229513_gastruloids_seurat[["RNA"]],
#       layer = "data"
#     )
#   )
# ))

# 2.2.2 Highly variable genes on log-normalised RNA (VST, top 3000) ----
GSE229513_gastruloids_seurat <- Seurat::FindVariableFeatures(
  object = GSE229513_gastruloids_seurat,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 3000
)
export_genes <- Seurat::VariableFeatures(
  GSE229513_gastruloids_seurat,
  assay = "RNA"
)

GSE229513_gastruloids_seurat <- subset(
  x = GSE229513_gastruloids_seurat,
  features = export_genes
)

# summarise_seurat_assays_layers(GSE229513_gastruloids_seurat)

# ==========================================================================
# 3. Save log-transformed counts for each time point ----
# ==========================================================================

time_points <- levels(Seurat::Idents(GSE229513_gastruloids_seurat))
tsv_files <- character(length(time_points))

# 3.1 Build one GRNBoost-ready TSV per time point ----
for (i in seq_along(time_points)) {
  # 3.1.1 Subset one time point ----
  tp <- time_points[[i]]
  tp_SeuratObject <- subset(
    x = GSE229513_gastruloids_seurat,
    idents = tp
  )
  # `dim(seurat)` = DefaultAssay features x cells (e.g. integrated ~2944 x n).

  # 3.1.2 Extract log-normalised RNA expression (integration gene set) ----
  tp_expression_matrix <- SeuratObject::GetAssayData(
    tp_SeuratObject,
    assay = "RNA",
    layer = "data"
  )
  # Object already subset to `export_genes` in §2.2.
  message(
    "The dimension of the RNA `data` matrix is: ",
    "Number of genes: ",
    dim(tp_expression_matrix)[1],
    ", Number of uniquely tagged cells: ",
    dim(tp_expression_matrix)[2],
    ", for time point ",
    tp,
    "."
  )

  # 3.1.3 Convert to GRNBoost format (cells as rows) ----
  tp_expression_table <- as.data.frame(as.matrix(
    Matrix::t(tp_expression_matrix)
  ))
  tp_expression_table <- cbind(
    cell_id = rownames(tp_expression_table),
    tp_expression_table
  )

  # 3.1.4 Save the time point TSV ----
  tsv_file <- file.path(
    output_dir,
    paste0(study, "_time_point_", tp, "_log_counts.tsv")
  )
  readr::write_tsv(
    tp_expression_table,
    file = tsv_file
  )

  tsv_files[[i]] <- tsv_file
}

# 3.2 Bundle all time point TSV files into one archive ----
utils::tar(
  tarfile = file.path(
    output_dir,
    paste0(
      study,
      "_",
      today,
      "_log_counts_after_hv_selection.tar.gz"
    )
  ),
  files = tsv_files,
  compression = "gzip"
)
# clean up temporary csv files

unlink(
  list.files(output_dir, pattern = "_log_counts.tsv", full.names = TRUE),
  force = TRUE
)

message(
  "\n\n",
  "Conversion complete. Output file: ",
  file.path(
    output_dir,
    paste0(
      study,
      "_",
      today,
      "_log_counts_after_hv_selection.tar.gz"
    )
  )
)

# cp -v \
#   "/mnt/DATA_11TB/projects/dtoo_project/GastroDeconv2FateMap/outputs/gene_regulatory_networks/suppinger_2026-05-13_log_counts_for_grnboost_time_points.tar.gz" \
#   "/mnt/DATA_11TB/projects/dtoo_project/gastruloids/02_programming/Mohamad_Al_charif_DTOO/GRNITE/Data/Gastruloids/"
