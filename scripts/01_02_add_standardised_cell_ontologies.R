# ==========================================================================
# 0. Settings ----
# ==========================================================================

study_today <- format(Sys.Date(), "%Y-%m-%d")
intermediate_dir <- "data/intermediate"
mapping_dir <- file.path(intermediate_dir, "ontology_mapping")

dir.create(mapping_dir, recursive = TRUE, showWarnings = FALSE)

suppinger_rds <- file.path(
  intermediate_dir,
  "suppinger_single_cell_2026-05-06.rds"
)
luque_rds <- file.path(
  intermediate_dir,
  "luque_single_cell_merged_2026-06-07.rds"
)

suppinger_unmapped_h5ad <- file.path(mapping_dir, "suppinger_unmapped.h5ad")
luque_unmapped_h5ad <- file.path(mapping_dir, "luque_unmapped.h5ad")
suppinger_harmonized_h5ad <- file.path(mapping_dir, "suppinger_harmonized.h5ad")
luque_harmonized_h5ad <- file.path(mapping_dir, "luque_harmonized.h5ad")

suppinger_out_rds <- file.path(
  intermediate_dir,
  glue::glue("suppinger_single_cell_annotated_{study_today}.rds")
)
luque_out_rds <- file.path(
  intermediate_dir,
  glue::glue("luque_single_cell_annotated_{study_today}.rds")
)

utf8_obs <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.character(x)) {
    return(enc2utf8(x))
  }
  x
}

ascii_obs <- function(x) {
  iconv(as.character(x), from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
}

prepare_export_sce <- function(sce, dataset_label) {
  cells <- colnames(sce)
  cd <- S4Vectors::DataFrame(
    cell_type_original = ascii_obs(
      SummarizedExperiment::colData(sce)$celltypeannotation
    ),
    timepoint_original = as.character(
      SummarizedExperiment::colData(sce)$timepoints
    ),
    dataset = dataset_label,
    row.names = cells
  )
  export <- sce
  SummarizedExperiment::colData(export) <- cd
  export
}

prepare_export_seurat <- function(obj, dataset_label) {
  cells <- colnames(obj)
  export <- obj
  export@meta.data <- data.frame(
    cell_type_original = ascii_obs(obj$luque_cluster_annotation),
    timepoint_original = as.character(obj$timepoint),
    dataset = dataset_label,
    row.names = cells
  )
  export
}

ontology_cols <- c(
  "cell_type_ontology_name",
  "cell_type_ontology_term_id",
  "cell_type_harmonized",
  "cell_state_standardized",
  "cell_type_mapping_relation",
  "cell_type_mapping_note",
  "annotation_artefact",
  "cell_line_original",
  "cell_line_ontology_name",
  "cell_line_ontology_term_id"
)

merge_ontology_metadata <- function(obj, h5ad_path) {
  harm <- anndataR::read_h5ad(h5ad_path, as = "SingleCellExperiment")
  harm_obs <- as.data.frame(SummarizedExperiment::colData(harm))
  cells <- colnames(obj)
  missing <- setdiff(cells, rownames(harm_obs))
  if (length(missing) > 0L) {
    stop(
      length(missing),
      " cells in the R object are missing from harmonised H5AD metadata."
    )
  }
  for (col in ontology_cols) {
    obj[[col]] <- harm_obs[cells, col, drop = TRUE]
  }
  obj
}

# ==========================================================================
# 1. Load objects ----
# ==========================================================================

suppinger_single_cell <- readRDS(suppinger_rds)
luque_single_cell_seurat <- readRDS(luque_rds)

# Merged Luque object colnames are numeric indices; RNA counts use cell barcodes.
realign_luque_seurat <- function(obj) {
  counts <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
  meta <- obj@meta.data
  rownames(meta) <- colnames(counts)
  Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
}

if (
  !identical(
    colnames(luque_single_cell_seurat),
    colnames(SeuratObject::GetAssayData(
      luque_single_cell_seurat,
      assay = "RNA",
      layer = "counts"
    ))
  )
) {
  luque_single_cell_seurat <- realign_luque_seurat(luque_single_cell_seurat)
}

# ==========================================================================
# 2. Export unmapped H5AD files (anndataR) ----
# ==========================================================================

if (!inherits(suppinger_single_cell, "SingleCellExperiment")) {
  suppinger_single_cell <- SingleCellExperiment::SingleCellExperiment(
    assays = SummarizedExperiment::assays(suppinger_single_cell),
    rowData = SummarizedExperiment::rowData(suppinger_single_cell),
    colData = SummarizedExperiment::colData(suppinger_single_cell),
    metadata = S4Vectors::metadata(suppinger_single_cell)
  )
}

SummarizedExperiment::colData(suppinger_single_cell)$cell_type_original <-
  utf8_obs(
    SummarizedExperiment::colData(suppinger_single_cell)$celltypeannotation
  )
SummarizedExperiment::colData(suppinger_single_cell)$timepoint_original <-
  as.character(
    SummarizedExperiment::colData(suppinger_single_cell)$timepoints
  )
SummarizedExperiment::colData(suppinger_single_cell)$dataset <- "Suppinger_2026"

luque_single_cell_seurat$cell_type_original <- utf8_obs(
  luque_single_cell_seurat$luque_cluster_annotation
)
luque_single_cell_seurat$timepoint_original <-
  as.character(luque_single_cell_seurat$timepoint)
luque_single_cell_seurat$dataset <- "Luque_2024"

suppinger_export <- prepare_export_sce(
  suppinger_single_cell,
  "Suppinger_2026"
)
luque_export <- prepare_export_seurat(
  luque_single_cell_seurat,
  "Luque_2024"
)

suppinger_x_assay <- if (
  "counts" %in% SummarizedExperiment::assayNames(suppinger_export)
) {
  "counts"
} else {
  SummarizedExperiment::assayNames(suppinger_export)[1]
}

anndataR::write_h5ad(
  suppinger_export,
  suppinger_unmapped_h5ad,
  x_mapping = suppinger_x_assay,
  obs_mapping = TRUE,
  var_mapping = TRUE,
  layers_mapping = FALSE,
  mode = "w"
)

anndataR::write_h5ad(
  luque_export,
  luque_unmapped_h5ad,
  assay_name = "RNA",
  x_mapping = "counts",
  obs_mapping = TRUE,
  var_mapping = TRUE,
  obsm_mapping = FALSE,
  mode = "w"
)

message("Exported unmapped H5AD files to ", mapping_dir)

# ==========================================================================
# 4. Merge harmonised metadata and save Seurat objects ----
# ==========================================================================

suppinger_harmonized_seurat <- Seurat::as.Seurat(suppinger_single_cell)
suppinger_harmonized_seurat <- merge_ontology_metadata(
  suppinger_harmonized_seurat,
  suppinger_harmonized_h5ad
)

luque_harmonized_seurat <- merge_ontology_metadata(
  luque_single_cell_seurat,
  luque_harmonized_h5ad
)

saveRDS(suppinger_harmonized_seurat, suppinger_out_rds)
saveRDS(luque_harmonized_seurat, luque_out_rds)

message("Annotation with BioOnty achieved")


suppinger_annotation <- readr::read_csv(
  "data/intermediate/suppinger_cell_type_annotations.csv"
)

tinytable::tt(
  suppinger_annotation |>
    dplyr::select(
      -.data$cell_barcode,
      -.data$dataset,
      .data$timepoint_original,
      .data$annotation_artefact
    ) |>
    dplyr::distinct(
      .data$cell_type_original,
      .data$cell_type_ontology_name,
      .keep_all = TRUE
    ),
  caption = "Unique cell types in Suppinger 2026"
)

luque_annotation <- readr::read_csv(
  "data/intermediate/luque_cell_type_annotations.csv"
)
