# ==========================================================================
# 0. Libraries and Filename settings ----
# ==========================================================================

luque_merged <- readRDS(
  file = "./data/intermediate/luque_single_cell_merged_2026-06-07.rds"
)
luque_merged <- subset(luque_merged, idents = "120h")
dim(luque_merged)


# summarise phenotype data

phenotype_data <- luque_merged@meta.data
tinytable::tt(
  phenotype_data |>
    dplyr::filter(.data$timepoint == "120h") |>
    dplyr::group_by( .data$Morphotype, .data$luque_cluster_annotation, .data$Sample.barcode) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(.data$Morphotype, .data$Sample.barcode, .data$luque_cluster_annotation),
  caption = "Number of cells per time point and celltype annotation"
)
