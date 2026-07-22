# ==========================================================================
# 0. Libraries and Filename settings ----
# ==========================================================================

library(Seurat)
library(tinytable)
library(patchwork)
library(ggplot2)
library(gridExtra)

study <- "luque"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/single-cell"
source("./R/utils.R")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ==========================================================================
# 1. Load and merge Seurat objects across time points ----
# ==========================================================================

# 1.1: Define timepoint paths and labels ----

timepoint_paths <- c(
  h48 = "./data/raw/GSE250136/GSM7974412_df48_final.rds.gz",
  h72 = "./data/raw/GSE250136/GSM7974413_df72_final.rds.gz",
  h96 = "./data/raw/GSE250136/GSM7974414_df96_final.rds.gz",
  h120 = "./data/raw/GSE250136/GSE250136_df120_final.rds.gz"
)
timepoint_labels <- c(
  h48 = "48h",
  h72 = "72h",
  h96 = "96h",
  h120 = "120h"
)

splitted_seurat_objects <- stats::setNames(
  lapply(names(timepoint_paths), function(id) {
    obj <- read_double_gz_rds(timepoint_paths[[id]])
    obj$timepoint <- timepoint_labels[[id]]
    obj
  }),
  names(timepoint_paths)
)


# 1.2: Identify common genes across all time points ----
common_genes <- Reduce(
  intersect,
  lapply(
    splitted_seurat_objects,
    function(obj) {
      rownames(
        SeuratObject::GetAssayData(
          object = obj,
          assay = "RNA",
          layer = "counts"
        )
      )
    }
  )
)
message(
  "Retaining ",
  length(common_genes),
  " shared genes across all time points."
)
splitted_seurat_objects <- lapply(
  splitted_seurat_objects,
  function(obj) subset(obj, features = common_genes)
)

saveRDS(
  splitted_seurat_objects,
  file = paste0(
    "./data/intermediate/",
    study,
    "_single_cell_splitted_",
    today,
    ".rds"
  )
)


# 1.3: Merge SCT objects first, then append 48h so the SCT assay is preserved where available ----
objs_with_sct <- splitted_seurat_objects[c("h72", "h96", "h120")]
GSE250136_merged <- merge(
  x = objs_with_sct$h72,
  y = objs_with_sct[c("h96", "h120")],
  add.cell.ids = c("h72", "h96", "h120"),
  project = "GSE250136_timecourse",
  merge.data = FALSE
)
GSE250136_merged <- merge(
  x = GSE250136_merged,
  y = splitted_seurat_objects$h48,
  add.cell.ids = c("", "h48"),
  project = "GSE250136_timecourse",
  merge.data = FALSE
)

# Keep a strict guardrail for downstream pseudo-bulk DESeq2:
# raw RNA counts must remain available as integer-like values after merge.
merged_rna_counts <- SeuratObject::LayerData(
  object = GSE250136_merged,
  assay = "RNA",
  layer = "counts"
)
if (is.null(merged_rna_counts)) {
  stop("Merged object is missing RNA/counts; raw-count DE analysis would fail.")
}
if (!all(abs(merged_rna_counts@x - round(merged_rna_counts@x)) < 1e-8)) {
  stop("Merged RNA/counts are not integer-like; check merge strategy.")
}

# 1.4: Finalise the merged object, adding timepoint and celltype annotations ----
SeuratObject::DefaultAssay(GSE250136_merged) <- "RNA"
Seurat::Idents(GSE250136_merged) <- "timepoint"

# 1.5 Add celltype annotation to the global merged object ----

GSE250136_merged@meta.data <- GSE250136_merged@meta.data |>
  tibble::rownames_to_column("cell_id") |>
  dplyr::mutate(
    seurat_clusters = as.character(.data$seurat_clusters)
  ) |>
  dplyr::inner_join(
    cluster_mapping,
    by = c("timepoint", "seurat_clusters")
  )

saveRDS(
  GSE250136_merged,
  file = paste0(
    "./data/intermediate/",
    study,
    "_single_cell_merged_",
    today,
    ".rds"
  )
)


GSE250136_merged <- readRDS(
  file = "./data/intermediate/luque_single_cell_merged_2026-06-07.rds"
)
phenotype_data <- GSE250136_merged@meta.data
tinytable::tt(
  phenotype_data |>
    dplyr::filter(.data$timepoint == "120h") |>
    dplyr::group_by(
      .data$Morphotype,
      .data$luque_cluster_annotation,
      .data$Sample.barcode
    ) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(
      .data$Morphotype,
      .data$luque_cluster_annotation,
      .data$Sample.barcode
    ),
  caption = "Number of cells per Morphotype, luque_cluster_annotation and Sample.barcode for the 120h time point"
)

### morphotype annotations per cell annotations per time point
tinytable::tt(
  phenotype_data |>
    dplyr::group_by(.data$timepoint, .data$Morphotype) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(.data$timepoint, .data$Morphotype),
  caption = "Number of morphotype annotations per time point"
)

source("./R/utils.R")
summarise_seurat_assays_layers(GSE250136_merged)

# ==========================================================================
# 2. Join luque_cluster_annotation onto per-time-point objects
# ==========================================================================

# 2.1 Add luque_cluster_annotation onto per-time-point objects ----

cluster_mapping <- readr::read_csv(
  "./data/intermediate/mapping_seurat_villaronga.csv",
  show_col_types = FALSE
) |>
  dplyr::mutate(
    seurat_clusters = as.character(.data$seurat_clusters)
  )

.add_luque_cluster_annotation <- function(object, timepoint_label) {
  meta <- object@meta.data |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::mutate(
      timepoint = timepoint_label,
      seurat_clusters = as.character(.data$seurat_clusters)
    ) |>
    dplyr::inner_join(
      cluster_mapping,
      by = c("timepoint", "seurat_clusters")
    )

  matched_cells <- meta$cell_id
  object <- subset(object, cells = matched_cells)
  annotation <- stats::setNames(
    meta$luque_cluster_annotation,
    meta$cell_id
  )
  object$luque_cluster_annotation <- unname(
    annotation[colnames(object)]
  )
  object
}

splitted_seurat_objects <- lapply(
  names(timepoint_labels),
  function(tp_id) {
    .add_luque_cluster_annotation(
      object = splitted_seurat_objects[[tp_id]],
      timepoint_label = timepoint_labels[[tp_id]]
    )
  }
)
names(splitted_seurat_objects) <- names(timepoint_labels)


# ==========================================================================
# 3. Export PCA + UMAP for each time point (one PDF page per time point)
# ==========================================================================

# merge() drops reductions; author PCA/UMAP live on each per-time-point object.
# Metadata (luque_cluster_annotation, Sample.barcode, Morphotype) matches
# GSE250136_merged after inner join with cluster_mapping.
annotation_ids <- sort(unique(unlist(lapply(
  splitted_seurat_objects,
  function(obj) as.character(obj$luque_cluster_annotation)
))))
annotation_colours <- stats::setNames(
  grDevices::hcl.colors(length(annotation_ids), palette = "Zissou1"),
  annotation_ids
)

morphotype_colours <- stats::setNames(
  c("neural_bias" = "#E41A1C", "TLS" = "#377EB8"),
  c("neural_bias", "TLS")
)

morphotype_shapes <- stats::setNames(
  c("neural_bias" = 16, "TLS" = 17),
  c("neural_bias", "TLS")
)

.get_meta_column <- function(object, column, fill = NA_character_) {
  meta <- object@meta.data
  if (!column %in% colnames(meta)) {
    return(rep(fill, nrow(meta)))
  }
  as.character(meta[[column]])
}

.make_projection_plot <- function(
  object,
  reduction,
  title,
  shape_by = c("Sample.barcode", "Morphotype"),
  show_shape_legend = FALSE
) {
  shape_by <- match.arg(shape_by)
  coords <- Seurat::Embeddings(object, reduction = reduction)
  plot_df <- cbind(
    as.data.frame(coords),
    luque_cluster_annotation = .get_meta_column(
      object,
      "luque_cluster_annotation"
    ),
    shape_group = .get_meta_column(object, shape_by)
  )

  colour_levels <- sort(unique(
    plot_df$luque_cluster_annotation[!is.na(plot_df$luque_cluster_annotation)]
  ))
  colour_values <- annotation_colours[colour_levels]
  colour_values[is.na(colour_values)] <- "#999999"

  shape_levels <- sort(unique(plot_df$shape_group[!is.na(plot_df$shape_group)]))
  if (shape_by == "Morphotype") {
    shape_values <- morphotype_shapes[shape_levels]
  } else {
    shape_values <- stats::setNames(
      rep(0:25, length.out = length(shape_levels)),
      shape_levels
    )
  }

  dim_x <- colnames(coords)[1L]
  dim_y <- colnames(coords)[2L]
  shape_guide <- if (show_shape_legend) {
    ggplot2::guide_legend(title = shape_by)
  } else {
    "none"
  }

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[[dim_x]],
      y = .data[[dim_y]],
      colour = .data$luque_cluster_annotation,
      shape = .data$shape_group
    )
  ) +
    ggplot2::geom_point(size = 1.2, alpha = 0.7) +
    ggplot2::scale_colour_manual(
      values = colour_values,
      name = "Luque cluster"
    ) +
    ggplot2::scale_shape_manual(
      values = shape_values,
      name = shape_by,
      drop = FALSE
    ) +
    ggplot2::guides(shape = shape_guide) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}

page_title_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
)

plot_list <- vector(mode = "list", length = length(timepoint_labels))
names(plot_list) <- unname(timepoint_labels)

for (tp_id in names(timepoint_labels)) {
  tp_label <- timepoint_labels[[tp_id]]
  tp_object <- splitted_seurat_objects[[tp_id]]

  if (tp_label == "120h") {
    combined_plot <- (.make_projection_plot(
      object = tp_object,
      reduction = "pca",
      title = "PCA",
      shape_by = "Sample.barcode"
    ) |
      .make_projection_plot(
        object = tp_object,
        reduction = "umap",
        title = "UMAP",
        shape_by = "Sample.barcode"
      )) /
      (.make_projection_plot(
        object = tp_object,
        reduction = "pca",
        title = "PCA",
        shape_by = "Morphotype",
        show_shape_legend = TRUE
      ) |
        .make_projection_plot(
          object = tp_object,
          reduction = "umap",
          title = "UMAP",
          shape_by = "Morphotype",
          show_shape_legend = TRUE
        ))
  } else {
    combined_plot <-
      .make_projection_plot(
        object = tp_object,
        reduction = "pca",
        title = "PCA",
        shape_by = "Sample.barcode"
      ) |
      .make_projection_plot(
        object = tp_object,
        reduction = "umap",
        title = "UMAP",
        shape_by = "Sample.barcode"
      )
  }

  combined_plot <- combined_plot +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = glue::glue("Time point: {tp_label}"),
      theme = page_title_theme
    ) &
    ggplot2::theme(legend.position = "bottom")
  plot_list[[tp_label]] <- patchwork::patchworkGrob(combined_plot)
}

all_timepoints_pages <- gridExtra::marrangeGrob(
  grobs = plot_list,
  nrow = 1,
  ncol = 1,
  top = NULL
)

# 3.2: Save the plots ----

plot_prefix <- file.path(
  output_dir,
  paste0(study, "_projection_plots_", today, ".pdf")
)
ggplot2::ggsave(
  filename = plot_prefix,
  plot = all_timepoints_pages,
  width = 14,
  height = 10,
  units = "in",
  dpi = 500
)

# ==========================================================================
# 4. 120h UMAP: one page per luque_cluster_annotation
# ==========================================================================

.make_120h_celltype_umap <- function(object, celltype_label) {
  coords <- Seurat::Embeddings(object, reduction = "umap")
  plot_df <- cbind(
    as.data.frame(coords),
    luque_cluster_annotation = as.character(object$luque_cluster_annotation),
    Morphotype = as.character(object$Morphotype),
    Sample.barcode = as.character(object$Sample.barcode)
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(.data$luque_cluster_annotation == celltype_label)

  morphotype_levels <- sort(unique(plot_df$Morphotype[
    !is.na(plot_df$Morphotype)
  ]))
  morphotype_values <- morphotype_colours[morphotype_levels]
  morphotype_values[is.na(morphotype_values)] <- "#999999"

  barcode_levels <- sort(unique(plot_df$Sample.barcode[
    !is.na(plot_df$Sample.barcode)
  ]))
  barcode_shapes <- stats::setNames(
    rep(0:25, length.out = length(barcode_levels)),
    barcode_levels
  )

  dim_x <- colnames(coords)[1L]
  dim_y <- colnames(coords)[2L]

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[[dim_x]],
      y = .data[[dim_y]],
      colour = .data$Morphotype,
      shape = .data$Sample.barcode
    )
  ) +
    ggplot2::geom_point(size = 1.5, alpha = 0.8) +
    ggplot2::scale_colour_manual(
      values = morphotype_values,
      name = "Morphotype",
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = barcode_shapes,
      name = "Sample.barcode",
      drop = FALSE
    ) +
    ggplot2::labs(
      title = glue::glue("120h - {celltype_label}"),
      subtitle = glue::glue("{nrow(plot_df)} cells")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.key.size = grid::unit(1.1, "lines"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12)
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        title = "Morphotype",
        override.aes = list(size = 4)
      ),
      shape = ggplot2::guide_legend(
        title = "Sample.barcode",
        ncol = 4,
        override.aes = list(size = 4)
      )
    )
}

h120_object <- splitted_seurat_objects$h120
celltype_levels <- sort(unique(as.character(
  h120_object$luque_cluster_annotation
)))

celltype_plot_list <- lapply(celltype_levels, function(ct) {
  p <- .make_120h_celltype_umap(h120_object, ct)
  ggplot2::ggplotGrob(p)
})
names(celltype_plot_list) <- celltype_levels

all_celltype_pages <- gridExtra::marrangeGrob(
  grobs = celltype_plot_list,
  nrow = 1,
  ncol = 1,
  top = NULL
)

celltype_pdf <- file.path(
  output_dir,
  paste0(study, "_120h_celltype_umap_", today, ".pdf")
)
ggplot2::ggsave(
  filename = celltype_pdf,
  plot = all_celltype_pages,
  width = 12,
  height = 9,
  units = "in",
  dpi = 500
)
