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
  h48  = "./data/raw/GSE250136/GSM7974412_df48_final.rds.gz",
  h72  = "./data/raw/GSE250136/GSM7974413_df72_final.rds.gz",
  h96  = "./data/raw/GSE250136/GSM7974414_df96_final.rds.gz",
  h120 = "./data/raw/GSE250136/GSE250136_df120_final.rds.gz"
)
timepoint_labels <- c(
  h48  = "48h",
  h72  = "72h",
  h96  = "96h",
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
  "Retaining ", length(common_genes),
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

# merge() keeps shared gene names but drops RNA@meta.features (0 columns on merged
# object). Recover VST feature metadata from the pre-merge objects instead.
# feature_info <- dplyr::bind_rows(lapply(
#   splitted_seurat_objects,
#   function(obj) {
#     meta <- obj@assays$RNA@meta.features
#     data.frame(
#       gene = rownames(meta),
#       meta,
#       timepoint = obj$timepoint[1L],
#       row.names = NULL
#     )
#   }
# ))
# feature_info <- feature_info |>
#   dplyr::filter(.data$gene %in% common_genes)

# tinytable::tt(
#   feature_info |>
#     utils::head(),
#   caption = paste(
#     "Gene feature information"
#   )
# )

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

# 1.4: Finalise the merged object, adding timepoint and celltype annotations ----
SeuratObject::DefaultAssay(GSE250136_merged) <- "RNA"
Seurat::Idents(GSE250136_merged) <- "timepoint"
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

dvc add data/intermediate/luque_single_cell_merged_2026-06-07.rds \
        data/intermediate/luque_single_cell_splitted_2026-06-07.rds
dvc push

dim(GSE250136_merged)
summarise_seurat_assays_layers(GSE250136_merged)



phenotype_data <- GSE250136_merged@meta.data
tinytable::tt(
  phenotype_data |>
    utils::head(),
  caption = "Phenotype data for the merged single-cell object"
)

tinytable::tt(
  phenotype_data |>
    dplyr::group_by(.data$timepoint, .data$seurat_clusters) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(.data$timepoint, dplyr::desc(.data$n_cells)),
  caption = "Number of cells per time point and Seurat cluster"
)

### morphotype annotations per cell annotations per time point
tinytable::tt(
  phenotype_data |>
    dplyr::group_by(.data$timepoint, .data$Morphotype) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(.data$timepoint, .data$Morphotype),
  caption = "Number of morphotype annotations per time point"
)

tinytable::tt(
  phenotype_data |>
    dplyr::group_by(.data$timepoint, .data$Morphotype, .data$Sample.barcode) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(.data$timepoint, .data$Sample.barcode, .data$Morphotype)  |> 
    tidyr::drop_na(),
  caption = "Number of morphotype annotations per time point and sample barcode"
)


# ==========================================================================
# 2. Export PCA + UMAP for each time point (one PDF page per time point)
# ==========================================================================

# merge() drops reductions; author PCA/UMAP live on each per-time-point object.
# Metadata (seurat_clusters, Sample.barcode, Morphotype) matches GSE250136_merged.
cluster_ids <- sort(unique(as.character(GSE250136_merged$seurat_clusters)))
cluster_colours <- stats::setNames(
  grDevices::hcl.colors(length(cluster_ids), palette = "Zissou1"),
  cluster_ids
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
    show_shape_legend = FALSE) {
  shape_by <- match.arg(shape_by)
  coords <- Seurat::Embeddings(object, reduction = reduction)
  plot_df <- cbind(
    as.data.frame(coords),
    seurat_clusters = .get_meta_column(object, "seurat_clusters"),
    shape_group = .get_meta_column(object, shape_by)
  )

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
      colour = .data$seurat_clusters,
      shape = .data$shape_group
    )
  ) +
    ggplot2::geom_point(size = 1.2, alpha = 0.7) +
    ggplot2::scale_colour_manual(
      values = cluster_colours,
      name = "Seurat cluster"
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
    combined_plot <- (
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
    ) / (
      .make_projection_plot(
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
        )
    )
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

# 2.2: Save the plots ----

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
