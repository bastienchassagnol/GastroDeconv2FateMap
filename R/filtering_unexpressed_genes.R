#' Filter pseudo-bulk genes for DESeq2 by cell type
#'
#' @description
#' Aggregates raw counts by cell type and biological replicate, analogous to
#' \code{Seurat::AggregateExpression()}, then applies
#' \code{edgeR::filterByExpr()} within that cell type. This keeps the tested
#' gene universe aligned with cell-type-specific pseudo-bulk
#' \code{Morphotype} contrasts. Data-adaptive methods such as HTSFilter or IHW
#' can be added downstream, but are not applied here.
#'
#' @param seurat_obj A \pkg{Seurat} object.
#' @param celltype_col,phenotype_col,barcode_col Metadata columns containing
#'   cell-type labels, condition labels, and biological replicate identifiers.
#' @param assay Assay used for raw counts. Default \code{"RNA"}.
#' @param cell_mask Optional logical vector aligned with cells in
#'   \code{seurat_obj}; only selected cells are used.
#' @param min_cells_per_sample Minimum cells required for a cell-type x barcode
#'   pseudo-bulk sample to be retained.
#' @param min_count,min_total_count,large_n,min_prop Parameters passed to
#'   \code{edgeR::filterByExpr()}. The default \code{min_count = 5} is slightly
#'   less stringent than bulk RNA-seq defaults because cell-type pseudo-bulk
#'   libraries are often smaller.
#'
#' @return Named list. Names are cell types and values are selected gene names.
#'
#' @export
filter_deseq2_unexpressed_genes <- function(
    seurat_obj,
    celltype_col = "luque_cluster_annotation",
    phenotype_col = "Morphotype",
    barcode_col = "Sample.barcode",
    assay = "RNA",
    cell_mask = NULL,
    min_cells_per_sample = 10L,
    min_count = 5L,
    min_total_count = 15L,
    large_n = 10L,
    min_prop = 0.7
) {
  .check_filter_inputs(
    seurat_obj = seurat_obj,
    meta_cols = c(celltype_col, phenotype_col, barcode_col),
    assay = assay,
    cell_mask = cell_mask
  )

  meta <- seurat_obj@meta.data
  meta$._cell_name <- .metadata_cell_names(meta, seurat_obj, assay = assay)
  counts <- SeuratObject::LayerData(
    object = seurat_obj,
    assay = assay,
    layer = "counts"
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  if (!is.null(cell_mask)) {
    cell_mask[is.na(cell_mask)] <- FALSE
    meta <- meta[cell_mask, , drop = FALSE]
  }
  meta <- meta[
    !is.na(meta[[celltype_col]]) &
      !is.na(meta[[phenotype_col]]) &
      !is.na(meta[[barcode_col]]),
    ,
    drop = FALSE
  ]

  celltypes <- sort(unique(as.character(meta[[celltype_col]])))
  keep_by_celltype <- stats::setNames(vector("list", length(celltypes)), celltypes)

  for (ct in celltypes) {
    meta_ct <- meta[as.character(meta[[celltype_col]]) == ct, , drop = FALSE]
    meta_ct$aggregate_id <- paste(
      meta_ct[[celltype_col]],
      meta_ct[[barcode_col]],
      meta_ct[[phenotype_col]],
      sep = "__"
    )
    sample_sizes <- table(meta_ct$aggregate_id)
    keep_groups <- names(sample_sizes)[sample_sizes >= min_cells_per_sample]
    meta_ct <- meta_ct[meta_ct$aggregate_id %in% keep_groups, , drop = FALSE]

    phenotypes <- unique(as.character(meta_ct[[phenotype_col]]))
    if (nrow(meta_ct) == 0L || length(phenotypes) < 2L) {
      keep_by_celltype[[ct]] <- character()
      next
    }

    counts_ct <- counts[, meta_ct$._cell_name, drop = FALSE]
    group <- factor(meta_ct$aggregate_id)
    group_design <- Matrix::sparse.model.matrix(~ 0 + group)
    colnames(group_design) <- levels(group)
    pb_counts <- as.matrix(counts_ct %*% group_design)

    pb_meta <- unique(meta_ct[, c("aggregate_id", phenotype_col), drop = FALSE])
    rownames(pb_meta) <- pb_meta$aggregate_id
    pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]
    pb_meta[[phenotype_col]] <- stats::relevel(
      factor(pb_meta[[phenotype_col]]),
      ref = "TLS"
    )

    design <- stats::model.matrix(
      stats::as.formula(paste0("~ ", phenotype_col)),
      data = pb_meta
    )
    keep <- edgeR::filterByExpr(
      y = pb_counts,
      design = design,
      min.count = min_count,
      min.total.count = min_total_count,
      large.n = large_n,
      min.prop = min_prop
    )
    keep_by_celltype[[ct]] <- rownames(pb_counts)[keep]
  }

  keep_by_celltype
}


#' Filter single-cell genes for MAST by cell type
#'
#' @description
#' Implements a sample-aware detection filter before MAST: a gene is retained
#' for a cell type when it is detected in at least
#' \code{min_cell_detection} of cells for enough biological replicates in at
#' least one \code{Morphotype}. This avoids strong global filters that mostly
#' keep housekeeping genes.
#'
#' @inheritParams filter_deseq2_unexpressed_genes
#' @param layer Assay layer used for detection. Default \code{"data"}; use
#'   \code{"counts"} for raw-count detection.
#' @param min_cell_detection Minimum within-sample cell detection fraction.
#' @param min_sample_fraction Fraction of biological replicates per phenotype
#'   that must pass \code{min_cell_detection}.
#' @param min_total_positive_cells Minimum number of positive cells in the
#'   tested cell type.
#' @param min_cells_per_group Minimum cells required in a
#'   barcode x phenotype x cell-type group before that group contributes to the
#'   replicate count.
#'
#' @return Named list. Names are cell types and values are selected gene names.
#'
#' @export
filter_mast_unexpressed_genes <- function(
    seurat_obj,
    celltype_col = "luque_cluster_annotation",
    phenotype_col = "Morphotype",
    barcode_col = "Sample.barcode",
    assay = "RNA",
    layer = "data",
    cell_mask = NULL,
    min_cell_detection = 0.05,
    min_sample_fraction = 0.50,
    min_total_positive_cells = 20L,
    min_cells_per_group = 10L
) {
  .check_filter_inputs(
    seurat_obj = seurat_obj,
    meta_cols = c(celltype_col, phenotype_col, barcode_col),
    assay = assay,
    cell_mask = cell_mask
  )

  meta <- seurat_obj@meta.data
  meta$._cell_name <- .metadata_cell_names(meta, seurat_obj, assay = assay)
  expr <- SeuratObject::LayerData(
    object = seurat_obj,
    assay = assay,
    layer = layer
  )
  if (!inherits(expr, "dgCMatrix")) {
    expr <- Matrix::Matrix(expr, sparse = TRUE)
  }

  if (!is.null(cell_mask)) {
    cell_mask[is.na(cell_mask)] <- FALSE
    meta <- meta[cell_mask, , drop = FALSE]
    expr <- expr[, cell_mask, drop = FALSE]
  }
  meta <- meta[match(colnames(expr), meta$._cell_name), , drop = FALSE]
  valid <- !is.na(meta[[celltype_col]]) &
    !is.na(meta[[phenotype_col]]) &
    !is.na(meta[[barcode_col]])
  meta <- meta[valid, , drop = FALSE]
  expr <- expr[, valid, drop = FALSE]

  celltypes <- sort(unique(as.character(meta[[celltype_col]])))
  keep_by_celltype <- stats::setNames(vector("list", length(celltypes)), celltypes)

  for (ct in celltypes) {
    idx_ct <- as.character(meta[[celltype_col]]) == ct
    meta_ct <- meta[idx_ct, , drop = FALSE]
    expr_ct <- expr[, idx_ct, drop = FALSE]

    if (ncol(expr_ct) == 0L) {
      keep_by_celltype[[ct]] <- character()
      next
    }

    meta_ct$group <- paste(
      meta_ct[[barcode_col]],
      meta_ct[[phenotype_col]],
      sep = "__"
    )
    group_sizes <- table(meta_ct$group)
    valid_groups <- names(group_sizes)[group_sizes >= min_cells_per_group]
    keep_cell <- meta_ct$group %in% valid_groups
    meta_ct <- meta_ct[keep_cell, , drop = FALSE]
    expr_ct <- expr_ct[, keep_cell, drop = FALSE]

    phenotypes <- sort(unique(as.character(meta_ct[[phenotype_col]])))
    if (ncol(expr_ct) == 0L || length(phenotypes) < 2L) {
      keep_by_celltype[[ct]] <- character()
      next
    }

    det <- expr_ct > 0
    group <- factor(meta_ct$group)
    group_design <- Matrix::sparse.model.matrix(~ 0 + group)
    colnames(group_design) <- levels(group)

    detected_cells_by_group <- det %*% group_design
    n_cells_by_group <- Matrix::colSums(group_design)
    detection_fraction <- sweep(
      detected_cells_by_group,
      2L,
      n_cells_by_group,
      FUN = "/"
    )

    group_meta <- data.frame(
      group = levels(group),
      stringsAsFactors = FALSE
    )
    group_meta$sample_id <- sub("__.*$", "", group_meta$group)
    group_meta[[phenotype_col]] <- sub("^.*__", "", group_meta$group)

    pass_group <- detection_fraction >= min_cell_detection
    keep_by_phenotype <- lapply(phenotypes, function(ph) {
      cols_ph <- which(group_meta[[phenotype_col]] == ph)
      min_samples_ph <- ceiling(
        min_sample_fraction * length(unique(group_meta$sample_id[cols_ph]))
      )
      Matrix::rowSums(pass_group[, cols_ph, drop = FALSE]) >= min_samples_ph
    })

    keep_detection <- Reduce(`|`, keep_by_phenotype)
    keep_positive_cells <- Matrix::rowSums(det) >= min_total_positive_cells
    keep_by_celltype[[ct]] <- rownames(expr_ct)[keep_detection & keep_positive_cells]
  }

  keep_by_celltype
}


#' Validate common gene-filter inputs
#'
#' @param seurat_obj A \pkg{Seurat} object.
#' @param assay Assay used to resolve layer cell names.
#' @param meta_cols Required metadata columns.
#' @param assay Required assay name.
#' @param cell_mask Optional cell mask.
#'
#' @return \code{NULL}, invisibly.
#'
#' @keywords internal
.check_filter_inputs <- function(seurat_obj, meta_cols, assay, cell_mask = NULL) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.", call. = FALSE)
  }
  missing <- setdiff(meta_cols, colnames(seurat_obj@meta.data))
  if (length(missing) > 0L) {
    stop(
      "Missing metadata column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop("Assay `", assay, "` not found in Seurat object.", call. = FALSE)
  }
  if (!is.null(cell_mask) && length(cell_mask) != ncol(seurat_obj)) {
    stop("`cell_mask` length must match the number of cells.", call. = FALSE)
  }
  invisible(NULL)
}


#' Resolve metadata rows to Seurat cell names
#'
#' @param meta Seurat metadata.
#' @param seurat_obj A \pkg{Seurat} object.
#'
#' @return Character vector aligned with \code{meta}.
#'
#' @keywords internal
.metadata_cell_names <- function(meta, seurat_obj, assay) {
  assay_cells <- colnames(SeuratObject::LayerData(
    object = seurat_obj,
    assay = assay,
    layer = "counts"
  ))
  if (!is.null(rownames(meta)) && all(rownames(meta) %in% assay_cells)) {
    return(rownames(meta))
  }
  if ("cell_id" %in% colnames(meta) && all(meta$cell_id %in% assay_cells)) {
    return(as.character(meta$cell_id))
  }
  stop(
    "Could not resolve Seurat cell names from metadata row names or `cell_id`.",
    call. = FALSE
  )
}
