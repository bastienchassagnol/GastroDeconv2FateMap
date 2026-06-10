# ============================================================================
# Pseudo-bulk simulation from single-cell Seurat objects ----
# ============================================================================
#
# Implements strategies 1–3 and generative count models (log-normal and
# Splatter-style negative binomial) described in
# vignettes/pseudo-bulk_generation.qmd. All bulk profiles are raw summed
# counts (no normalisation).
# ============================================================================


# ============================================================================
# Internal helpers ----
# ============================================================================

#' @keywords internal
.get_sc_counts_and_meta <- function(
    seurat_obj,
    assay = "RNA"
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.", call. = FALSE)
  }
  counts <- SeuratObject::GetAssayData(
    object = seurat_obj,
    assay = assay,
    layer = "counts"
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  meta <- seurat_obj@meta.data
  list(counts = counts, meta = meta, assay = assay)
}


#' @keywords internal
.validate_pseudo_bulk_columns <- function(
    meta,
    phenotype_col,
    barcode_col
) {
  missing <- setdiff(c(phenotype_col, barcode_col), colnames(meta))
  if (length(missing) > 0L) {
    stop(
      "Missing metadata column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyNA(meta[[phenotype_col]]) || anyNA(meta[[barcode_col]])) {
    stop(
      "Phenotype and barcode columns must not contain NA values.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' @keywords internal
.balanced_phenotype_labels <- function(
    n_samples,
    phenotypes
) {
  if (n_samples %% 2L != 0L) {
    stop(
      "`n_samples` must be even for balanced morphotype allocation.",
      call. = FALSE
    )
  }
  if (length(phenotypes) < 2L) {
    stop("At least two phenotypes are required.", call. = FALSE)
  }
  n_per_pheno <- n_samples %/% length(phenotypes)
  rep(phenotypes, each = n_per_pheno)
}


#' @keywords internal
.default_cells_per_sample <- function(meta, barcode_col) {
  as.integer(stats::median(
    table(meta[[barcode_col]]),
    na.rm = TRUE
  ))
}


#' @keywords internal
.get_feature_metadata <- function(seurat_obj, assay = "RNA") {
  assay_obj <- seurat_obj[[assay]]
  feat <- assay_obj@meta.features
  if (is.null(feat) || nrow(feat) == 0L) {
    feat <- S4Vectors::DataFrame(
      gene = rownames(assay_obj),
      row.names = rownames(assay_obj)
    )
  }
  feat
}


#' @keywords internal
.build_pseudo_bulk_se <- function(
    count_matrix,
    col_data,
    feature_data
) {

  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = count_matrix),
    rowData = feature_data,
    colData = S4Vectors::DataFrame(col_data)
  )
}


#' @keywords internal
.aggregate_cells_to_bulk <- function(counts, cell_idx) {
  if (length(cell_idx) == 0L) {
    stop("No cells selected for aggregation.", call. = FALSE)
  }
  Matrix::rowSums(counts[, cell_idx, drop = FALSE])
}


#' @keywords internal
.infer_lognormal_gene_params <- function(
    counts,
    cell_idx,
    max_cells = 2000L
) {
  idx <- cell_idx
  if (length(idx) > max_cells) {
    idx <- sample(idx, max_cells)
  }
  log_x <- log1p(as.matrix(counts[, idx, drop = FALSE]))
  list(
    mean_log = rowMeans(log_x),
    sd_log = matrixStats::rowSds(log_x)
  )
}


#' @keywords internal
.mom_nb_size <- function(x) {
  m <- mean(x)
  if (m <= 0) {
    return(Inf)
  }
  v <- stats::var(x)
  if (is.na(v) || v <= m) {
    return(Inf)
  }
  m^2 / (v - m)
}


#' @keywords internal
.infer_nb_gene_params <- function(
    counts,
    cell_idx,
    max_cells = 2000L
) {
  idx <- cell_idx
  if (length(idx) > max_cells) {
    idx <- sample(idx, max_cells)
  }
  x <- as.matrix(counts[, idx, drop = FALSE])
  mu <- rowMeans(x)
  size <- apply(x, 1L, .mom_nb_size)
  size[!is.finite(size) | size <= 0] <- 1e6
  list(mu = mu, size = size)
}


#' @keywords internal
.infer_splatter_params <- function(counts, cell_idx, max_cells = 2000L) {
  idx <- cell_idx
  if (length(idx) > max_cells) {
    idx <- sample(idx, max_cells)
  }
  mat <- as.matrix(counts[, idx, drop = FALSE])
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = mat)
  )
  splatter::splatEstimate(sce)
}


#' @keywords internal
.infer_hierarchical_params <- function(
    counts,
    meta,
    phenotype_col,
    barcode_col,
    model = c("lognormal", "negative_binomial")
) {
  model <- match.arg(model)
  phenotypes <- sort(unique(as.character(meta[[phenotype_col]])))
  barcodes <- sort(unique(as.character(meta[[barcode_col]])))

  biological <- stats::setNames(vector("list", length(phenotypes)), phenotypes)
  for (ph in phenotypes) {
    idx <- which(meta[[phenotype_col]] == ph)
    biological[[ph]] <- if (model == "lognormal") {
      .infer_lognormal_gene_params(counts, idx)
    } else {
      list(
        gene = .infer_nb_gene_params(counts, idx),
        splatter = .infer_splatter_params(counts, idx)
      )
    }
  }

  lib_sizes <- stats::setNames(
    vapply(
      barcodes,
      function(bc) {
        idx <- which(meta[[barcode_col]] == bc)
        sum(Matrix::colSums(counts[, idx, drop = FALSE]))
      },
      numeric(1L)
    ),
    barcodes
  )
  global_lib <- mean(lib_sizes)
  technical <- list(
    library_size = lib_sizes,
    library_scale = lib_sizes / global_lib
  )

  list(biological = biological, technical = technical)
}


#' @keywords internal
.simulate_lognormal_cells <- function(
    bio_params,
    tech_scale,
    n_cells,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  mu <- bio_params$mean_log
  sd <- pmax(bio_params$sd_log, 0.01)
  n_genes <- length(mu)
  z <- stats::rnorm(
    n = n_genes * n_cells,
    mean = rep(mu, each = n_cells),
    sd = rep(sd, each = n_cells)
  )
  mat <- matrix(z, nrow = n_genes, ncol = n_cells)
  out <- pmax(0, round(expm1(mat) * tech_scale))
  matrix(out, nrow = n_genes, ncol = n_cells)
}


#' @keywords internal
.simulate_nb_cells <- function(
    bio_params,
    tech_scale,
    n_cells,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  mu <- bio_params$gene$mu * tech_scale
  size <- bio_params$gene$size
  n_genes <- length(mu)
  draws <- stats::rnbinom(
    n = n_genes * n_cells,
    size = rep(size, each = n_cells),
    mu = rep(mu, each = n_cells)
  )
  matrix(draws, nrow = n_genes, ncol = n_cells)
}


#' @keywords internal
.pick_anchor_barcode <- function(meta, phenotype_col, barcode_col, phenotype) {
  barcodes <- unique(as.character(meta[[barcode_col]][
    meta[[phenotype_col]] == phenotype
  ]))
  sample(barcodes, size = 1L)
}


# ============================================================================
# Strategy 1: direct barcode aggregation ----
# ============================================================================

#' Aggregate pseudo-bulk samples by barcode (organoid)
#'
#' @description
#' Strategy 1 from \code{vignettes/pseudo-bulk_generation.qmd}: treat each
#' barcode as one bulk sample by summing raw single-cell counts within
#' barcodes, analogous to \code{\link[Seurat]{AggregateExpression}} with
#' \code{slot = "counts"} and no normalisation.
#'
#' @param seurat_obj A \pkg{Seurat} object with raw counts in the RNA assay.
#' @param phenotype_col Metadata column for organoid-level phenotype labels
#'   (e.g. \code{"Morphotype"}).
#' @param barcode_col Metadata column for barcode / organoid identifiers
#'   (e.g. \code{"Sample.barcode"}).
#' @param assay Assay name from which raw counts are read. Default
#'   \code{"RNA"}.
#'
#' @return A \pkg{SummarizedExperiment} with one column per distinct barcode.
#'   \code{colData} contains \code{barcode_id}, the phenotype label,
#'   \code{library_depth} (column sum of counts), and
#'   \code{simulation_method = "barcode_aggregation"}.
#'
#' @seealso [simulate_bootstrap_samples()], [simulate_generative_models()],
#'   \code{vignettes/pseudo-bulk_generation.qmd}
#'
#' @export
aggregate_barcode_pseudo_bulk <- function(
    seurat_obj,
    phenotype_col = "Morphotype",
    barcode_col = "Sample.barcode",
    assay = "RNA"
) {
  sc <- .get_sc_counts_and_meta(seurat_obj, assay = assay)
  .validate_pseudo_bulk_columns(sc$meta, phenotype_col, barcode_col)

  barcodes <- sort(unique(as.character(sc$meta[[barcode_col]])))
  n_bar <- length(barcodes)
  bulk_mat <- matrix(
    0,
    nrow = nrow(sc$counts),
    ncol = n_bar,
    dimnames = list(rownames(sc$counts), barcodes)
  )

  barcode_id <- character(n_bar)
  phenotype <- character(n_bar)
  library_depth <- numeric(n_bar)

  for (i in seq_along(barcodes)) {
    bc <- barcodes[[i]]
    idx <- which(sc$meta[[barcode_col]] == bc)
    bulk_mat[, i] <- .aggregate_cells_to_bulk(sc$counts, idx)
    barcode_id[[i]] <- bc
    phenotype[[i]] <- as.character(sc$meta[[phenotype_col]][idx[[1L]]])
    library_depth[[i]] <- sum(bulk_mat[, i])
  }

  col_data <- data.frame(
    barcode_id = barcode_id,
    phenotype = phenotype,
    library_depth = library_depth,
    simulation_method = "barcode_aggregation",
    row.names = barcodes,
    stringsAsFactors = FALSE
  )
  names(col_data)[names(col_data) == "phenotype"] <- phenotype_col

  .build_pseudo_bulk_se(
    count_matrix = bulk_mat,
    col_data = col_data,
    feature_data = .get_feature_metadata(seurat_obj, assay = assay)
  )
}


# ============================================================================
# Strategies 2–3: bootstrap pseudo-bulk ----
# ============================================================================

#' Simulate pseudo-bulk samples by cell resampling (bootstrap)
#'
#' @description
#' Strategies 2 and 3 from \code{vignettes/pseudo-bulk_generation.qmd}.
#' \describe{
#'   \item{\code{replicate_type = "biological"}}{
#'     Phenotype-stratified bootstrap: resample cells from all organoids
#'     sharing the target phenotype (strategy 2).}
#'   \item{\code{replicate_type = "technical"}}{
#'     Organoid-preserving bootstrap: resample cells within a single barcode
#'     before aggregation (strategy 3).}
#' }
#'
#' @param seurat_obj A \pkg{Seurat} object with raw counts.
#' @param phenotype_col Metadata column for phenotype labels.
#' @param barcode_col Metadata column for barcode identifiers.
#' @param n_samples Total number of pseudo-bulk samples to generate. Must be
#'   even; half are assigned to each phenotype level.
#' @param replicate_type \code{"biological"} (phenotype-level pooling) or
#'   \code{"technical"} (within-barcode resampling).
#' @param cells_per_sample Number of cells aggregated per pseudo-bulk sample.
#'   Defaults to the median barcode cell count in the object.
#' @param assay Assay from which raw counts are read.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A \pkg{SummarizedExperiment} with raw summed counts and balanced
#'   phenotypes in \code{colData}.
#'
#' @seealso [aggregate_barcode_pseudo_bulk()], [simulate_generative_models()]
#'
#' @export
simulate_bootstrap_samples <- function(
    seurat_obj,
    phenotype_col = "Morphotype",
    barcode_col = "Sample.barcode",
    n_samples = 100L,
    replicate_type = c("biological", "technical"),
    cells_per_sample = NULL,
    assay = "RNA",
    seed = NULL
) {
  replicate_type <- match.arg(replicate_type)
  sc <- .get_sc_counts_and_meta(seurat_obj, assay = assay)
  .validate_pseudo_bulk_columns(sc$meta, phenotype_col, barcode_col)

  phenotypes <- sort(unique(as.character(sc$meta[[phenotype_col]])))
  pheno_labels <- .balanced_phenotype_labels(n_samples, phenotypes)

  if (is.null(cells_per_sample)) {
    cells_per_sample <- .default_cells_per_sample(sc$meta, barcode_col)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  bulk_mat <- matrix(
    0,
    nrow = nrow(sc$counts),
    ncol = n_samples,
    dimnames = list(
      rownames(sc$counts),
      paste0("bootstrap_", seq_len(n_samples))
    )
  )

  barcode_id <- character(n_samples)
  library_depth <- numeric(n_samples)

  for (i in seq_len(n_samples)) {
    ph <- pheno_labels[[i]]
    if (replicate_type == "biological") {
      pool <- which(sc$meta[[phenotype_col]] == ph)
      barcode_id[[i]] <- "pooled"
    } else {
      bc <- .pick_anchor_barcode(
        sc$meta, phenotype_col, barcode_col, ph
      )
      pool <- which(
        sc$meta[[barcode_col]] == bc &
          sc$meta[[phenotype_col]] == ph
      )
      barcode_id[[i]] <- bc
    }
    sampled <- sample(
      pool,
      size = cells_per_sample,
      replace = TRUE
    )
    bulk_mat[, i] <- .aggregate_cells_to_bulk(sc$counts, sampled)
    library_depth[[i]] <- sum(bulk_mat[, i])
  }

  method_tag <- if (replicate_type == "biological") {
    "bootstrap_biological"
  } else {
    "bootstrap_technical"
  }

  col_data <- data.frame(
    barcode_id = barcode_id,
    library_depth = library_depth,
    simulation_method = method_tag,
    replicate_type = replicate_type,
    row.names = colnames(bulk_mat),
    stringsAsFactors = FALSE
  )
  col_data[[phenotype_col]] <- pheno_labels

  .build_pseudo_bulk_se(
    count_matrix = bulk_mat,
    col_data = col_data,
    feature_data = .get_feature_metadata(seurat_obj, assay = assay)
  )
}


# ============================================================================
# Generative models: log-normal and Splatter-style NB ----
# ============================================================================

#' Simulate pseudo-bulk samples from generative count models
#'
#' @description
#' Splatter- and HADACA3-inspired generative simulation (strategies 7 and 10
#' in \code{vignettes/pseudo-bulk_generation.qmd}). Parameters are inferred at
#' two levels:
#' \enumerate{
#'   \item **Biological** (phenotype): per-gene means and dispersions within
#'     each \code{phenotype_col} level (\code{splatter::splatEstimate} for the
#'     negative-binomial path).
#'   \item **Technical** (barcode): per-organoid library-size scaling factors.
#' }
#' Simulated single cells are drawn from the inferred model, then summed to
#' sample-level raw bulk counts.
#'
#' @param seurat_obj A \pkg{Seurat} object with raw counts.
#' @param phenotype_col Metadata column for phenotype labels.
#' @param barcode_col Metadata column for barcode identifiers.
#' @param n_samples Total pseudo-bulk samples (even; balanced across
#'   phenotypes).
#' @param model \code{"lognormal"} for log-normal gene counts, or
#'   \code{"negative_binomial"} for Splatter-style negative-binomial counts.
#' @param cells_per_sample Cells simulated and summed per pseudo-bulk sample.
#' @param assay Assay from which raw counts are read.
#' @param seed Optional random seed.
#'
#' @return A \pkg{SummarizedExperiment} with raw summed counts.
#'
#' @details
#' Negative-binomial gene-wise dispersions use method-of-moments estimates per
#' phenotype, complemented by \code{splatter::splatEstimate} on a subsample of
#' cells for library-size structure. Log-normal parameters are estimated on
#' \code{log1p(counts)} per phenotype.
#'
#' @seealso [aggregate_barcode_pseudo_bulk()], [simulate_bootstrap_samples()],
#'   \url{https://oshlacklab.com/splatter/articles/splatter.html}
#'
#' @export
simulate_generative_models <- function(
    seurat_obj,
    phenotype_col = "Morphotype",
    barcode_col = "Sample.barcode",
    n_samples = 100L,
    model = c("lognormal", "negative_binomial"),
    cells_per_sample = NULL,
    assay = "RNA",
    seed = NULL
) {
  model <- match.arg(model)
  sc <- .get_sc_counts_and_meta(seurat_obj, assay = assay)
  .validate_pseudo_bulk_columns(sc$meta, phenotype_col, barcode_col)

  phenotypes <- sort(unique(as.character(sc$meta[[phenotype_col]])))
  pheno_labels <- .balanced_phenotype_labels(n_samples, phenotypes)

  if (is.null(cells_per_sample)) {
    cells_per_sample <- .default_cells_per_sample(sc$meta, barcode_col)
  }

  params <- .infer_hierarchical_params(
    counts = sc$counts,
    meta = sc$meta,
    phenotype_col = phenotype_col,
    barcode_col = barcode_col,
    model = model
  )

  bulk_mat <- matrix(
    0,
    nrow = nrow(sc$counts),
    ncol = n_samples,
    dimnames = list(
      rownames(sc$counts),
      paste0(model, "_", seq_len(n_samples))
    )
  )

  barcode_id <- character(n_samples)
  library_depth <- numeric(n_samples)

  for (i in seq_len(n_samples)) {
    ph <- pheno_labels[[i]]
    bc <- .pick_anchor_barcode(
      sc$meta, phenotype_col, barcode_col, ph
    )
    tech_scale <- params$technical$library_scale[[bc]]
    bio <- params$biological[[ph]]

    cell_seed <- if (is.null(seed)) NULL else seed + i
    cell_mat <- if (model == "lognormal") {
      .simulate_lognormal_cells(
        bio_params = bio,
        tech_scale = tech_scale,
        n_cells = cells_per_sample,
        seed = cell_seed
      )
    } else {
      .simulate_nb_cells(
        bio_params = bio,
        tech_scale = tech_scale,
        n_cells = cells_per_sample,
        seed = cell_seed
      )
    }

    bulk_mat[, i] <- rowSums(cell_mat)
    barcode_id[[i]] <- bc
    library_depth[[i]] <- sum(bulk_mat[, i])
  }

  method_tag <- if (model == "lognormal") {
    "generative_lognormal"
  } else {
    "generative_negative_binomial"
  }

  col_data <- data.frame(
    barcode_id = barcode_id,
    library_depth = library_depth,
    simulation_method = method_tag,
    generative_model = model,
    row.names = colnames(bulk_mat),
    stringsAsFactors = FALSE
  )
  col_data[[phenotype_col]] <- pheno_labels

  .build_pseudo_bulk_se(
    count_matrix = bulk_mat,
    col_data = col_data,
    feature_data = .get_feature_metadata(seurat_obj, assay = assay)
  )
}
