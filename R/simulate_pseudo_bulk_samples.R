# ============================================================================
# Pseudo-bulk simulation from single-cell Seurat objects ----
# ============================================================================
#
# Implements strategies 1–3 and generative count models (log-normal and
# Splatter-style negative binomial) described in
# vignettes/pseudo_bulk_generation.qmd. All bulk profiles are raw summed
# counts (no normalisation).
# ============================================================================




# ============================================================================
# Internal helpers ----
# ============================================================================

#' Balanced phenotype labels for simulated samples
#'
#' @description
#' Allocates \code{n_samples} across phenotype levels as evenly as possible.
#' When \code{n_samples} is not divisible by the number of phenotypes, the
#' first \code{n_samples %% K} levels receive one extra sample each.
#'
#' @param n_samples Total number of pseudo-bulk samples to generate.
#' @param phenotypes Character vector of phenotype levels.
#'
#' @return Character vector of length \code{n_samples} with (near-)balanced
#'   phenotype labels.
#'
#' @keywords internal
#' @importFrom forcats fct_inorder
.balanced_phenotype_labels <- function(
    n_samples,
    phenotypes
) {
  if (length(phenotypes) < 2L) {
    stop("At least two phenotypes are required.", call. = FALSE)
  }
  n_ph <- length(phenotypes)
  base_n <- n_samples %/% n_ph
  extra <- n_samples %% n_ph
  # floor(n/K) per level; distribute +1 to the first `extra` levels
  counts <- rep(base_n, n_ph)
  if (extra > 0L) {
    counts[seq_len(extra)] <- counts[seq_len(extra)] + 1L
  }
  labels <- rep(phenotypes, times = counts)
  as.character(forcats::fct_inorder(factor(labels, levels = phenotypes)))
}


#' Default number of cells per pseudo-bulk sample
#'
#' @description Returns the median barcode cell count in \code{meta}.
#'
#' @param meta Cell metadata \code{data.frame}.
#' @param barcode_col Column name for barcode identifiers.
#'
#' @return Integer scalar.
#'
#' @keywords internal
.default_cells_per_sample <- function(meta, barcode_col) {
  as.integer(stats::median(
    table(meta[[barcode_col]]),
    na.rm = TRUE
  ))
}


#' Extract gene feature metadata from a Seurat assay
#'
#' @description Returns \code{assay@meta.features}, or a minimal
#'   \code{S4Vectors::DataFrame} of gene names when feature metadata is empty.
#'
#' @param seurat_obj A \pkg{Seurat} object.
#' @param assay Assay name.
#'
#' @return A \code{DataFrame} with one row per gene.
#'
#' @keywords internal
.get_feature_metadata <- function(seurat_obj, assay = "RNA") {
  assay_obj <- seurat_obj[[assay]]
  feat <- tryCatch(assay_obj@meta.features, error = function(e) NULL)
  if (is.null(feat) || nrow(feat) == 0L) {
    gene_names <- rownames(assay_obj)
    if (is.null(gene_names)) {
      gene_names <- rownames(
        SeuratObject::LayerData(assay_obj, layer = "counts")
      )
    }
    feat <- S4Vectors::DataFrame(
      gene = gene_names,
      row.names = gene_names
    )
  }
  feat
}


#' Build a SummarizedExperiment from pseudo-bulk counts
#'
#' @description Assembles raw count matrix, sample metadata, and gene
#'   annotations into a \pkg{SummarizedExperiment} object.
#'
#' @param count_matrix Numeric matrix (genes x samples).
#' @param col_data Sample metadata \code{data.frame}.
#' @param feature_data Gene metadata \code{DataFrame}.
#'
#' @return A \pkg{SummarizedExperiment} with assay \code{"counts"}.
#'
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


#' Infer log-normal gene parameters from single cells
#'
#' @title Log-normal parameter estimation on \code{log1p} counts
#' @description
#' Estimates per-gene mean and standard deviation on the
#' \eqn{\log(1+x)} scale within a cell subset. Used for the generative
#' log-normal pseudo-bulk path (strategy 7).
#'
#' @details
#' For gene \eqn{g} and cells \eqn{c \in \mathcal{C}}, with subsampled
#' counts \eqn{x_{gc}}:
#' \deqn{\hat{\mu}_g = \frac{1}{|\mathcal{C}|}\sum_{c \in \mathcal{C}}
#'   \log(1 + x_{gc})}
#' \deqn{\hat{\sigma}_g = \mathrm{sd}_{c \in \mathcal{C}}
#'   \left(\log(1 + x_{gc})\right)}
#'
#' Simulated counts are drawn as
#' \eqn{x_{gc}^{\mathrm{sim}} = \max\{0, \lfloor e^{\tilde{z}_{gc}} - 1 \rfloor\}}
#' with \eqn{\tilde{z}_{gc} \sim \mathcal{N}(\hat{\mu}_g, \hat{\sigma}_g^2)}.
#'
#' This stabilises zeros before log transformation, following common
#' pseudo-bulk pipelines (strategy 7 in
#' \code{vignettes/pseudo_bulk_generation.qmd}) and mean-pseudo-bulk
#' approaches using \code{log2(counts + epsilon)}.
#'
#' @param counts Sparse or dense count matrix (genes x cells).
#' @param cell_idx Integer indices of cells in the estimation subset.
#' @param max_cells Maximum cells used for estimation (random subsample).
#'
#' @return List with \code{mean_log} and \code{sd_log} numeric vectors.
#'
#' @references
#' \url{https://oshlacklab.com/splatter/articles/splatter.html} (log-normal
#' library-size modelling motivation).
#'
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
  # z_gc = log(1 + x_gc); estimate mu_g, sigma_g per gene
  log_x <- log1p(as.matrix(counts[, idx, drop = FALSE]))
  list(
    mean_log = rowMeans(log_x),
    sd_log = matrixStats::rowSds(log_x)
  )
}


#' Method-of-moments negative-binomial size from counts
#'
#' @title Per-gene NB dispersion via method of moments
#' @description
#' Computes the negative-binomial \code{size} (dispersion) parameter from a
#' vector of counts using method-of-moments, as in HADACA3 pseudo-single-cell
#' generation.
#'
#' @details
#' For counts with mean \eqn{m} and variance \eqn{v}, under
#' \eqn{X \sim \mathrm{NB}(\mu=m, \mathrm{size}=\theta)}:
#' \deqn{\mathrm{Var}(X) = m + \frac{m^2}{\theta}
#'   \quad\Rightarrow\quad
#'   \hat{\theta} = \frac{m^2}{v - m}}
#' when \eqn{v > m}. Otherwise \eqn{\hat{\theta} = \infty} (Poisson limit).
#'
#' @param x Numeric vector of counts for one gene.
#'
#' @return Dispersion \code{size} (theta); may be \code{Inf}.
#'
#' @references
#' HADACA3 in silico pseudo-bulk workflow
#' (\url{https://github.com/bioinfo-LIG/hadaca3_framework}), which fits
#' \code{MASS::glm.nb} per gene; method-of-moments is used here for speed
#' at scale.
#'
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


#' Infer negative-binomial gene parameters from single cells
#'
#' @title Per-gene NB mean and dispersion estimation
#' @description
#' Estimates per-gene mean count and NB dispersion within a cell subset,
#' following the HADACA3 / Splatter count-modelling rationale.
#'
#' @details
#' For gene \eqn{g} and cells \eqn{c \in \mathcal{C}}:
#' \deqn{\hat{\mu}_g = \frac{1}{|\mathcal{C}|}\sum_{c \in \mathcal{C}} x_{gc}}
#' \deqn{\hat{\theta}_g = \frac{\hat{\mu}_g^2}{
#'   \widehat{\mathrm{Var}}(x_{g\cdot}) - \hat{\mu}_g}}
#'
#' Simulated counts use \code{stats::rnbinom} parameterised by
#' \eqn{(\mu_g, \mathrm{size}=\theta_g)}.
#'
#' @param counts Sparse or dense count matrix (genes x cells).
#' @param cell_idx Integer indices of cells in the estimation subset.
#' @param max_cells Maximum cells used for estimation (random subsample).
#'
#' @return List with \code{mu} and \code{size} numeric vectors.
#'
#' @references
#' Zappia et al. (2017) Splatter: simulation of single-cell RNA sequencing
#' data. \emph{Genome Biology} 18, 174.
#' \doi{10.1186/s13059-017-1303-1}
#'
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
  # theta_g = mu_g^2 / (Var_g - mu_g) per gene
  size <- apply(x, 1L, .mom_nb_size)
  size[!is.finite(size) | size <= 0] <- 1e6
  list(mu = mu, size = size)
}


#' Estimate Splatter simulation parameters from single cells
#'
#' @title Splatter \code{splatEstimate} on a phenotype cell subset
#' @description
#' Wraps \code{splatter::splatEstimate} to infer global scRNA-seq simulation
#' parameters, including library-size location and scale.
#'
#' @details
#' Splatter models total library size \eqn{L_c} per cell (often log-normal)
#' and gene-wise means as functions of \eqn{L_c}. \code{splatEstimate}
#' returns a \code{SplatParams} object with, among others, library-size
#' \eqn{(\mathrm{loc}_L, \mathrm{scale}_L)} used to draw technical depth
#' variation at simulation time.
#'
#' @param counts Sparse or dense count matrix (genes x cells).
#' @param cell_idx Integer indices of cells in the estimation subset.
#' @param max_cells Maximum cells passed to \code{splatEstimate}.
#'
#' @return A \code{splatter::SplatParams} object.
#'
#' @references
#' \url{https://oshlacklab.com/splatter/reference/splatEstimate.html};
#' Zappia et al. (2017) \doi{10.1186/s13059-017-1303-1}.
#'
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


#' Infer biological and technical simulation parameters
#'
#' @title Hierarchical parameter estimation (phenotype + barcode)
#' @description
#' Estimates generative-model parameters at the biological level
#' (phenotype) and, for the log-normal model only, at the technical level
#' (barcode library-size scaling).
#'
#' @details
#' **Biological level** (per phenotype \eqn{p}):
#' \itemize{
#'   \item Log-normal: \eqn{(\hat{\mu}_g, \hat{\sigma}_g)} on
#'     \eqn{\log(1+x)} scale.
#'   \item Negative binomial: per-gene \eqn{(\hat{\mu}_g, \hat{\theta}_g)}
#'     plus \code{splatter::splatEstimate} for library-size parameters.
#' }
#'
#' **Technical level** (log-normal only): per barcode \eqn{b},
#' \deqn{s_b = \frac{\sum_{c: b(c)=b} L_c}{\overline{\sum_{c} L_c}}}
#' where \eqn{L_c = \sum_g x_{gc}.} This factor is applied when simulating
#' log-normal cells anchored to barcode \eqn{b}.
#'
#' For the negative-binomial path, barcode library scaling is \emph{not}
#' applied separately: Splatter already estimates library-size variation
#' via \code{splatEstimate}, and an extra scaling factor would be redundant.
#'
#' @param counts Count matrix (genes x cells).
#' @param meta Cell metadata.
#' @param phenotype_col Phenotype column name.
#' @param barcode_col Barcode column name.
#' @param model \code{"lognormal"} or \code{"negative_binomial"}.
#'
#' @return List with \code{biological} and \code{technical} components.
#'
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

  # --- Biological level: per-phenotype gene parameters --------------------
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

  # --- Technical level: barcode library scaling (log-normal only) ---------
  technical <- NULL
  if (model == "lognormal") {
    barcodes <- sort(unique(as.character(meta[[barcode_col]])))
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
  }

  list(biological = biological, technical = technical)
}


#' Draw per-cell library-size factors from Splatter estimates
#'
#' @description
#' Samples cell-level library-depth multipliers from the log-normal library
#' model inferred by \code{splatter::splatEstimate}.
#'
#' @details
#' With Splatter library parameters \eqn{(\mathrm{loc}_L, \mathrm{scale}_L)}:
#' \deqn{\log L_c \sim \mathcal{N}(\mathrm{loc}_L, \mathrm{scale}_L^2)}
#' Factors are normalised to mean 1 within the simulated cell batch.
#'
#' @param n_cells Number of cells to simulate.
#' @param splatter_params A \code{SplatParams} object.
#'
#' @return Numeric vector of length \code{n_cells}.
#'
#' @keywords internal
.draw_splatter_library_factors <- function(n_cells, splatter_params) {
  loc <- splatter_params@lib.loc
  scale <- splatter_params@lib.scale
  lib_norm <- splatter_params@lib.norm
  # splatEstimate sets lib.norm = TRUE when library sizes look Gaussian on
  # the count scale (see splatter NOTE after splatEstimate).
  if (isTRUE(lib_norm) && loc > 1e4) {
    sizes <- stats::rnorm(n_cells, mean = loc, sd = scale)
  } else if (isTRUE(lib_norm)) {
    sizes <- stats::rlnorm(n_cells, meanlog = loc, sdlog = scale)
  } else {
    sizes <- stats::rnorm(n_cells, mean = loc, sd = scale)
  }
  sizes <- pmax(sizes, 1)
  sizes / mean(sizes)
}


#' Simulate log-normal single-cell counts
#'
#' @description
#' Draws \code{n_cells} count vectors from per-gene log-normal models,
#' optionally scaled by a barcode library factor.
#'
#' @details
#' For each cell \eqn{c} and gene \eqn{g}:
#' \deqn{z_{gc} \sim \mathcal{N}(\hat{\mu}_g, \hat{\sigma}_g^2), \quad
#'   x_{gc} = \max\left\{0, \left\lfloor e^{z_{gc}} \cdot s_b - 1 \right\rfloor
#'   \right\}}
#'
#' @param bio_params List from \code{.infer_lognormal_gene_params()}.
#' @param tech_scale Barcode library-size scaling factor.
#' @param n_cells Number of cells to simulate.
#' @param seed Optional random seed.
#'
#' @return Numeric matrix (genes x cells).
#'
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
  # z_gc ~ N(mu_g, sigma_g); x_gc = max(0, round(expm1(z) * tech_scale))
  z <- stats::rnorm(
    n = n_genes * n_cells,
    mean = rep(mu, each = n_cells),
    sd = rep(sd, each = n_cells)
  )
  mat <- matrix(z, nrow = n_genes, ncol = n_cells)
  out <- pmax(0, round(expm1(mat) * tech_scale))
  matrix(out, nrow = n_genes, ncol = n_cells)
}


#' Simulate negative-binomial single-cell counts
#'
#' @description
#' Draws \code{n_cells} count vectors from per-gene NB models with
#' cell-specific library depth from Splatter estimates.
#'
#' @details
#' For cell \eqn{c} with library factor \eqn{\ell_c} (from Splatter) and
#' gene \eqn{g}:
#' \deqn{x_{gc} \sim \mathrm{NB}\!\left(\mu = \hat{\mu}_g \cdot \ell_c,\;
#'   \mathrm{size} = \hat{\theta}_g\right)}
#'
#' @param bio_params List with \code{gene} and \code{splatter} components.
#' @param n_cells Number of cells to simulate.
#' @param seed Optional random seed.
#'
#' @return Numeric matrix (genes x cells).
#'
#' @keywords internal
.simulate_nb_cells <- function(
    bio_params,
    n_cells,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  mu_g <- bio_params$gene$mu
  size_g <- bio_params$gene$size
  n_genes <- length(mu_g)
  # l_c ~ Splatter library model; x_gc ~ NB(mu_g * l_c, size_g)
  lib_factors <- .draw_splatter_library_factors(
    n_cells,
    bio_params$splatter
  )
  mat <- matrix(0, nrow = n_genes, ncol = n_cells)
  for (j in seq_len(n_cells)) {
    mat[, j] <- stats::rnbinom(
      n = n_genes,
      size = size_g,
      mu = mu_g * lib_factors[[j]]
    )
  }
  mat
}


#' Pick a random anchor barcode for a phenotype
#'
#' @description
#' Samples one barcode identifier among organoids carrying the target
#' phenotype. Used to tag simulated samples and, for log-normal simulation,
#' to select a technical library-scaling factor.
#'
#' @param meta Cell metadata.
#' @param phenotype_col Phenotype column name.
#' @param barcode_col Barcode column name.
#' @param phenotype Target phenotype level.
#'
#' @return Character scalar (barcode id).
#'
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
#' Strategy 1 from \code{vignettes/pseudo_bulk_generation.qmd}: treat each
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
#' @param cell_mask Optional logical vector aligned with \code{seurat_obj}
#'   cells; when provided, only these cells are aggregated.
#'
#' @return A \pkg{SummarizedExperiment} with one column per distinct barcode.
#'   \code{colData} contains \code{barcode_id}, the phenotype label,
#'   \code{library_depth} (column sum of counts), and
#'   \code{simulation_method = "barcode_aggregation"}.
#'
#' @seealso `simulate_bootstrap_samples()`, `simulate_generative_models()`,
#'   \code{vignettes/pseudo_bulk_generation.qmd}
#'
#' @importFrom Seurat GetAssay
#' @importFrom glue glue
#' @export
aggregate_barcode_pseudo_bulk <- function(
    seurat_obj,
    phenotype_col = "Morphotype",
    barcode_col = "Sample.barcode",
    assay = "RNA",
    layer = "counts",
    cell_mask = NULL
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop(glue::glue("`seurat_obj` must be a Seurat object."), call. = FALSE)
  }
  meta <- seurat_obj@meta.data
  missing <- setdiff(c(phenotype_col, barcode_col), colnames(meta))
  if (length(missing) > 0L) {
    stop(
      "Missing metadata column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop(glue::glue("Assay `{assay}` not found in Seurat object."), call. = FALSE)
  }
  counts <- SeuratObject::LayerData(
    object = seurat_obj,
    assay = assay,
    layer = layer
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  if (!identical(colnames(counts), rownames(meta))) {
    missing_cells <- setdiff(colnames(counts), rownames(meta))
    if (length(missing_cells) > 0L) {
      stop(
        "Assay cell names must be present in metadata row names.",
        call. = FALSE
      )
    }
    meta <- meta[colnames(counts), , drop = FALSE]
  }
  if (!is.null(cell_mask)) {
    if (length(cell_mask) != nrow(meta)) {
      stop("`cell_mask` length must match number of cells.", call. = FALSE)
    }
    meta <- meta[cell_mask, , drop = FALSE]
    counts <- counts[, cell_mask, drop = FALSE]
  }
  if (anyNA(meta[[phenotype_col]]) || anyNA(meta[[barcode_col]])) {
    stop(
      "Phenotype and barcode columns must not contain NA values.",
      call. = FALSE
    )
  }

  barcodes <- sort(unique(as.character(meta[[barcode_col]])))
  n_bar <- length(barcodes)
  bulk_mat <- matrix(
    0,
    nrow = nrow(counts),
    ncol = n_bar,
    dimnames = list(rownames(counts), barcodes)
  )

  barcode_id <- character(n_bar)
  phenotype <- character(n_bar)
  library_depth <- numeric(n_bar)

  for (i in seq_along(barcodes)) {
    bc <- barcodes[[i]]
    idx <- which(meta[[barcode_col]] == bc)
    # Y_gb = sum_{c: b(c)=b} x_gc
    bulk_mat[, i] <- Matrix::rowSums(counts[, idx, drop = FALSE])
    barcode_id[[i]] <- bc
    phenotype[[i]] <- as.character(meta[[phenotype_col]][idx[[1L]]])
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


#' Aggregate pseudo-bulk samples by cell type and barcode
#'
#' @description
#' Extension of `aggregate_barcode_pseudo_bulk()`: sum raw single-cell counts
#' within each cell type x barcode combination, yielding one pseudo-bulk
#' column per organoid within each annotated cell type. Uses the same raw-count
#' aggregation principle as Strategy 1 in
#' \code{vignettes/pseudo_bulk_generation.qmd}.
#'
#' @param seurat_obj A \pkg{Seurat} object with raw counts in the RNA assay.
#' @param celltype_col Metadata column for cell-type labels
#'   (e.g. \code{"luque_cluster_annotation"}).
#' @param phenotype_col Metadata column for organoid-level phenotype labels
#'   (e.g. \code{"Morphotype"}).
#' @param barcode_col Metadata column for barcode / organoid identifiers
#'   (e.g. \code{"Sample.barcode"}).
#' @param assay Assay name from which raw counts are read. Default
#'   \code{"RNA"}.
#' @param cell_mask Optional logical vector aligned with \code{seurat_obj}
#'   cells; when provided, only these cells are aggregated.
#'
#' @return A \pkg{SummarizedExperiment} with one column per distinct
#'   cell type x barcode pair. \code{colData} contains \code{celltype},
#'   \code{barcode_id}, the phenotype label, \code{n_cells},
#'   \code{library_depth}, and
#'   \code{simulation_method = "celltype_barcode_aggregation"}.
#'
#' @seealso `aggregate_barcode_pseudo_bulk()`, `simulate_bootstrap_samples()`,
#'   `simulate_generative_models()`, \code{vignettes/pseudo_bulk_generation.qmd}
#'
#' @importFrom Seurat GetAssay
#' @importFrom glue glue
#' @export
aggregate_celltype_barcode_pseudo_bulk <- function(
    seurat_obj,
    celltype_col = "luque_cluster_annotation",
    phenotype_col = "Morphotype",
    barcode_col = "Sample.barcode",
    assay = "RNA",
    cell_mask = NULL
) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop(glue::glue("`seurat_obj` must be a Seurat object."), call. = FALSE)
  }
  meta <- seurat_obj@meta.data
  missing <- setdiff(
    c(celltype_col, phenotype_col, barcode_col),
    colnames(meta)
  )
  if (length(missing) > 0L) {
    stop(
      "Missing metadata column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop(glue::glue("Assay `{assay}` not found in Seurat object."), call. = FALSE)
  }
  counts <- SeuratObject::LayerData(
    object = seurat_obj,
    assay = assay,
    layer = "counts"
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }
  if (!is.null(cell_mask)) {
    if (length(cell_mask) != nrow(meta)) {
      stop("`cell_mask` length must match number of cells.", call. = FALSE)
    }
    meta <- meta[cell_mask, , drop = FALSE]
    counts <- counts[, cell_mask, drop = FALSE]
  }
  if (anyNA(meta[[celltype_col]]) ||
      anyNA(meta[[phenotype_col]]) ||
      anyNA(meta[[barcode_col]])) {
    stop(
      "Cell type, phenotype and barcode columns must not contain NA values.",
      call. = FALSE
    )
  }

  meta$._pb_group <- interaction(
    meta[[celltype_col]],
    meta[[barcode_col]],
    drop = TRUE,
    sep = "|"
  )
  pb_groups <- sort(unique(as.character(meta$._pb_group)))
  n_pb <- length(pb_groups)
  bulk_mat <- matrix(
    0,
    nrow = nrow(counts),
    ncol = n_pb,
    dimnames = list(rownames(counts), pb_groups)
  )

  celltype <- character(n_pb)
  barcode_id <- character(n_pb)
  phenotype <- character(n_pb)
  n_cells <- integer(n_pb)
  library_depth <- numeric(n_pb)

  for (i in seq_along(pb_groups)) {
    grp <- pb_groups[[i]]
    idx <- which(meta$._pb_group == grp)
  # Y_gb = sum_{c: group(c)=g} x_gc
    bulk_mat[, i] <- Matrix::rowSums(counts[, idx, drop = FALSE])
    celltype[[i]] <- as.character(meta[[celltype_col]][idx[[1L]]])
    barcode_id[[i]] <- as.character(meta[[barcode_col]][idx[[1L]]])
    phenotype[[i]] <- as.character(meta[[phenotype_col]][idx[[1L]]])
    n_cells[[i]] <- length(idx)
    library_depth[[i]] <- sum(bulk_mat[, i])
  }

  col_data <- data.frame(
    celltype = celltype,
    barcode_id = barcode_id,
    phenotype = phenotype,
    n_cells = n_cells,
    library_depth = library_depth,
    simulation_method = "celltype_barcode_aggregation",
    row.names = pb_groups,
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
#' Strategies 2 and 3 from \code{vignettes/pseudo_bulk_generation.qmd}.
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
#' @param n_samples Total number of pseudo-bulk samples to generate.
#'   Phenotypes are allocated as evenly as possible across samples.
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
#' @seealso `aggregate_barcode_pseudo_bulk()`, `simulate_generative_models()`
#'
#' @importFrom Seurat GetAssay
#' @importFrom glue glue
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

  if (!inherits(seurat_obj, "Seurat")) {
    stop(glue::glue("`seurat_obj` must be a Seurat object."), call. = FALSE)
  }
  meta <- seurat_obj@meta.data
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

  if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop(glue::glue("Assay `{assay}` not found in Seurat object."), call. = FALSE)
  }
  counts <- SeuratObject::LayerData(
    object = seurat_obj,
    assay = assay,
    layer = "counts"
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }

  phenotypes <- sort(unique(as.character(meta[[phenotype_col]])))
  pheno_labels <- .balanced_phenotype_labels(n_samples, phenotypes)

  if (is.null(cells_per_sample)) {
    cells_per_sample <- .default_cells_per_sample(meta, barcode_col)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  bulk_mat <- matrix(
    0,
    nrow = nrow(counts),
    ncol = n_samples,
    dimnames = list(
      rownames(counts),
      paste0("bootstrap_", seq_len(n_samples))
    )
  )

  barcode_id <- character(n_samples)
  library_depth <- numeric(n_samples)

  for (i in seq_len(n_samples)) {
    ph <- pheno_labels[[i]]

    # --- Sample cells C_i from the appropriate pool -----------------------
    if (replicate_type == "biological") {
      # C_i ~ Sample{c : y_{b(c)} = y_i}, |C_i| = N
      pool <- which(meta[[phenotype_col]] == ph)
      barcode_id[[i]] <- "pooled"
    } else {
      bc <- .pick_anchor_barcode(meta, phenotype_col, barcode_col, ph)
      # C_i ~ Sample{c : b(c)=b}, |C_i| = N
      pool <- which(
        meta[[barcode_col]] == bc &
          meta[[phenotype_col]] == ph
      )
      barcode_id[[i]] <- bc
    }

    sampled <- sample(pool, size = cells_per_sample, replace = TRUE)
    # Y_gi = sum_{c in C_i} x_gc
    bulk_mat[, i] <- Matrix::rowSums(counts[, sampled, drop = FALSE])
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
#' in \code{vignettes/pseudo_bulk_generation.qmd}). Parameters are inferred at
#' the biological (phenotype) level; technical library scaling is applied
#' only for the log-normal model. Simulated cells are summed to sample-level
#' raw bulk counts.
#'
#' @param seurat_obj A \pkg{Seurat} object with raw counts.
#' @param phenotype_col Metadata column for phenotype labels.
#' @param barcode_col Metadata column for barcode identifiers.
#' @param n_samples Total pseudo-bulk samples (balanced across phenotypes).
#' @param model \code{"lognormal"} for log-normal gene counts, or
#'   \code{"negative_binomial"} for Splatter-style negative-binomial counts.
#' @param cells_per_sample Cells simulated and summed per pseudo-bulk sample.
#' @param assay Assay from which raw counts are read.
#' @param seed Optional random seed.
#'
#' @return A \pkg{SummarizedExperiment} with raw summed counts.
#'
#' @details
#' **Log-normal path**: per-gene \eqn{(\hat{\mu}_g, \hat{\sigma}_g)} on
#' \eqn{\log(1+x)} scale per phenotype, with barcode library scaling
#' \eqn{s_b} at simulation time.
#'
#' **Negative-binomial path**: per-gene \eqn{(\hat{\mu}_g, \hat{\theta}_g)}
#' via method-of-moments (HADACA3-style) and per-cell library factors from
#' \code{splatter::splatEstimate}. No extra barcode scaling is applied,
#' because Splatter already models library-size variation.
#'
#' Bulk aggregation: \deqn{Y_{gi} = \sum_{c \in \mathcal{C}_i} x_{gc}^{\mathrm{sim}}}
#'
#' @seealso `aggregate_barcode_pseudo_bulk()`, `simulate_bootstrap_samples()`,
#'   \url{https://oshlacklab.com/splatter/articles/splatter.html}
#'
#' @importFrom Seurat GetAssay
#' @importFrom glue glue
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

  if (!inherits(seurat_obj, "Seurat")) {
    stop(glue::glue("`seurat_obj` must be a Seurat object."), call. = FALSE)
  }
  meta <- seurat_obj@meta.data
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

  if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop(glue::glue("Assay `{assay}` not found in Seurat object."), call. = FALSE)
  }
  counts <- SeuratObject::LayerData(
    object = seurat_obj,
    assay = assay,
    layer = "counts"
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- Matrix::Matrix(counts, sparse = TRUE)
  }

  phenotypes <- sort(unique(as.character(meta[[phenotype_col]])))
  pheno_labels <- .balanced_phenotype_labels(n_samples, phenotypes)

  if (is.null(cells_per_sample)) {
    cells_per_sample <- .default_cells_per_sample(meta, barcode_col)
  }

  # --- Parameter inference: biological (+ technical for log-normal) -------
  params <- .infer_hierarchical_params(
    counts = counts,
    meta = meta,
    phenotype_col = phenotype_col,
    barcode_col = barcode_col,
    model = model
  )

  bulk_mat <- matrix(
    0,
    nrow = nrow(counts),
    ncol = n_samples,
    dimnames = list(
      rownames(counts),
      paste0(model, "_", seq_len(n_samples))
    )
  )

  barcode_id <- character(n_samples)
  library_depth <- numeric(n_samples)

  for (i in seq_len(n_samples)) {
    ph <- pheno_labels[[i]]
    bc <- .pick_anchor_barcode(meta, phenotype_col, barcode_col, ph)
    bio <- params$biological[[ph]]
    cell_seed <- if (is.null(seed)) NULL else seed + i

    # --- Simulate single cells, then aggregate to bulk --------------------
    cell_mat <- if (model == "lognormal") {
      tech_scale <- params$technical$library_scale[[bc]]
      .simulate_lognormal_cells(
        bio_params = bio,
        tech_scale = tech_scale,
        n_cells = cells_per_sample,
        seed = cell_seed
      )
    } else {
      # NB: library depth from Splatter only (no extra barcode scaling)
      .simulate_nb_cells(
        bio_params = bio,
        n_cells = cells_per_sample,
        seed = cell_seed
      )
    }

    # Y_gi = sum_{c in C_i} x_gc^sim
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
