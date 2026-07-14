# ============================================================================
# 1. Parameters: screening libraries, file settings, and loading ----
# ============================================================================
#
# Critical framing (see also the review at the end of the corresponding
# chat/vignette): there is no independently collected bulk cohort for this
# study. `naive_pseudo_bulk_raw_counts.rds` (barcode-level sums of the same
# 120h single cells; see `scripts/05_01_generate_pseudo_bulk.R`) is used here
# as the "bulk" reference required by SigBridgeR's phenotype-guided screens.

# nohup Rscript --no-save --no-restore scripts/06_04_SigBridgeR_bulk_guided_single_cell_analysis.R \
#   > "logs/luque_$(date +%F)_SigBridgeR-analysis.log" 2>&1 &

## 1.1 Screening algorithm libraries and file settings ----
#
# Dependencies are called with `package::function()` to keep startup quiet in
# batch logs (no attach messages from Bioconductor/tidyverse stacks).

study <- "luque"
output_dir <- "outputs/phenotype-aware"
source("./R/utils.R")
source("./R/simulate_pseudo_bulk_samples.R")
# source("./R/degas_uv_environment.R") # disabled with DEGAS

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

PATHS <- list(
  sc_merged = "./data/intermediate/luque_single_cell_merged_2026-06-07.rds",
  bulk_pseudo_bulk = "./data/intermediate/luque_naive_pseudo_bulk_raw_counts.rds"
  # degas_project = file.path(".uv", "degas") # disabled with DEGAS
)

COL <- list(
  cell_annotation = "luque_cluster_annotation", # pre-existing single-cell labels
  bulk_barcode = "Sample.barcode",
  bulk_morphotype = "Morphotype",
  sc_input_assay = "RNA",
  sc_analysis_assay = "SCT" # already SCT-filtered/normalised; no re-preprocessing
)

CFG <- list(
  phenotype_name = "NeuralBiasedMorphotype",
  # Binary encoding requested by the user: TLS is the reference class (0),
  # neural_bias is the positive class (1). Positive screen calls are thus
  # interpreted as "associated with neural-biased differentiation".
  reference_label = "TLS",
  positive_label = "neural_bias",
  reference_code = 0L,
  positive_code = 1L,
  phenotype_class = "binary",

  seed = 42L,
  cv_folds = 5L,
  bootstrap_n = 100L,
  fdr = 0.05,
  log2fc = 0.5,
  penalty_mix = 0.5,
  selection_fraction = 0.20,
  min_common_genes = 500L,
  min_samples_per_class = 5L,
  # Multi-core budget (disabled — caused segfaults / pthread errors on this data):
  # cpu_budget = min(
  #   8L,
  #   max(1L, parallelly::availableCores() %/% 2L)
  # ),
  cpu_budget = 1L

  # degas_python_version = "3.8.18", # disabled with DEGAS
  # degas_python_packages = c(
  #   "tensorflow==2.4.1",
  #   "protobuf==3.20.3",
  #   "numpy==1.19.5"
  # )
)
# Session-wide single-thread caps (OpenMP / BLAS backends).
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)
# SigBridgeR package options (seed, verbosity); see Other_Function_Details vignette.
SigBridgeR::setFuncOption(seed = CFG$seed)

# Backend thread budget before any TensorFlow import (DEGAS).
# Multi-core setup (disabled):
# SigBridgeR::setThreads(
#   threads = CFG$cpu_budget,
#   cheapr = NULL,
#   dt = 1L,
#   openmp = 1L,
#   tf_config = list(
#     inter_op = max(1L, CFG$cpu_budget %/% 4L),
#     intra_op = CFG$cpu_budget
#   ),
#   verbose = TRUE
# )
SigBridgeR::setThreads(
  threads = 1L,
  cheapr = NULL,
  dt = 1L,
  openmp = 1L,
  tf_config = list(inter_op = 1L, intra_op = 1L),
  verbose = TRUE
)

SCREEN <- list(
  shared = list(
    assay = COL$sc_analysis_assay, # single-cell assay/layer scored by every method
    seed = CFG$seed, # common RNG seed passed to each backend
    cv_folds = CFG$cv_folds, # cross-validation folds (LP-SGL; Scissor reliability test)
    fdr = CFG$fdr, # FDR/DE significance threshold shared across DE-based steps
    log2fc = CFG$log2fc, # minimum bulk log2 fold-change to call a DE gene
    penalty_mix = CFG$penalty_mix, # L1 vs network-regularisation weight (`alpha`)
    selection_fraction = CFG$selection_fraction # target top/bottom fraction of scored cells
  ),
  scPP = list(
    ref_group = CFG$reference_code, # bulk class treated as the AUCell baseline (TLS)
    Log2FC_cutoff = CFG$log2fc, # bulk DE threshold defining the marker gene set
    estimate_cutoff = 0.2, # minimum AUCell enrichment estimate to keep a marker
    probs = CFG$selection_fraction, # AUCell score quantile defining Positive/Negative cells
    # parallel = TRUE # scPP internal AUCell parallel ranking (disabled)
    parallel = FALSE
  ),
  LP_SGL = list(
    family = "logit", # binomial link for the sparse-group-lasso regression
    resolution = 0.6, # Leiden resolution defining the cell communities/groups
    alpha = CFG$penalty_mix, # group-lasso vs L1 sparsity balance
    nfold = CFG$cv_folds, # cross-validation folds for the lasso path
    dge_analysis = list(
      run = TRUE, # derive bulk DE genes before fitting LP-SGL
      logFC_threshold = CFG$log2fc, # DE log2FC cutoff feeding the DE gene set
      pval_threshold = CFG$fdr # DE p-value cutoff feeding the DE gene set
    )
  ),
  # Scissor = list(
  #   family = "binomial",
  #   alpha = CFG$penalty_mix,
  #   cutoff = CFG$selection_fraction,
  #   reliability_test = list(
  #     run = TRUE,
  #     n = CFG$bootstrap_n,
  #     nfold = CFG$cv_folds
  #   ),
  #   intermediate_dir = file.path(output_dir, "Scissor_luque_res")
  # ),
  # DEGAS = list(
  #   degas_params = list(
  #     DEGAS.model_type = "ClassClass",
  #     DEGAS.architecture = "Standard",
  #     DEGAS.ff_depth = 3L,
  #     DEGAS.bag_depth = CFG$bootstrap_n,
  #     DEGAS.train_steps = 2000L,
  #     DEGAS.scbatch_sz = 200L,
  #     DEGAS.patbatch_sz = 50L,
  #     DEGAS.hidden_feats = 50L,
  #     DEGAS.lambda1 = CFG$penalty_mix,
  #     DEGAS.lambda2 = CFG$penalty_mix,
  #     DEGAS.lambda3 = CFG$penalty_mix,
  #     DEGAS.seed = CFG$seed
  #   ),
  #   sc_data.pheno_colname = COL$cell_annotation,
  #   select_fraction = CFG$selection_fraction,
  #   intermediate_dir = file.path(output_dir, "DEGAS_res")
  # ),
  consensus = list(
    method_weights = c(
      scPP = 1 / 2,
      LP_SGL = 1 / 2
      # DEGAS = 1 / 3 # disabled
      # Scissor = 0.25 # disabled
    ),
    ties.method = "random" # tie-breaking rule when the weighted vote is exactly split
  )
)

METHOD_COLUMNS <- c(
  scPP = "scPP",
  LP_SGL = "LP_SGL"
  # DEGAS = "DEGAS" # disabled
  # Scissor = "scissor" # disabled
)

## 1.2 Loading ----

luque_merged <- readRDS(PATHS$sc_merged)
if (
  !identical(
    colnames(luque_merged),
    colnames(SeuratObject::GetAssayData(
      luque_merged,
      assay = "RNA",
      layer = "counts"
    ))
  )
) {
  luque_merged <- realign_luque_seurat(luque_merged)
}

# 120h only, and only cells whose organoid barcode already carries a
# Morphotype call and an existing cluster annotation (see Section 2/3 for why
# this matters: every retained cell already has a known ground-truth label).
sc_120h <- luque_merged[
  ,
  Seurat::WhichCells(
    luque_merged,
    expression = timepoint == "120h" &
      !is.na(Morphotype) &
      !is.na(luque_cluster_annotation)
  )
]

naive_pseudo_bulk_raw_counts <- readRDS(PATHS$bulk_pseudo_bulk)

message(
  glue::glue(
    "Loaded {ncol(sc_120h)} single cells (120h) and ",
    "{ncol(naive_pseudo_bulk_raw_counts)} pseudo-bulk samples."
  )
)


# ============================================================================
# 2. Preprocessing: bulk, single-cell, and phenotype ----
# ============================================================================

## 2.1 Binary phenotype (reference = TLS -> 0, positive = neural_bias -> 1) ----

morphotype_raw <- SummarizedExperiment::colData(naive_pseudo_bulk_raw_counts) |>
  as.data.frame() |>
  dplyr::pull(.data[[COL$bulk_morphotype]]) |>
  as.character() |>
  stats::setNames(colnames(naive_pseudo_bulk_raw_counts))

phenotype_mapped <- SigBridgeR::PhenoMap(
  morphotype_raw,
  morphotype_raw == CFG$positive_label ~ CFG$positive_code,
  morphotype_raw == CFG$reference_label ~ CFG$reference_code
) |>
  as.integer() |>
  stats::setNames(names(morphotype_raw))

stopifnot(
  length(unique(phenotype_mapped)) == 2L,
  all(table(phenotype_mapped) >= CFG$min_samples_per_class)
)
table(Morphotype = morphotype_raw, Encoded = phenotype_mapped)

## 2.2 Bulk preprocessing (pseudo-bulk counts -> filtered log2-CPM) ----

# `BulkPreProcess()` performs sample/gene QC on raw counts; expression is
# then converted to log2-CPM (no per-gene z-scoring) for screening, matching
# SigBridgeR's own examples.
bulk_counts_filtered <- SigBridgeR::BulkPreProcess(
  data = list(
    count_matrix = as.matrix(SummarizedExperiment::assay(
      naive_pseudo_bulk_raw_counts,
      "counts"
    )),
    sample_info = data.frame(
      sample = names(phenotype_mapped),
      condition = phenotype_mapped
    )
  ),
  gene_symbol_conversion = FALSE, # mouse gene symbols already present
  min_count_threshold = 10L,
  min_gene_expressed = 3L,
  min_correlation = -1, # keep the QC report; do not drop samples by correlation
  seed = CFG$seed,
  verbose = TRUE
)

matched_bulk <- edgeR::cpm(
  bulk_counts_filtered,
  log = TRUE,
  prior.count = 1,
  normalized.lib.sizes = TRUE
)

## 2.3 Single-cell preprocessing ----

# No `SCTransform()`/re-clustering: the object already carries a filtered SCT
# assay and biologist-reviewed `luque_cluster_annotation` labels (see
# `scripts/01_02_add_standardised_cell_ontologies.R`); they are reused as-is.
SeuratObject::DefaultAssay(sc_120h) <- COL$sc_analysis_assay

# scPP ranks cells from the RNA `data` layer (AUCell); after `realign_luque_seurat()`
# RNA often has counts only. Other methods use the SCT assay from `screen_common`.
if (!"data" %in% SeuratObject::Layers(sc_120h[[COL$sc_input_assay]])) {
  sc_120h <- Seurat::NormalizeData(
    sc_120h,
    assay = COL$sc_input_assay,
    verbose = FALSE
  )
}

# LP-SGL expects a `{assay}_snn` graph (here `SCT_snn`). Built-in `RunPCA()` is
# numerically unstable on this object (inflated PC scores), so PCA uses
# `RSpectra::svds()` on scaled SCT `data` before neighbours and UMAP.
add_stable_pca <- function(
    obj,
    assay,
    nfeatures = 2000L,
    npcs = 30L
) {
  if (length(SeuratObject::VariableFeatures(obj)) == 0L) {
    obj <- Seurat::FindVariableFeatures(
      obj,
      assay = assay,
      nfeatures = nfeatures,
      verbose = FALSE
    )
  }
  vf <- SeuratObject::VariableFeatures(obj)
  mat <- SeuratObject::LayerData(obj, assay = assay, layer = "data")[
    vf,
    ,
    drop = FALSE
  ]
  if (inherits(mat, "dgCMatrix")) {
    mat <- as.matrix(mat)
  }
  vars <- matrixStats::rowVars(mat)
  keep <- is.finite(vars) & vars > 0
  mat <- mat[keep, , drop = FALSE]
  mat <- t(scale(t(mat)))
  mat[!is.finite(mat)] <- 0
  cell_matrix <- scale(t(mat), center = TRUE, scale = FALSE)
  pc <- RSpectra::svds(cell_matrix, k = npcs)
  emb <- pc$u %*% diag(pc$d, ncol = npcs)
  rownames(emb) <- colnames(obj)
  colnames(emb) <- paste0("PC_", seq_len(npcs))
  obj[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb,
    stdev = pc$d / sqrt(max(1, nrow(cell_matrix) - 1L)),
    assay = assay,
    key = "PC_"
  )
  obj
}

snn_graph <- paste0(COL$sc_analysis_assay, "_snn")
if (!snn_graph %in% names(sc_120h@graphs)) {
  sc_120h <- add_stable_pca(
    sc_120h,
    assay = COL$sc_analysis_assay,
    npcs = 30L
  )
  sc_120h <- Seurat::FindNeighbors(
    sc_120h,
    reduction = "pca",
    dims = 1:30,
    verbose = FALSE
  )
}

# UMAP panels in Section 4 reuse the PCA neighbourhood graph built above.
if (!"umap" %in% SeuratObject::Reductions(sc_120h)) {
  sc_120h <- Seurat::RunUMAP(
    sc_120h,
    dims = 1:30,
    reduction = "pca",
    verbose = FALSE
  )
}

## 2.4 Final alignment and validation ----

phenotype <- SigBridgeR::PhenoPreProcess(
  bulk = matched_bulk,
  phenotype = phenotype_mapped,
  phenotype_class = CFG$phenotype_class,
  verbose = TRUE
)
matched_bulk <- matched_bulk[, names(phenotype), drop = FALSE]

sc_expression <- SeuratObject::LayerData(
  sc_120h,
  assay = COL$sc_analysis_assay,
  layer = "data"
)
common_genes <- intersect(rownames(matched_bulk), rownames(sc_expression))
matched_bulk <- matched_bulk[common_genes, , drop = FALSE]
sc_120h <- sc_120h[common_genes, ]

input_validation <- tibble::tibble(
  item = c(
    "Single cells",
    "Bulk (pseudo-bulk) samples",
    "Shared genes",
    "Positive samples",
    "Reference samples"
  ),
  value = c(
    ncol(sc_120h),
    ncol(matched_bulk),
    length(common_genes),
    sum(phenotype == CFG$positive_code),
    sum(phenotype == CFG$reference_code)
  )
)
tinytable::tt(input_validation, caption = "Validated SigBridgeR inputs")

message(
  glue::glue(
    "Preprocessing complete: {ncol(matched_bulk)} bulk samples, ",
    "{ncol(sc_120h)} single cells, {length(common_genes)} shared genes."
  )
)


# ============================================================================
# 3. Screening algorithms and strategies ----
# ============================================================================
#
# Two built-in SigBridgeR screen methods active here (scPP, LP-SGL).
# Scissor and DEGAS are disabled.

# ## 3.1 Dedicated DEGAS Python environment (uv, not Conda) ----
#
# degas_python <- create_degas_uv_env(
#   project_dir = PATHS$degas_project,
#   python_version = CFG$degas_python_version,
#   packages = CFG$degas_python_packages,
#   verbose = TRUE
# )
# dir.create(
#   SCREEN$DEGAS$intermediate_dir,
#   recursive = TRUE,
#   showWarnings = FALSE
# )

## 3.2 Run the built-in screens sequentially ----
#
# SigBridgeR already registers scPP, LP-SGL, Scissor, and DEGAS as built-in
# `screen_method` values; hyper-parameters are passed directly to `Screen()`
# without custom executor wrappers or `Register()`.

# --- Scissor helpers (disabled) ---
# preprocessCore::normalize.quantiles() spawns pthread workers per column and
# fails (pthread_create 22) on wide matrices (genes × bulk + cells).
# safe_normalize_quantiles <- function(x, copy = TRUE, keep.names = FALSE) { ... }
# with_safe_scissor_quantile_norm <- function(expr) { ... }

# --- DEGAS helpers (disabled) ---
# make_degas_train_fn <- function() { ... }

screen_common <- list(
  matched_bulk = matched_bulk,
  sc_data = sc_120h,
  phenotype = phenotype,
  label_type = CFG$phenotype_name,
  phenotype_class = CFG$phenotype_class,
  assay = SCREEN$shared$assay,
  verbose = TRUE
)

screen_method_args <- list(
  scPP = list(
    screen_method = "scPP",
    assay = COL$sc_input_assay, # ScPP reads RNA `data`; SCT is for the other methods
    ref_group = SCREEN$scPP$ref_group,
    Log2FC_cutoff = SCREEN$scPP$Log2FC_cutoff,
    estimate_cutoff = SCREEN$scPP$estimate_cutoff,
    probs = SCREEN$scPP$probs,
    parallel = SCREEN$scPP$parallel # FALSE: single-thread AUCell
  ),
  LP_SGL = list(
    screen_method = "LP_SGL",
    family = SCREEN$LP_SGL$family,
    resolution = SCREEN$LP_SGL$resolution,
    alpha = SCREEN$LP_SGL$alpha,
    nfold = SCREEN$LP_SGL$nfold,
    dge_analysis = SCREEN$LP_SGL$dge_analysis
  )
  # Scissor = list(
  #   screen_method = "Scissor",
  #   family = SCREEN$Scissor$family,
  #   alpha = SCREEN$Scissor$alpha,
  #   cutoff = SCREEN$Scissor$cutoff,
  #   reliability_test = SCREEN$Scissor$reliability_test
  # ),
  # DEGAS = list(
  #   screen_method = "DEGAS",
  #   sc_data.pheno_colname = SCREEN$DEGAS$sc_data.pheno_colname,
  #   select_fraction = SCREEN$DEGAS$select_fraction,
  #   tmp_dir = SCREEN$DEGAS$intermediate_dir,
  #   save_cache = SCREEN$DEGAS$intermediate_dir,
  #   degas_params = c(
  #     SCREEN$DEGAS$degas_params,
  #     list(DEGAS.pyloc = degas_python)
  #   )
  # )
)

screen_checkpoint_path <- function(method) {
  file.path(output_dir, glue::glue("screen_{method}.rds"))
}

screen_results <- stats::setNames(
  vector("list", length(screen_method_args)),
  names(screen_method_args)
)
for (method in names(screen_method_args)) {
  checkpoint <- screen_checkpoint_path(method)
  if (file.exists(checkpoint)) {
    screen_results[[method]] <- readRDS(checkpoint)
    message(glue::glue("Loaded checkpoint for {method}."))
  }
}

methods_to_run <- setdiff(names(screen_method_args), names(purrr::compact(screen_results)))

for (method in methods_to_run) {
  args <- screen_method_args[[method]]
  out <- do.call(SigBridgeR::Screen, c(screen_common, args))

  screen_results[[method]] <- out
  saveRDS(out, screen_checkpoint_path(method))
  message(glue::glue("Screening completed for {method}."))
}
saveRDS(screen_results, file.path(output_dir, "screen_results.rds"))

### 3.2.1 Screening QC across the active methods ----

screen_qc <- purrr::imap_dfr(METHOD_COLUMNS, function(column, method) {
  tab <- screen_results[[method]]$scRNA_data[[]][[column]] |>
    table(useNA = "ifany")
  tibble::tibble(
    method = method,
    status = names(tab),
    n_cells = as.integer(tab)
  )
})
readr::write_csv(screen_qc, file.path(output_dir, "screening_cell_counts.csv"))


# ============================================================================
# 4. Visualisations ----
# ============================================================================

## 4.1 Merge results and equally weighted consensus vote ----

# Collapses any missing/empty calls to `Neutral` so all methods share the same
# three-level Positive/Negative/Neutral scale.
normalise_screen_status <- function(x) {
  x <- dplyr::if_else(
    is.na(x) | x %in% c("", "Other"),
    "Neutral",
    as.character(x)
  )
  factor(x, levels = c("Negative", "Neutral", "Positive"))
}

# `MergeResult()` calls `as.data.table(..., keep.rownames = "cell_id")`; drop any
# pre-existing `cell_id` metadata column to avoid duplicate-name errors.
strip_obs_cell_id <- function(screen_obj) {
  if ("cell_id" %in% colnames(screen_obj$scRNA_data[[]])) {
    screen_obj$scRNA_data$cell_id <- NULL
  }
  screen_obj
}

merged <- SigBridgeR::MergeResult(
  strip_obs_cell_id(screen_results$scPP),
  strip_obs_cell_id(screen_results$LP_SGL)
  # strip_obs_cell_id(screen_results$DEGAS) # disabled
  # strip_obs_cell_id(screen_results$Scissor) # disabled
)
for (column in unname(METHOD_COLUMNS)) {
  merged[[column]] <- normalise_screen_status(merged[[]][[column]])
}

SigBridgeR::setFuncOption(seed = CFG$seed)
merged$WeightedVote <- merged[[]][, unname(METHOD_COLUMNS), drop = FALSE] |>
  stats::setNames(names(METHOD_COLUMNS)) |>
  SigBridgeR::WeightedVote(
    weights = SCREEN$consensus$method_weights[names(METHOD_COLUMNS)],
    ties.method = SCREEN$consensus$ties.method
  ) |>
  normalise_screen_status()

saveRDS(merged, file.path(output_dir, "merged_weighted_vote.rds"))

## 4.2 UpSet plot of positive-cell intersections ----

# `ScreenUpset()` fails with chk >= 0.10 (`FUN.VALUE` must be a vector); use
# the same ggupset logic directly.
build_screen_upset <- function(
    screened_seurat,
    screen_type,
    order_by = c("freq", "degree"),
    n_intersections = 20L
) {
  meta_data <- screened_seurat[[]]
  max_comb <- length(screen_type)
  all_combinations <- unlist(
    lapply(seq_len(max_comb), function(k) {
      if (length(screen_type) >= k) {
        combs <- utils::combn(screen_type, k, simplify = FALSE)
        stats::setNames(
          combs,
          vapply(combs, function(comb) {
            if (length(comb) == 1L) {
              comb
            } else {
              paste(comb, collapse = " & ")
            }
          }, character(1))
        )
      }
    }),
    recursive = FALSE
  )
  positive_matrix <- as.matrix(
    meta_data[, screen_type, drop = FALSE] == "Positive"
  )
  counts <- vapply(all_combinations, function(sets) {
    row_matches <- Matrix::rowSums(
      positive_matrix[, sets, drop = FALSE]
    ) == length(sets)
    sum(row_matches, na.rm = TRUE)
  }, numeric(1))
  intersection_data <- tibble::tibble(
    intersection = names(all_combinations),
    sets = all_combinations,
    count = counts
  )
  ggplot2::ggplot(
    intersection_data,
    ggplot2::aes(x = sets, y = count)
  ) +
    ggplot2::geom_col(fill = "#4E79A7", alpha = 0.9, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.5) +
    ggupset::scale_x_upset(
      order_by = order_by,
      sets = screen_type,
      n_intersections = n_intersections
    ) +
    ggupset::theme_combmatrix() +
    ggplot2::theme_minimal()
}

upset_result <- list(
  plot = build_screen_upset(
    screened_seurat = merged,
    screen_type = unname(METHOD_COLUMNS),
    order_by = c("freq", "degree"),
    n_intersections = 20L
  )
)
ggplot2::ggsave(
  filename = file.path(output_dir, glue::glue("{study}_screening_upset.pdf")),
  plot = upset_result$plot,
  width = 12,
  height = 7
)

## 4.3 Stacked-bar fractions by the existing cluster annotation ----

fraction_result <- SigBridgeR::ScreenFractionPlot(
  screened_seurat = merged,
  group_by = COL$cell_annotation,
  screen_type = c(unname(METHOD_COLUMNS), "WeightedVote"),
  show_null = TRUE,
  show_plot = FALSE,
  plot_title = "Phenotype-guided cell fractions",
  ncol = 2L,
  verbose = TRUE
)
ggplot2::ggsave(
  filename = file.path(
    output_dir,
    glue::glue("{study}_screening_fractions_by_cluster.pdf")
  ),
  plot = fraction_result$combined_plot,
  width = 14,
  height = 10
)

## 4.4 UMAP: existing annotation, two methods, and consensus (shared legend) ----

# Discrete screening-call palette shared by the method panels and the
# consensus panel (all use the same Negative/Neutral/Positive scale).
status_colours <- c(
  Negative = "#386C9B",
  Neutral = "#CECECE",
  Positive = "#FF3333"
)

# Builds one screening-status UMAP panel (no legend); the legend is extracted
# once and shared across all such panels (see below).
build_status_umap <- function(column, title) {
  Seurat::DimPlot(
    merged,
    reduction = "umap",
    group.by = column,
    cols = status_colours,
    raster = TRUE
  ) +
    ggplot2::ggtitle(title)
}

original_panel <- Seurat::DimPlot(
  merged,
  reduction = "umap",
  group.by = COL$cell_annotation,
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggplot2::ggtitle("Existing luque_cluster_annotation") +
  Seurat::NoLegend()

status_panels <- c(
  purrr::imap(METHOD_COLUMNS, build_status_umap),
  list(
    consensus = build_status_umap("WeightedVote", "Equally weighted consensus")
  )
)

# One shared Negative/Neutral/Positive legend for the status panels
# (two methods + consensus), extracted with `cowplot::get_legend()`.
shared_status_legend <- cowplot::get_legend(status_panels[[1L]])
status_panels_no_legend <- lapply(status_panels, function(p) {
  p + Seurat::NoLegend()
})

umap_grid <- cowplot::plot_grid(
  plotlist = c(list(original_panel), status_panels_no_legend),
  ncol = 3L
)
umap_page <- cowplot::plot_grid(
  cowplot::ggdraw() +
    cowplot::draw_label(
      glue::glue("{study} | neural-biased morphotype screening"),
      fontface = "bold",
      size = 14
    ),
  umap_grid,
  shared_status_legend,
  ncol = 1L,
  rel_heights = c(0.06, 1, 0.08)
)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    glue::glue("{study}_screening_umap_panels.pdf")
  ),
  plot = umap_page,
  width = 18,
  height = 10
)

message("SigBridgeR bulk-guided single-cell screening completed.")
