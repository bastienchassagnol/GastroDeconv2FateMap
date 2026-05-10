library(omnideconv)
source("R/utils.R")
# ==========================================================================
# 0. Command-line configuration; and parse arguments ----
# ==========================================================================

# run this instruction in non-interactive way
# both stdout and stderr are redirected to the log file
# nohup Rscript --no-save --no-restore scripts/02_02_omnideconv.R --config_path configs/omnideconv.yml 2>&1 | tee "logs/suppinger_$(date +%F)_omnideconv.log"
# parse arguments
parser <- argparse::ArgumentParser()
parser$add_argument("--config_path", type = "character", required = FALSE)
args <- parser$parse_args()

config_path <- args$config_path
if (is.null(config_path)) {
  config_path <- "configs/omnideconv.yml"
}
if (!file.exists(config_path)) {
  stop(glue::glue("Configuration file not found: {config_path}"))
}
# 0.1. Load configuration ----

config <- yaml::read_yaml(config_path)

study <- config$run$study
if (is.null(study) || identical(study, "")) {
  study <- "suppinger"
}
today <- config$run$date
if (is.null(today) || identical(today, "")) {
  today <- format(Sys.Date(), "%Y-%m-%d")
}
# 0.3. Load input data ----
bulk_omics_path <- config$input$bulk_omics
single_cell_omics_path <- config$input$single_cell_omics
if (is.null(bulk_omics_path) || is.null(single_cell_omics_path)) {
  stop("Both `input.bulk_omics` and `input.single_cell_omics` are required.")
}
bulk_time_point_column <- config$columns$bulk_time_point
if (is.null(bulk_time_point_column) || identical(bulk_time_point_column, "")) {
  bulk_time_point_column <- "time_point_id"
}

# 0.4. Load metadata columns ----
bulk_treatment_column <- config$columns$bulk_treatment
if (is.null(bulk_treatment_column) || identical(bulk_treatment_column, "")) {
  bulk_treatment_column <- "treatment_status"
}
bulk_control_value <- config$filters$bulk_treatment_value
if (is.null(bulk_control_value) || identical(bulk_control_value, "")) {
  bulk_control_value <- "control"
}

single_cell_time_point_column <- config$columns$single_cell_time_point
if (
  is.null(single_cell_time_point_column) ||
    identical(single_cell_time_point_column, "")
) {
  single_cell_time_point_column <- "timepoints"
}
cell_type_column <- config$columns$cell_type
if (is.null(cell_type_column) || identical(cell_type_column, "")) {
  cell_type_column <- "celltypeannotation"
}
batch_column <- config$columns$batch
if (is.null(batch_column) || identical(batch_column, "")) {
  batch_column <- "batch"
}

# 0.4. Deconvolution parameters ----
random_seed <- as.integer(config$run$random_seed)
if (is.na(random_seed)) {
  random_seed <- 123L
}
set.seed(random_seed)
max_cells_per_celltype <- as.integer(config$sampling$max_cells_per_celltype)
if (is.na(max_cells_per_celltype)) {
  max_cells_per_celltype <- 200L
}

method_names <- config$deconvolution$methods
if (is.null(method_names) || length(method_names) == 0L) {
  method_names <- c("BayesPrism", "DWLS", "MuSiC", "Scaden", "SCDC")
}
unknown_methods <- setdiff(
  method_names,
  names(omnideconv::deconvolution_methods)
)
if (length(unknown_methods) > 0L) {
  stop(glue::glue(
    "Unknown deconvolution method(s): {paste(unknown_methods, collapse = ', ')}"
  ))
}
deconv_methods <- omnideconv::deconvolution_methods[method_names]


normalize_results <- config$deconvolution$normalize_results
if (is.null(normalize_results)) {
  normalize_results <- TRUE
}
normalize_results <- isTRUE(normalize_results)

verbose <- config$run$verbose
if (is.null(verbose)) {
  verbose <- TRUE
}
verbose <- isTRUE(verbose)


# 0.5. Scaden-specific parameters ----
scaden_model_path <- config$scaden$model_path
if (is.null(scaden_model_path) || identical(scaden_model_path, "")) {
  scaden_model_path <- ".models/scaden_model_{study}"
}
scaden_model_path <- normalizePath(
  glue::glue(scaden_model_path),
  mustWork = FALSE
)
if (!dir.exists(dirname(scaden_model_path))) {
  dir.create(dirname(scaden_model_path), recursive = TRUE)
}
scaden_temp_dir <- config$scaden$temp_dir
if (is.null(scaden_temp_dir) || identical(scaden_temp_dir, "")) {
  scaden_temp_dir <- tempdir()
}
scaden_temp_dir <- normalizePath(
  glue::glue(scaden_temp_dir),
  mustWork = FALSE
)
if (!dir.exists(scaden_temp_dir)) {
  dir.create(scaden_temp_dir, recursive = TRUE)
}
scaden_batch_size <- as.integer(config$scaden$batch_size)
if (is.na(scaden_batch_size)) {
  scaden_batch_size <- 128L
}
scaden_learning_rate <- as.numeric(config$scaden$learning_rate)
if (is.na(scaden_learning_rate)) {
  scaden_learning_rate <- 1e-04
}
scaden_steps <- as.integer(config$scaden$steps)
if (is.na(scaden_steps)) {
  scaden_steps <- 500L
}
scaden_var_cutoff <- as.numeric(config$scaden$var_cutoff)
if (is.na(scaden_var_cutoff)) {
  scaden_var_cutoff <- 1
}
scaden_cells <- as.integer(config$scaden$cells)
if (is.na(scaden_cells)) {
  scaden_cells <- 50L
}
scaden_samples <- as.integer(config$scaden$samples)
if (is.na(scaden_samples)) {
  scaden_samples <- 1000L
}

# 0.5b. DWLS-specific parameters ----
dwls_method <- config$dwls$signature_method
if (is.null(dwls_method) || identical(dwls_method, "")) {
  dwls_method <- "mast_optimized"
}
dwls_submethod <- config$dwls$deconvolution_submethod
if (is.null(dwls_submethod) || identical(dwls_submethod, "")) {
  dwls_submethod <- "DampenedWLS"
}
# 0.6. Output file ----
output_file <- config$output$file
if (is.null(output_file) || identical(output_file, "")) {
  output_file <- "outputs/deconvolution/{study}_deconvolution_results_{today}.rds"
}
output_file <- glue::glue(output_file)
if (!dir.exists(dirname(output_file))) {
  dir.create(dirname(output_file), recursive = TRUE)
}

# ==========================================================================
# 1. Load preprocessed objects ----
# ==========================================================================

bulk_omics <- readRDS(bulk_omics_path)
single_cell_omics <- readRDS(single_cell_omics_path)

# 1.1. Pre-loop filtering: drop early-treatment bulk samples ----
bulk_omics <- bulk_omics[,
  SummarizedExperiment::colData(bulk_omics)[[bulk_treatment_column]] ==
    bulk_control_value
]
# Shared time points between control bulk and single-cell
shared_timepoints <- intersect(
  unique(as.character(
    SummarizedExperiment::colData(bulk_omics)[[bulk_time_point_column]]
  )),
  unique(as.character(
    SummarizedExperiment::colData(single_cell_omics)[[
      single_cell_time_point_column
    ]]
  ))
)
message("Shared time points: ", paste(shared_timepoints, collapse = ", "))

# 1.2. Balanced single-cell subsampling ----
# Sample up to MAX_CELLS_PER_CELLTYPE cells per cell type × time point,
# without replacement; groups smaller than the cap return all their cells.
sampled_barcodes <- SummarizedExperiment::colData(single_cell_omics) |>
  as.data.frame() |>
  tibble::rownames_to_column("barcode") |>
  dplyr::select(
    barcode,
    dplyr::all_of(c(single_cell_time_point_column, cell_type_column))
  ) |>
  dplyr::group_by(
    dplyr::across(
      dplyr::all_of(c(single_cell_time_point_column, cell_type_column))
    )
  ) |>
  dplyr::slice_sample(
    n = max_cells_per_celltype,
    replace = FALSE
  ) |>
  dplyr::ungroup() |>
  dplyr::pull(barcode)

single_cell_omics_sampled <- single_cell_omics[, sampled_barcodes]

# ==========================================================================
# 4. Deconvolution ----
# ==========================================================================

# Outer map over time points; inner map over methods.
# list_rbind at each level adds the grouping variable as a column
# ("time_point" and "deconvolution_algorithm" respectively).
deconv_results <- purrr::map(
  purrr::set_names(shared_timepoints),
  \(tp) {
    message("\n── Deconvolution at time point: ", tp, " ──────────────────────")

    # --- 4a. Subset bulk to current time point (genes × samples matrix) ----
    bulk_tp <- bulk_omics[,
      SummarizedExperiment::colData(bulk_omics)[[bulk_time_point_column]] == tp
    ]
    bulk_matrix <- as.matrix(SummarizedExperiment::assay(bulk_tp, "counts"))

    # --- 4b. Subset (already balanced) single-cell to current time point ---
    sc_tp <- single_cell_omics_sampled[,
      SummarizedExperiment::colData(single_cell_omics_sampled)[[
        single_cell_time_point_column
      ]] ==
        tp
    ]

    # --- 4c. Run each deconvolution method ---------------------------------
    # Scaden requires a two-step build + predict workflow; all other methods
    # go through omnideconv:::deconvolute directly.
    purrr::map(
      deconv_methods,
      \(method) {
        message(
          "\n── With the following deconvolution algorithm: ",
          method,
          " ──"
        )
        tryCatch(
          {
            # leverage the fact that in R the if codition is an expression
            # rather than statatement without returning any value
            # correct the if condition to return the same output as the else statement
            # see fucntional paradigm here: https://adv-r.hadley.nz/fp.html
            # same logic for the purrr package
            raw_result <- if (method == "scaden") {
              scaden_path <- omnideconv::build_model_scaden(
                single_cell_object = as.matrix(
                  SummarizedExperiment::assay(sc_tp, "counts")
                ),
                cell_type_annotations = SummarizedExperiment::colData(
                  sc_tp
                )[[cell_type_column]],
                bulk_gene_expression = bulk_matrix,
                model_path = scaden_model_path,
                temp_dir = scaden_temp_dir,
                batch_size = scaden_batch_size,
                learning_rate = scaden_learning_rate,
                steps = scaden_steps,
                var_cutoff = scaden_var_cutoff,
                cells = scaden_cells,
                samples = scaden_samples,
                dataset_name = paste0("scaden_", study),
                verbose = verbose
              )
              omnideconv::deconvolute_scaden(
                signature = scaden_path,
                bulk_gene_expression = bulk_matrix,
                verbose = verbose
              )
            } else if (method == "dwls") {
              # DWLS::DEAnalysisMASTOptimized uses matrix arithmetic (`+`);
              # SingleCellExperiment is not numeric.
              sc_mat <- as.matrix(
                SummarizedExperiment::assay(sc_tp, "counts")
              )
              rownames(sc_mat) <- SummarizedExperiment::rownames(sc_tp)
              colnames(sc_mat) <- SummarizedExperiment::colnames(sc_tp)
              signature_matrix <- omnideconv::build_model_dwls(
                single_cell_object = sc_mat,
                cell_type_annotations = SummarizedExperiment::colData(sc_tp)[[
                  cell_type_column
                ]],
                dwls_method = dwls_method,
                verbose = verbose
              )
              raw_dwls <- omnideconv:::deconvolute_dwls(
                bulk_gene_expression = bulk_matrix,
                signature = signature_matrix,
                dwls_submethod = dwls_submethod,
                verbose = verbose
              )
              normalise_cell_estimates(raw_dwls)
            } else {
              if (
                method == "music" &&
                  dplyr::n_distinct(
                    SummarizedExperiment::colData(sc_tp)[[batch_column]]
                  ) <
                    2L
              ) {
                warning("MuSiC requires at least two single-cell subjects.")
              }

              omnideconv:::deconvolute(
                bulk_gene_expression = bulk_matrix,
                method = method,
                single_cell_object = sc_tp,
                batch_ids = SummarizedExperiment::colData(
                  sc_tp
                )[[batch_column]],
                cell_type_column_name = cell_type_column,
                assay_name = "counts",
                normalize_results = normalize_results,
                verbose = verbose
              )
            }

            # Tidy proportion matrix → long-format tibble
            tidy_result <- as.data.frame(raw_result) |>
              tibble::rownames_to_column("sample_id") |>
              tidyr::pivot_longer(
                cols = -sample_id,
                names_to = "cell_type",
                values_to = "cell_proportion"
              )
            tidy_result
          },
          error = function(e) {
            message(
              "  ERROR in ",
              method,
              " at ",
              tp,
              ": ",
              conditionMessage(e)
            )
            NULL
          }
        )
      }
    ) |>
      purrr::compact() |>
      purrr::list_rbind(names_to = "deconvolution_algorithm")
  }
) |>
  purrr::list_rbind(names_to = "time_point") |>
  dplyr::select(
    deconvolution_algorithm,
    time_point,
    sample_id,
    cell_type,
    cell_proportion
  )

# ==========================================================================
# 4. Save deconvolution results ----
# ==========================================================================

saveRDS(
  deconv_results,
  file = output_file
)
