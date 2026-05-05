library(omnideconv)

omnideconv::bseqsc_config(file = "./R/CIBERSORT.R", error = TRUE)

# ==========================================================================
# 0. Settings ----
# ==========================================================================

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
scaden_model_path <- normalizePath(
  paste0("./.models/scaden_model_", study),
  mustWork = FALSE
)

deconv_methods <- c(
  "bayesprism",
  "dwls",
  "scdc"
)

# ==========================================================================
# 1. Load preprocessed objects ----
# ==========================================================================

bulk_omics <- readRDS(
  "data/intermediate/suppinger_bulk_summarized_experiment_2026-04-30.rds"
)
single_cell_omics <- readRDS(
  "data/intermediate/suppinger_single_cell_2026-04-30.rds"
)

# ==========================================================================
# 2. Pre-loop filtering: drop early-treatment bulk samples ----
# ==========================================================================

bulk_omics_ctl <- bulk_omics[,
  bulk_omics$treatment_status == "control"
]

# Shared time points between control bulk and single-cell
shared_timepoints <- intersect(
  levels(bulk_omics_ctl$time_point_id),
  levels(single_cell_omics$timepoints)
)

message("Shared time points: ", paste(shared_timepoints, collapse = ", "))

# ==========================================================================
# 3. Deconvolution loop ----
# ==========================================================================

deconv_results <- vector(mode = "list", length = length(shared_timepoints))
names(deconv_results) <- shared_timepoints

for (tp in shared_timepoints) {
  message("\n── Deconvolution at time point: ", tp, " ──────────────────────")

  # --- 3a. Subset bulk to current time point (genes × samples matrix) ----
  bulk_tp <- bulk_omics_ctl[,
    bulk_omics_ctl$time_point_id == tp
  ]
  bulk_matrix <- SummarizedExperiment::assay(bulk_tp, "counts")

  # --- 3b. Subset single-cell to current time point ----------------------
  sc_tp <- single_cell_omics[,
    single_cell_omics$timepoints == tp
  ]

  # --- 3c. Run each deconvolution method ---------------------------------
  # batch_ids <- rep_len(c("1", "2"), ncol(sc_tp))
  for (method in deconv_methods[2:3]) {
    message("\n── With the following deconvolution algorithm: ", method, " ──")
    deconv_results[[tp]][[method]] <- tryCatch(
      omnideconv:::deconvolute(
        bulk_gene_expression = bulk_matrix,
        method = method,
        single_cell_object = sc_tp[,1:200],
        cell_type_column_name = "celltypeannotation",
        assay_name = "counts",
        normalize_results = TRUE,
        verbose = TRUE
      ), 
      error = function(e) {
        message("  ERROR in ", method, " at ", tp, ": ", conditionMessage(e))
      }
    )
  }
}
# ==========================================================================
# 4. Save results ----
# ==========================================================================

saveRDS(
  deconv_results,
  file = paste0(
    "outputs/deconvolution/",
    study,
    "_deconvolution_results_",
    today,
    ".rds"
  )
)



# scaden model building and deconvolution
# if (!dir.exists(dirname(scaden_model_path))) {
#   dir.create(dirname(scaden_model_path), recursive = TRUE)
# }
# scaden_model_path_second <- omnideconv::build_model_scaden(
#   single_cell_object = as.matrix(SummarizedExperiment::assay(
#     sc_tp[, 1:200],
#     "counts"
#   )),
#   cell_type_annotations = SummarizedExperiment::colData(
#     sc_tp
#   )$celltypeannotation[1:200],
#   bulk_gene_expression = bulk_matrix,
#   model_path = scaden_model_path,
#   temp_dir = tempdir(),
#   batch_size = 128,
#   learning_rate = 1e-04,
#   steps = 500,
#   var_cutoff = 1,
#   cells = 50,
#   samples = 1000,
#   dataset_name = paste0("scaden_", study),
#   verbose = TRUE
# )

# test_scaden <- omnideconv::deconvolute_scaden(
#   signature = scaden_model_path_second,
#   bulk_gene_expression = bulk_matrix,
#   verbose = TRUE
# )
