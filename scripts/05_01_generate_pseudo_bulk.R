# ==========================================================================
# 0. Libraries and filename settings ----
# ==========================================================================

library(Seurat)
library(SummarizedExperiment)
library(tinytable)
library(ggplot2)

study <- "luque"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "outputs/phenotype-aware"
source("./R/utils.R")
source("./R/simulate_pseudo_bulk_samples.R")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

n_simulated <- 100L
set.seed(42)


# ==========================================================================
# 1. Load single-cell object ----
# ==========================================================================

GSE250136_120h <- readRDS(
  file = "./data/intermediate/luque_single_cell_merged_2026-06-07.rds"
) |>
  subset(idents = "120h")

dim(GSE250136_120h)

phenotype_data <- GSE250136_120h@meta.data
tinytable::tt(
  phenotype_data |>
    dplyr::group_by(.data$timepoint, .data$Morphotype, .data$Sample.barcode) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(.data$timepoint, .data$Sample.barcode, .data$Morphotype) |>
    tidyr::drop_na(),
  caption = "Number of morphotype annotations per time point and sample barcode"
)


# ==========================================================================
# 2. Strategy 1 — direct barcode aggregation ----
# ==========================================================================

naive_pseudo_bulk <- aggregate_barcode_pseudo_bulk(
  seurat_obj = GSE250136_120h,
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode"
)

saveRDS(
  naive_pseudo_bulk,
  file = glue::glue(
    "./data/intermediate/{study}_naive_pseudo_bulk_{today}.rds"
  )
)


# ==========================================================================
# 3. Strategy 2 — phenotype-stratified bootstrap ----
# ==========================================================================

pb_boot_bio <- simulate_bootstrap_samples(
  seurat_obj = GSE250136_120h,
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  n_samples = n_simulated,
  replicate_type = "biological",
  seed = 42L
)

message(
  "Strategy 2 (bootstrap biological): ",
  ncol(pb_boot_bio),
  " samples; phenotypes: ",
  paste(table(SummarizedExperiment::colData(pb_boot_bio)$Morphotype),
        collapse = ", ")
)


# ==========================================================================
# 4. Strategy 3 — organoid-preserving bootstrap ----
# ==========================================================================

pb_boot_tech <- simulate_bootstrap_samples(
  seurat_obj = GSE250136_120h,
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  n_samples = n_simulated,
  replicate_type = "technical",
  seed = 42L
)

message(
  "Strategy 3 (bootstrap technical): ",
  ncol(pb_boot_tech),
  " samples; phenotypes: ",
  paste(table(SummarizedExperiment::colData(pb_boot_tech)$Morphotype),
        collapse = ", ")
)


# ==========================================================================
# 5. Strategy 7 — generative log-normal pseudo-bulk ----
# ==========================================================================

pb_gen_lognorm <- simulate_generative_models(
  seurat_obj = GSE250136_120h,
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  n_samples = n_simulated,
  model = "lognormal",
  seed = 42L
)

message(
  "Strategy 7 (generative log-normal): ",
  ncol(pb_gen_lognorm),
  " samples; phenotypes: ",
  paste(table(SummarizedExperiment::colData(pb_gen_lognorm)$Morphotype),
        collapse = ", ")
)


# ==========================================================================
# 6. Strategy 10 — Splatter-style negative-binomial pseudo-bulk ----
# ==========================================================================

pb_gen_nb <- simulate_generative_models(
  seurat_obj = GSE250136_120h,
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  n_samples = n_simulated,
  model = "negative_binomial",
  seed = 42L
)

message(
  "Strategy 10 (generative negative binomial): ",
  ncol(pb_gen_nb),
  " samples; phenotypes: ",
  paste(table(SummarizedExperiment::colData(pb_gen_nb)$Morphotype),
        collapse = ", ")
)


# ==========================================================================
# 7. Save SummarizedExperiment objects ----
# ==========================================================================

pseudo_bulk_outputs <- list(
  barcode_aggregation = pb_barcode,
  bootstrap_biological = pb_boot_bio,
  bootstrap_technical = pb_boot_tech,
  generative_lognormal = pb_gen_lognorm,
  generative_negative_binomial = pb_gen_nb
)

saveRDS(
  pseudo_bulk_outputs,
  file = glue::glue(
    "./data/intermediate/{study}_pseudo_bulk_simulations_{today}.rds"
  )
)

message("Saved pseudo-bulk simulations to data/intermediate/")
