# ============================================================================
# LLM-guided interpretation of clusterProfiler enrichment results
# ============================================================================
#
# Reference:
# https://yulab-smu.top/biomedical-knowledge-mining-book/interpretation.html
#
# Loads the compact RDS from `scripts/06_03_enrichment_analyses.R` and runs
# `interpret_agent()` on every cell type with non-empty enrichment results.
#
# API credentials: create a local `.env` at the project root (git-ignored):
#   OPENAI_API_KEY="sk-..."
#   LLM_MODEL="gpt-5.4-mini"
# `aisdk` reads these automatically; ChatGPT Pro is separate from API billing.
#
# nohup Rscript --no-save --no-restore \
#   scripts/06_03_enrichment_analyses_llm_guided.R \
#   > "logs/GSE250136_$(date +%F)_llm_interpretation.log" 2>&1 &
#
# ============================================================================

# ============================================================================
# 0. Parameters ----
# ============================================================================

study <- "GSE250136"
de_method <- "pseudobulk_deseq2"
output_prefix <- paste(study, de_method, sep = "_")
enrichment_dir <- "outputs/biological-exploration/enrichment-analyses"
output_dir <- file.path(enrichment_dir, "llm-interpretation")
results_rds_file <- file.path(
  enrichment_dir,
  paste0(output_prefix, "_clusterprofiler_enrichment_results.rds")
)

n_pathways <- 50L
add_ppi <- TRUE
interpret_task <- "interpretation"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(aisdk)
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(glue)
  library(ggplot2)
})

message("LLM model: ", aisdk::get_model())

# ============================================================================
# 1. Prior knowledge and biological context ----
# ============================================================================

luque_celltype_priors <- c(
  neural = paste(
    "Luque gastruloid atlas: neural / neuroectoderm-associated population.",
    "Expect neural induction, anterior-posterior patterning, and synaptic",
    "or neuronal differentiation programmes when comparing morphotypes."
  ),
  neuromesodermal.progenitors = paste(
    "Luque gastruloid atlas: neuromesodermal progenitor (NMP-like) state.",
    "Bipotent posterior progenitors linking spinal neural and paraxial",
    "mesoderm fates; Wnt/FGF-responsive axial progenitor biology."
  ),
  pluripotent = paste(
    "Luque gastruloid atlas: pluripotent / epiblast-like population.",
    "Pre-gastrulation competence, self-renewal, and exit from naïve",
    "pluripotency rather than terminal differentiation."
  ),
  somitic = paste(
    "Luque gastruloid atlas: somitic / paraxial mesoderm population.",
    "Segmentation clock, somitogenesis, and musculoskeletal mesoderm",
    "differentiation signatures."
  ),
  unknown = paste(
    "Luque gastruloid atlas: low-confidence or unassigned cluster.",
    "Interpret cautiously; prioritise robust pathway concordance over",
    "lineage naming."
  )
)

#' Build the LLM context string for one cell type.
#'
#' Combines study design, Luque prior knowledge, and top enrichment headlines
#' from `summary_tbl` so the agent focuses on mechanism rather than relabelling.
build_interpretation_context <- function(
  celltype_level,
  parameters,
  summary_tbl
) {
  prior <- luque_celltype_priors[[celltype_level]]
  if (is.null(prior)) {
    prior <- glue::glue(
      "Luque gastruloid cell type `{celltype_level}`; no curated prior text."
    )
  }

  ct_summary <- summary_tbl |>
    dplyr::filter(.data$celltype == celltype_level) |>
    dplyr::arrange(.data$method, .data$database)

  top_hits <- if (nrow(ct_summary) > 0L) {
    ct_summary |>
      dplyr::mutate(
        line = glue::glue(
          "{method}/{database}: {top_term} (FDR={signif(top_FDR, 3)})"
        )
      ) |>
      dplyr::pull(.data$line)
  } else {
    character()
  }

  glue::glue(
    "Study: {parameters$study}. ",
    "Pseudo-bulk DESeq2 per Luque cell type at 120 h gastruloids. ",
    "Contrast: neural_bias vs TLS morphotype ",
    "(positive log2FC = higher in neural_bias). ",
    "Task: {interpret_task} only — explain mechanisms and implications; ",
    "do not perform cell-type assignment or relabelling. ",
    "Prior biological knowledge for this cluster: {prior}. ",
    if (length(top_hits) > 0L) {
      paste0(
        "Top enrichment headlines: ",
        paste(top_hits, collapse = "; "),
        "."
      )
    } else {
      ""
    }
  )
}

# ============================================================================
# 2. Helpers: enrichment objects and DESeq2 fold changes ----
# ============================================================================

has_enrichment_terms <- function(x) {
  !is.null(x) &&
    inherits(x, c("enrichResult", "gseaResult")) &&
    nrow(as.data.frame(x)) > 0L
}

collect_interpret_objects <- function(celltype_result) {
  objs <- list()
  for (method in c("ORA", "GSEA")) {
    method_list <- celltype_result[[method]]
    if (is.null(method_list)) {
      next
    }
    for (db in names(method_list)) {
      obj <- method_list[[db]]
      if (!has_enrichment_terms(obj)) {
        next
      }
      objs[[paste(method, db, sep = "_")]] <- obj
    }
  }
  objs
}

#' Extract a named log2 fold-change vector for one cell type from DESeq2 output.
#'
#' Reads `parameters$input_de_file`, keeps one value per gene symbol, and
#' returns the vector passed to `interpret_agent(gene_fold_change = ...)`.
load_gene_fold_changes <- function(
  input_de_file,
  celltype_level,
  fdr_cutoff = 0.10,
  lfc_cutoff = 0.5
) {
  df_degs <- readr::read_csv(input_de_file, show_col_types = FALSE)
  padj_col <- intersect(
    c("padj_celltype_BH", "padj_BH", "padj"),
    names(df_degs)
  )
  if (length(padj_col) == 0L) {
    stop("Could not find an adjusted p-value column in ", input_de_file)
  }
  padj_col <- padj_col[[1L]]

  df_ct <- df_degs |>
    dplyr::filter(
      .data$celltype == celltype_level,
      is.finite(.data$log2FoldChange)
    ) |>
    dplyr::transmute(
      gene_symbol = .data$gene,
      log2_fc = .data$log2FoldChange,
      padj = .data[[padj_col]]
    ) |>
    dplyr::group_by(.data$gene_symbol) |>
    dplyr::slice_max(
      order_by = abs(.data$log2_fc),
      n = 1L,
      with_ties = FALSE
    ) |>
    dplyr::ungroup()

  gene_fold_change <- stats::setNames(df_ct$log2_fc, df_ct$gene_symbol)
  gene_fold_change <- gene_fold_change[!is.na(names(gene_fold_change))]
  gene_fold_change <- gene_fold_change[!duplicated(names(gene_fold_change))]

  de_genes <- df_ct |>
    dplyr::filter(
      .data$padj < fdr_cutoff,
      abs(.data$log2_fc) >= lfc_cutoff
    ) |>
    dplyr::pull(.data$gene_symbol) |>
    unique()

  list(
    gene_fold_change = gene_fold_change,
    n_genes_total = length(gene_fold_change),
    n_de_genes = length(de_genes),
    de_genes = de_genes
  )
}

#' Run the clusterProfiler multi-agent interpreter for one cell type.
#'
#' Calls `interpret_agent()` with enrichment results, biological context,
#' DESeq2 fold changes, and STRING PPI data (`add_ppi = TRUE`).
interpret_one_celltype <- function(
  enrichment_payload,
  context,
  gene_fold_change
) {
  clusterProfiler::interpret_agent(
    x = enrichment_payload,
    context = context,
    n_pathways = n_pathways,
    add_ppi = add_ppi,
    gene_fold_change = gene_fold_change,
    verbose = TRUE
  )
}

save_interpretation_outputs <- function(
  celltype_level,
  interpret_res,
  output_dir,
  output_prefix
) {
  base_name <- paste0(
    output_prefix,
    "_llm_interpret_",
    celltype_level
  )

  report_txt <- file.path(output_dir, paste0(base_name, ".txt"))
  if (is.list(interpret_res) && !is.null(interpret_res$report)) {
    writeLines(interpret_res$report, report_txt)
  } else {
    capture.output(print(interpret_res), file = report_txt)
  }

  report_rds <- file.path(output_dir, paste0(base_name, ".rds"))
  saveRDS(interpret_res, report_rds)

  network_pdf <- file.path(output_dir, paste0(base_name, "_network.pdf"))
  network_saved <- FALSE

  network_saved <- tryCatch(
    {
      p <- plot(interpret_res)
      ggplot2::ggsave(
        filename = network_pdf,
        plot = p,
        width = 12,
        height = 10,
        units = "in",
        limitsize = FALSE
      )
      TRUE
    },
    error = function(e) {
      message(
        "Network plot failed for ",
        celltype_level,
        ": ",
        conditionMessage(e)
      )
      FALSE
    }
  )

  list(
    report_txt = report_txt,
    report_rds = report_rds,
    network_pdf = if (network_saved) network_pdf else NA_character_
  )
}

# ============================================================================
# 3. Load enrichment RDS ----
# ============================================================================

if (!file.exists(results_rds_file)) {
  stop(
    "Enrichment RDS not found: ",
    results_rds_file,
    ". Run scripts/06_03_enrichment_analyses.R first."
  )
}

enrichment_bundle <- readRDS(results_rds_file)
required_rds <- c(
  "parameters",
  "summary_tbl",
  "enrichment_long",
  "enrichment_results"
)
missing_rds <- setdiff(required_rds, names(enrichment_bundle))
if (length(missing_rds) > 0L) {
  stop(
    "RDS is missing components: ",
    paste(missing_rds, collapse = ", ")
  )
}

parameters <- enrichment_bundle$parameters
summary_tbl <- enrichment_bundle$summary_tbl
celltype_results <- enrichment_bundle$enrichment_results
celltypes_to_run <- names(celltype_results)

# ============================================================================
# 4. Run interpretation per cell type ----
# ============================================================================

interpretation_outputs <- list()

for (celltype_level in celltypes_to_run) {
  message("Preparing interpretation for: ", celltype_level)

  enrichment_payload <- collect_interpret_objects(
    celltype_results[[celltype_level]]
  )
  if (length(enrichment_payload) == 0L) {
    message("Skipping ", celltype_level, ": no non-empty enrichment objects.")
    next
  }

  # DESeq2 fold changes: named log2FC vector (neural_bias vs TLS) injected so
  # the Detective agent can infer pathway direction, not just enrichment.
  fc_bundle <- load_gene_fold_changes(
    input_de_file = parameters$input_de_file,
    celltype_level = celltype_level,
    fdr_cutoff = parameters$fdr_cutoff,
    lfc_cutoff = parameters$lfc_cutoff_ora
  )

  # Biological context: study contrast, Luque priors, and top pathway headlines
  # guiding the Storyteller toward mechanistic interpretation, not relabelling.
  context <- build_interpretation_context(
    celltype_level = celltype_level,
    parameters = parameters,
    summary_tbl = summary_tbl
  )

  # Multi-agent LLM call: Cleaner + Detective (PPI + fold changes) + Storyteller.
  interpret_res <- interpret_one_celltype(
    enrichment_payload = enrichment_payload,
    context = context,
    gene_fold_change = fc_bundle$gene_fold_change
  )

  saved <- save_interpretation_outputs(
    celltype_level = celltype_level,
    interpret_res = interpret_res,
    output_dir = output_dir,
    output_prefix = output_prefix
  )
  interpretation_outputs[[celltype_level]] <- saved
  message("Saved: ", saved$report_txt)
}

message(
  "LLM interpretation complete. ",
  length(interpretation_outputs),
  " cell type(s) processed."
)
