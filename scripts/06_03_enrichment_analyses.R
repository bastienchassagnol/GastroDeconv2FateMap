# ============================================================================
# 0. Library and hyperparameter setup ----
# ============================================================================

study <- "GSE250136"
# de_method <- "pseudobulk_deseq2"
de_method <- "mast"
technique <- paste("clusterprofiler_enrichment", de_method, sep = "_")
output_prefix <- paste(study, de_method, sep = "_")
de_input_dir <- "outputs/biological-exploration/DEA-analyses"
results_date <- "2026-07-09"
output_dir <- "outputs/biological-exploration/enrichment-analyses"
results_rds_file <- file.path(
  output_dir,
  paste0(output_prefix, "_clusterprofiler_enrichment_results.rds")
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ReactomePA)
  library(org.Mm.eg.db)
  library(msigdbr)
  library(GOSemSim)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(gridExtra)
})

set.seed(1)

fdr_cutoff <- 0.10
volcano_fdr_cutoff <- 0.05
lfc_cutoff_ora <- 0.5
go_ontologies <- c("BP")
msig_species <- "Mus musculus"
msig_collection <- "M2"
msig_subcollection <- "CP:REACTOME"
reactome_organism <- "mouse"
min_gs_size <- 10
max_gs_size <- 500
show_category <- 20
# input_de_file <- file.path(
#   de_input_dir,
#   paste0(
#     study,
#     "_pseudobulk_deseq2_120h_celltype_model_results_",
#     results_date,
#     ".csv"
#   )
# )
input_de_file <- file.path(
  de_input_dir,
  paste0(
    study,
    "_mast_120h_celltype_model_results_",
    results_date,
    ".csv"
  )
)

# ============================================================================
# 1. Standardise DE table to one row per gene and cell type ----
# ============================================================================

# nohup Rscript --no-save --no-restore \
#   scripts/06_03_enrichment_analyses.R \
#   > "logs/GSE250136_$(date +%F)_enrichment_analyses.log" 2>&1 &

df_degs <- readr::read_csv(input_de_file, show_col_types = FALSE)
required_de_cols <- c(
  "study",
  "technique",
  "analysis_level",
  "gene",
  "celltype",
  "log2FoldChange",
  "pvalue"
)
missing_de_cols <- setdiff(required_de_cols, names(df_degs))
if (length(missing_de_cols) > 0L) {
  stop("Missing required DE columns: ", paste(missing_de_cols, collapse = ", "))
}

padj_col <- intersect(c("padj_celltype_BH", "padj_BH", "padj"), names(df_degs))
if (length(padj_col) == 0L) {
  stop("Could not find an adjusted p-value column in the DE table.")
}
padj_col <- padj_col[[1L]]

de_long <- df_degs |>
  dplyr::transmute(
    study = .data$study,
    de_technique = .data$technique,
    analysis_level = .data$analysis_level,
    gene_symbol = .data$gene,
    celltype = .data$celltype,
    log2_fc = .data$log2FoldChange,
    pvalue = .data$pvalue,
    padj = .data[[padj_col]]
  ) |>
  dplyr::mutate(
    pvalue = tidyr::replace_na(.data$pvalue, 1),
    padj = tidyr::replace_na(.data$padj, 1),
    log2_fc = tidyr::replace_na(.data$log2_fc, 0),
    raw_p_for_rank = pmax(.data$pvalue, .Machine$double.xmin),
    # GSEA ranks use raw p-values as weights. FDR still defines ORA genes.
    rank_score = .data$log2_fc * (-log10(.data$raw_p_for_rank)),
    direction = dplyr::case_when(
      .data$log2_fc > 0 ~ "up_in_neural_bias",
      .data$log2_fc < 0 ~ "down_in_neural_bias",
      TRUE ~ "zero"
    ),
    is_de_for_ora = .data$padj < fdr_cutoff &
      abs(.data$log2_fc) >= lfc_cutoff_ora
  )

# ============================================================================
# 2. Map mouse gene symbols to Entrez IDs ----
# ============================================================================

# clusterProfiler, ReactomePA, and MSigDB Reactome are most stable here when
# keyed by Entrez IDs. One-to-many mappings are retained for ORA, while GSEA
# keeps the strongest absolute ranking score for each Entrez ID.
gene_map <- clusterProfiler::bitr(
  unique(de_long$gene_symbol),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db::org.Mm.eg.db
) |>
  tibble::as_tibble() |>
  dplyr::distinct(.data$SYMBOL, .data$ENTREZID) |>
  dplyr::mutate(ENTREZID = as.character(.data$ENTREZID))

de_annot <- de_long |>
  dplyr::inner_join(gene_map, by = c("gene_symbol" = "SYMBOL")) |>
  dplyr::filter(!is.na(.data$ENTREZID), is.finite(.data$rank_score))

message(
  "Mapped genes: ",
  dplyr::n_distinct(de_annot$ENTREZID),
  " / ",
  dplyr::n_distinct(de_long$gene_symbol)
)

# ============================================================================
# 3. Prepare MSigDB Reactome TERM2GENE and TERM2NAME tables ----
# ============================================================================

msig_raw <- msigdbr::msigdbr(
  db_species = "MM",
  species = msig_species,
  collection = msig_collection,
  subcollection = msig_subcollection
)

gene_col <- intersect(
  c("ncbi_gene", "entrez_gene", "entrez_gene_id"),
  colnames(msig_raw)
)[[1]]

msig_reactome <- list(
  term2gene = msig_raw |>
    dplyr::transmute(
      term = .data$gs_name,
      gene = as.character(.data[[gene_col]])
    ) |>
    dplyr::filter(!is.na(.data$gene), .data$gene != "") |>
    dplyr::distinct(),
  term2name = msig_raw |>
    dplyr::transmute(
      term = .data$gs_name,
      name = dplyr::coalesce(.data$gs_description, .data$gs_name)
    ) |>
    dplyr::distinct()
)

# ============================================================================
# 4. Enrichment utilities for result conversion and term similarity ----
# ============================================================================

has_terms <- function(x) {
  !is.null(x) &&
    inherits(x, c("enrichResult", "gseaResult")) &&
    nrow(as.data.frame(x)) > 0
}

object_to_tbl <- function(obj, celltype_level, method, database) {
  if (!has_terms(obj)) {
    return(tibble::tibble())
  }

  as.data.frame(obj) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      study = study,
      technique = technique,
      celltype = celltype_level,
      method = method,
      database = database,
      neglog10_padj = -log10(pmax(.data$p.adjust, .Machine$double.xmin)),
      sig_fdr = .data$p.adjust < fdr_cutoff,
      .before = 1
    )
}

add_term_similarity <- function(obj, database) {
  if (!has_terms(obj)) {
    return(obj)
  }

  out <- tryCatch(
    {
      if (stringr::str_starts(database, "GO_")) {
        ont <- stringr::str_remove(database, "^GO_")
        sem_data <- GOSemSim::godata(
          OrgDb = "org.Mm.eg.db",
          ont = ont,
          computeIC = TRUE
        )
        enrichplot::pairwise_termsim(obj, semData = sem_data)
      } else {
        enrichplot::pairwise_termsim(obj)
      }
    },
    error = function(e) {
      message(
        "Could not compute term similarity for ",
        database,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(out)) obj else out
}

use_adjusted_pvalue_for_plot <- function(obj) {
  if (inherits(obj, "list")) {
    return(lapply(obj, use_adjusted_pvalue_for_plot))
  }
  if (inherits(obj, "compareClusterResult")) {
    if ("p.adjust" %in% colnames(obj@compareClusterResult)) {
      obj@compareClusterResult$pvalue <- obj@compareClusterResult$p.adjust
    }
    return(obj)
  }
  if (!has_terms(obj) || !"p.adjust" %in% colnames(obj@result)) {
    return(obj)
  }
  obj@result$pvalue <- obj@result$p.adjust
  obj
}

# Collect all non-empty enrichment objects for one method across databases.
collect_method_objects <- function(
  method_list,
  method_label,
  needs_termsim = FALSE
) {
  objs <- list()
  for (db in names(method_list)) {
    obj <- method_list[[db]]
    if (!has_terms(obj)) {
      next
    }
    obj <- use_adjusted_pvalue_for_plot(obj)
    obj@result$database <- db
    obj@result$method <- method_label
    # manhattanplot colours terms by ONTOLOGY; use database labels when
    # combining GO, ReactomePA, and MSigDB in one panel.
    obj@result$ONTOLOGY <- db
    if (needs_termsim) {
      obj <- add_term_similarity(obj, db)
    }
    objs[[db]] <- obj
  }
  objs
}

# Merge multiple database-level objects when enrichplot needs one combined input.
merge_method_objects <- function(objs) {
  if (length(objs) == 0L) {
    return(NULL)
  }
  if (length(objs) == 1L) {
    return(objs[[1L]])
  }
  clusterProfiler::merge_result(objs)
}

# ============================================================================
# 5. Run ORA and GSEA per cell type ----
# ============================================================================

run_enrichment_one_celltype <- function(celltype_level) {
  message("Running enrichment for: ", celltype_level)

  # --------------------------------------------------------------------------
  # 5.1 Prepare one cell-type universe, ORA gene set, and GSEA ranked list ----
  # --------------------------------------------------------------------------
  df_celltype <- de_annot |>
    dplyr::filter(.data$celltype == celltype_level)

  universe <- df_celltype |>
    dplyr::distinct(.data$ENTREZID) |>
    dplyr::pull(.data$ENTREZID) |>
    as.character()

  ora_genes <- df_celltype |>
    dplyr::filter(.data$is_de_for_ora) |>
    dplyr::distinct(.data$ENTREZID) |>
    dplyr::pull(.data$ENTREZID) |>
    as.character()

  rank_tbl <- df_celltype |>
    dplyr::filter(is.finite(.data$rank_score), !is.na(.data$ENTREZID)) |>
    dplyr::group_by(.data$ENTREZID) |>
    dplyr::slice_max(
      order_by = abs(.data$rank_score),
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      rank_score = .data$rank_score +
        dplyr::dense_rank(.data$log2_fc) * 1e-12
    ) |>
    dplyr::arrange(dplyr::desc(.data$rank_score))

  gene_list <- rank_tbl$rank_score
  names(gene_list) <- rank_tbl$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)

  res <- list(ORA = list(), GSEA = list())

  # --------------------------------------------------------------------------
  # 5.2 GO over-representation and ranked enrichment analyses ----
  # --------------------------------------------------------------------------
  for (ont in go_ontologies) {
    db_name <- paste0("GO_", ont)
    message("Running GO ORA for: ", celltype_level, " / ", db_name)
    res$ORA[[db_name]] <- tryCatch(
      clusterProfiler::enrichGO(
        gene = ora_genes,
        universe = universe,
        OrgDb = org.Mm.eg.db::org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        readable = TRUE
      ),
      error = function(e) {
        message(
          "GO ORA failed for ",
          celltype_level,
          " / ",
          db_name,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )

    message("Running GO GSEA for: ", celltype_level, " / ", db_name)
    res$GSEA[[db_name]] <- tryCatch(
      clusterProfiler::gseGO(
        geneList = gene_list,
        OrgDb = org.Mm.eg.db::org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        minGSSize = min_gs_size,
        maxGSSize = max_gs_size,
        verbose = FALSE
      ),
      error = function(e) {
        message(
          "GO GSEA failed for ",
          celltype_level,
          " / ",
          db_name,
          ": ",
          conditionMessage(e)
        )
        NULL
      }
    )
  }

  # --------------------------------------------------------------------------
  # 5.3 ReactomePA pathway ORA and GSEA ----
  # --------------------------------------------------------------------------
  message("Running ReactomePA ORA for: ", celltype_level)
  res$ORA[["ReactomePA_Reactome"]] <- tryCatch(
    ReactomePA::enrichPathway(
      gene = ora_genes,
      universe = universe,
      organism = reactome_organism,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      readable = TRUE
    ),
    error = function(e) {
      message(
        "ReactomePA ORA failed for ",
        celltype_level,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )

  message("Running ReactomePA GSEA for: ", celltype_level)
  res$GSEA[["ReactomePA_Reactome"]] <- tryCatch(
    ReactomePA::gsePathway(
      geneList = gene_list,
      organism = reactome_organism,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      verbose = FALSE
    ),
    error = function(e) {
      message(
        "ReactomePA GSEA failed for ",
        celltype_level,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )

  # --------------------------------------------------------------------------
  # 5.4 Mouse MSigDB curated Reactome ORA and GSEA ----
  # --------------------------------------------------------------------------
  message("Running MSigDB Reactome ORA for: ", celltype_level)
  res$ORA[["MSigDB_M2_CP_REACTOME"]] <- tryCatch(
    clusterProfiler::enricher(
      gene = ora_genes,
      universe = universe,
      TERM2GENE = msig_reactome$term2gene,
      TERM2NAME = msig_reactome$term2name,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size
    ),
    error = function(e) {
      message(
        "MSigDB Reactome ORA failed for ",
        celltype_level,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )

  message("Running MSigDB Reactome GSEA for: ", celltype_level)
  res$GSEA[["MSigDB_M2_CP_REACTOME"]] <- tryCatch(
    clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = msig_reactome$term2gene,
      TERM2NAME = msig_reactome$term2name,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      verbose = FALSE
    ),
    error = function(e) {
      message(
        "MSigDB Reactome GSEA failed for ",
        celltype_level,
        ": ",
        conditionMessage(e)
      )
      NULL
    }
  )

  res
}

celltype_results <- stats::setNames(
  lapply(unique(de_annot$celltype), run_enrichment_one_celltype),
  unique(de_annot$celltype)
)

enrichment_long <- purrr::imap_dfr(celltype_results, function(res_ct, ct) {
  dplyr::bind_rows(
    purrr::imap_dfr(res_ct$ORA, ~ object_to_tbl(.x, ct, "ORA", .y)),
    purrr::imap_dfr(res_ct$GSEA, ~ object_to_tbl(.x, ct, "GSEA", .y))
  )
})

summary_tbl <- enrichment_long |>
  dplyr::group_by(.data$celltype, .data$method, .data$database) |>
  dplyr::summarise(
    n_terms_tested = dplyr::n(),
    n_terms_FDR = sum(.data$p.adjust < fdr_cutoff, na.rm = TRUE),
    top_term = .data$Description[which.min(.data$p.adjust)][1],
    top_FDR = min(.data$p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$celltype, .data$method, .data$database)

# ============================================================================
# 6. Native enrichplot visualisations saved with ggsave() ----
# ============================================================================

build_native_page <- function(
  celltype_level,
  plot_fun,
  title,
  needs_termsim = FALSE,
  methods = c("ORA", "GSEA"),
  combine_databases = TRUE
) {
  res_ct <- celltype_results[[celltype_level]]
  plot_list <- list()

  if (combine_databases) {
    for (method in methods) {
      objs <- collect_method_objects(
        method_list = res_ct[[method]],
        method_label = method,
        needs_termsim = needs_termsim
      )

      if (length(objs) == 0L) {
        message(
          "Skipping empty ",
          method,
          " results for ",
          title,
          " / ",
          celltype_level
        )
        next
      }

      panel_title <- paste(celltype_level, method, sep = " | ")
      plot_obj <- tryCatch(
        plot_fun(objs, method) +
          ggplot2::ggtitle(panel_title) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = 18, face = "bold")
          ),
        error = function(e) {
          message(
            "Plot failed for ",
            title,
            " / ",
            celltype_level,
            " / ",
            method,
            ": ",
            conditionMessage(e)
          )
          NULL
        }
      )

      if (!is.null(plot_obj)) {
        plot_list[[method]] <- plot_obj
      }
    }

    if (length(plot_list) == 0L) {
      message("Skipping ", title, " page for ", celltype_level, ".")
      return(NULL)
    }

    return(
      patchwork::wrap_plots(plot_list, ncol = length(plot_list)) +
        patchwork::plot_annotation(
          title = paste0(title, " | Cell type: ", celltype_level),
          subtitle = paste0(
            "Columns: ORA and GSEA; all databases combined within each method; ",
            "FDR focus = ",
            fdr_cutoff
          )
        ) &
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 24, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 16)
        )
    )
  }

  db_order <- unique(c(names(res_ct$ORA), names(res_ct$GSEA)))
  for (method in methods) {
    for (db in db_order) {
      obj <- res_ct[[method]][[db]]
      panel_title <- paste(celltype_level, method, db, sep = " | ")

      if (!has_terms(obj)) {
        next
      }

      obj <- use_adjusted_pvalue_for_plot(obj)
      if (needs_termsim) {
        obj <- add_term_similarity(obj, db)
      }

      plot_obj <- tryCatch(
        plot_fun(obj, method, db) +
          ggplot2::ggtitle(panel_title) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(size = 18, face = "bold")
          ),
        error = function(e) {
          message(
            "Plot failed for ",
            title,
            " / ",
            celltype_level,
            " / ",
            method,
            " / ",
            db,
            ": ",
            conditionMessage(e)
          )
          NULL
        }
      )

      if (!is.null(plot_obj)) {
        plot_list[[paste(method, db, sep = "_")]] <- plot_obj
      }
    }
  }

  if (length(plot_list) == 0L) {
    message("Skipping ", title, " page for ", celltype_level, ".")
    return(NULL)
  }

  patchwork::wrap_plots(plot_list, ncol = length(db_order)) +
    patchwork::plot_annotation(
      title = paste0(title, " | Cell type: ", celltype_level),
      subtitle = paste0(
        "Rows: ORA then GSEA; columns split by ontology/pathway source; ",
        "FDR focus = ",
        fdr_cutoff
      )
    ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 24, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 16)
    )
}

save_native_grid_pdf <- function(
  filename,
  plot_fun,
  title,
  needs_termsim = FALSE,
  methods = c("ORA", "GSEA"),
  combine_databases = TRUE,
  width = 26,
  height = 16
) {
  pages <- lapply(
    names(celltype_results),
    build_native_page,
    plot_fun = plot_fun,
    title = title,
    needs_termsim = needs_termsim,
    methods = methods,
    combine_databases = combine_databases
  ) |>
    purrr::compact()

  if (length(pages) == 0L) {
    message("No pages generated for ", title, "; skipping ", filename)
    return(invisible(NULL))
  }

  ggplot2::ggsave(
    filename = file.path(output_dir, filename),
    plot = gridExtra::marrangeGrob(pages, nrow = 1, ncol = 1),
    width = width,
    height = height,
    units = "in",
    limitsize = FALSE
  )
}

# Dotplot: overview of term strength and gene support across ORA/GSEA sources.
save_native_grid_pdf(
  filename = paste0(output_prefix, "_01_dotplot_ORA_GSEA.pdf"),
  title = "Dotplot",
  plot_fun = function(objs, method) {
    x <- merge_method_objects(objs)
    p <- if (inherits(x, "compareClusterResult")) {
      enrichplot::dotplot(
        x,
        showCategory = show_category,
        label_format = 50,
        split = "Cluster",
        color = "pvalue"
      )
    } else if ("ONTOLOGY" %in% colnames(as.data.frame(x))) {
      enrichplot::dotplot(
        x,
        showCategory = show_category,
        label_format = 50,
        split = "ONTOLOGY",
        color = "pvalue"
      ) +
        ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free")
    } else {
      enrichplot::dotplot(
        x,
        showCategory = show_category,
        label_format = 50,
        color = "pvalue"
      )
    }
    p +
      ggplot2::scale_size_continuous(range = c(5, 13)) +
      ggplot2::scale_colour_continuous(
        name = "BH adjusted p-value",
        labels = function(x) signif(x, 3)
      ) +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 15),
        axis.text.x = ggplot2::element_text(size = 14),
        axis.title = ggplot2::element_text(size = 16),
        legend.text = ggplot2::element_text(size = 13),
        legend.title = ggplot2::element_text(size = 14),
        strip.text = ggplot2::element_text(size = 15, face = "bold")
      )
  },
  width = 28,
  height = 18
)

# Manhattan plot: highlights strongest FDR-adjusted enrichment signals.
save_native_grid_pdf(
  filename = paste0(output_prefix, "_02_manhattanplot_ORA_GSEA.pdf"),
  title = "Manhattan plot",
  plot_fun = function(objs, method) {
    p <- enrichplot::manhattanplot(
      objs,
      showCategory = 10
    )

    p +
      ggplot2::scale_fill_viridis_d(name = "Database") +
      ggplot2::geom_hline(
        yintercept = -log10(fdr_cutoff),
        colour = "blue",
        linetype = "dashed",
        linewidth = 0.6
      ) +
      ggplot2::annotate(
        "text",
        x = -Inf,
        y = -log10(fdr_cutoff) + 0.1,
        label = paste0("BH padj = ", signif(fdr_cutoff, 3)),
        colour = "blue",
        hjust = -0.05,
        size = 5
      ) +
      ggplot2::theme_bw(base_size = 17) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 12),
        axis.title = ggplot2::element_text(size = 15),
        legend.text = ggplot2::element_text(size = 12),
        legend.title = ggplot2::element_text(size = 13)
      )
  },
  width = 30,
  height = 16
)

# Tree plot: clusters semantically similar enriched terms into larger themes.
save_native_grid_pdf(
  filename = paste0(output_prefix, "_03_treeplot_ORA_GSEA.pdf"),
  title = "Tree plot",
  plot_fun = function(x, method, database) {
    n_terms <- nrow(as.data.frame(x))
    n_cluster <- min(6L, max(1L, n_terms - 1L))
    enrichplot::treeplot(
      x,
      nCluster = n_cluster,
      label_format = 45,
      fontsize = 4,
      fontsize_cladelab = 5
    )
  },
  needs_termsim = TRUE,
  combine_databases = FALSE,
  width = 32,
  height = 22
)

# Semantic-space plot: shows term-similarity modules as networks or embeddings.
save_native_grid_pdf(
  filename = paste0(output_prefix, "_04_semantic_space_ORA_GSEA.pdf"),
  title = "Semantic-space plot",
  plot_fun = function(x, method, database) {
    n_terms <- nrow(as.data.frame(x))
    n_cluster <- min(6L, max(1L, n_terms - 1L))
    if ("ssplot" %in% getNamespaceExports("enrichplot")) {
      enrichplot::ssplot(x, nCluster = n_cluster) +
        ggplot2::labs(colour = "BH adjusted p-value")
    } else {
      enrichplot::emapplot(x, node_label = "group", layout = "kk") +
        ggplot2::labs(colour = "BH adjusted p-value")
    }
  },
  needs_termsim = TRUE,
  combine_databases = FALSE,
  width = 34,
  height = 24
)

# ============================================================================
# 7. Volcano-style enrichment plot ----
# ============================================================================

parse_gene_ids <- function(gene_id_string) {
  if (is.na(gene_id_string) || gene_id_string == "") {
    return(character())
  }
  unlist(strsplit(gene_id_string, "/", fixed = TRUE))
}

compute_ora_signed_mean_lfc <- function(gene_id_string, celltype_level) {
  ids <- parse_gene_ids(gene_id_string)
  df_celltype <- de_annot |>
    dplyr::filter(.data$celltype == celltype_level)

  if (length(intersect(ids, df_celltype$gene_symbol)) > 0) {
    mean(df_celltype$log2_fc[df_celltype$gene_symbol %in% ids], na.rm = TRUE)
  } else {
    mean(df_celltype$log2_fc[df_celltype$ENTREZID %in% ids], na.rm = TRUE)
  }
}

plot_volcano_one_celltype <- function(celltype_level, show_legend = FALSE) {
  dat <- enrichment_long |>
    dplyr::filter(.data$celltype == celltype_level) |>
    dplyr::group_by(.data$method, .data$database) |>
    dplyr::arrange(.data$p.adjust, .by_group = TRUE) |>
    dplyr::slice_head(n = 150) |>
    dplyr::ungroup() |>
    dplyr::rowwise() |>
    dplyr::mutate(
      volcano_x = dplyr::case_when(
        .data$method == "GSEA" ~ .data$NES,
        .data$method == "ORA" ~ compute_ora_signed_mean_lfc(
          .data$geneID,
          .data$celltype
        ),
        TRUE ~ NA_real_
      ),
      direction_label = dplyr::case_when(
        .data$volcano_x > 0 ~ "higher in neural_bias",
        .data$volcano_x < 0 ~ "higher in TLS",
        TRUE ~ "neutral"
      ),
      size_for_plot = dplyr::case_when(
        .data$method == "GSEA" ~ as.numeric(.data$setSize),
        .data$method == "ORA" ~ as.numeric(.data$Count),
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(is.finite(.data$volcano_x), is.finite(.data$neglog10_padj)) |>
    dplyr::group_by(.data$method, .data$database) |>
    dplyr::arrange(.data$p.adjust, .by_group = TRUE) |>
    dplyr::mutate(
      label_me = dplyr::row_number() <= 3 |
        .data$p.adjust < volcano_fdr_cutoff,
      label_text = paste0(
        stringr::str_trunc(.data$Description, 55),
        "\nFDR = ",
        signif(.data$p.adjust, 3)
      )
    ) |>
    dplyr::ungroup()

  if (nrow(dat) == 0L) {
    message("No enrichment results for volcano-style plot: ", celltype_level)
    return(NULL)
  }

  y_limit <- max(dat$neglog10_padj, na.rm = TRUE)
  y_limit <- max(y_limit, -log10(volcano_fdr_cutoff) + 0.5)
  fdr_line <- -log10(volcano_fdr_cutoff)

  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = .data$volcano_x,
      y = .data$neglog10_padj,
      colour = .data$direction_label,
      size = .data$size_for_plot
    )
  ) +
    ggplot2::geom_hline(
      yintercept = fdr_line,
      colour = "blue",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_vline(
      xintercept = c(-1, 1),
      colour = "darkgreen",
      linetype = "solid",
      linewidth = 0.5
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = fdr_line + 0.05 * y_limit,
      label = paste0("BH padj = ", signif(volcano_fdr_cutoff, 3)),
      colour = "blue",
      hjust = -0.05,
      size = 4
    ) +
    ggplot2::geom_point(alpha = 0.8) +
    ggrepel::geom_label_repel(
      data = dat |> dplyr::filter(.data$label_me),
      ggplot2::aes(label = .data$label_text),
      size = 3,
      max.overlaps = 40,
      label.size = 0.15,
      min.segment.length = 0
    ) +
    ggplot2::facet_wrap(method ~ database, ncol = 3, scales = "free") +
    ggplot2::scale_y_continuous(limits = c(0, y_limit * 1.1)) +
    ggplot2::scale_size_continuous(
      name = "Gene-set size / count",
      range = c(2, 7),
      labels = scales::label_number(accuracy = 1)
    ) +
    ggplot2::scale_colour_manual(
      values = c(
        "higher in neural_bias" = "#D55E00",
        "higher in TLS" = "#0072B2",
        "neutral" = "grey60"
      )
    ) +
    ggplot2::labs(
      title = bquote(
        "Enrichment volcano-style summary - " * bold(.(celltype_level))
      ),
      subtitle = "Green guides: NES/signed effect = -1 and 1",
      x = "Effect direction: NES for GSEA; signed mean log2FC for ORA",
      y = expression(-log[10]("FDR-adjusted p-value")),
      colour = "Direction"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        colour = "grey70"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = if (show_legend) "bottom" else "none"
    )
}

volcano_pages <- lapply(
  names(celltype_results),
  plot_volcano_one_celltype,
  show_legend = FALSE
) |>
  purrr::compact()

if (length(volcano_pages) > 0L) {
  volcano_legend_source <- plot_volcano_one_celltype(
    names(celltype_results)[[1]],
    show_legend = TRUE
  )
  volcano_legend <- cowplot::get_legend(volcano_legend_source)
  volcano_pages <- lapply(volcano_pages, function(page) {
    cowplot::plot_grid(
      page,
      volcano_legend,
      ncol = 1,
      rel_heights = c(1, 0.08)
    )
  })

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(output_prefix, "_06_volcano_ORA_GSEA.pdf")
    ),
    plot = gridExtra::marrangeGrob(volcano_pages, nrow = 1, ncol = 1),
    width = 28,
    height = 18,
    units = "in",
    limitsize = FALSE
  )
}

# ============================================================================
# 8. Save compact enrichment output ----
# ============================================================================

saveRDS(
  list(
    parameters = list(
      study = study,
      technique = technique,
      input_de_file = input_de_file,
      fdr_cutoff = fdr_cutoff,
      volcano_fdr_cutoff = volcano_fdr_cutoff,
      lfc_cutoff_ora = lfc_cutoff_ora,
      msig_collection = msig_collection,
      msig_subcollection = msig_subcollection,
      go_ontologies = go_ontologies,
      min_gs_size = min_gs_size,
      max_gs_size = max_gs_size
    ),
    summary_tbl = summary_tbl,
    enrichment_long = enrichment_long,
    enrichment_results = celltype_results
  ),
  results_rds_file
)

message("clusterProfiler enrichment workflow completed.")

# Minimal Manhattan plot reconstruction from the saved compact RDS.
if (FALSE) {
  enrichment_test <- readRDS(
    # "outputs/biological-exploration/enrichment-analyses/GSE250136_pseudobulk_deseq2_clusterprofiler_enrichment_results.rds"
    file.path(
      output_dir,
      paste0(output_prefix, "_clusterprofiler_enrichment_results.rds")
    )
  )

  celltype_for_test <- "neuromesodermal.progenitors"
  ora_objs <- collect_method_objects(
    enrichment_test$enrichment_results[[celltype_for_test]]$ORA,
    method_label = "ORA"
  )
  gsea_objs <- collect_method_objects(
    enrichment_test$enrichment_results[[celltype_for_test]]$GSEA,
    method_label = "GSEA"
  )

  class(ora_objs$GO_BP)
  class(gsea_objs$GO_BP)

  patchwork::wrap_plots(
    list(
      ORA = enrichplot::manhattanplot(
        ora_objs,
        showCategory = 10,
        split = "database"
      ),
      GSEA = enrichplot::manhattanplot(
        gsea_objs,
        showCategory = 10,
        split = "database"
      )
    ),
    ncol = 2
  ) +
    ggplot2::geom_hline(
      yintercept = -log10(enrichment_test$parameters$fdr_cutoff),
      colour = "blue",
      linetype = "dashed",
      linewidth = 0.6
    )
}

# temp_enrichment <- readRDS(
#   "outputs/biological-exploration/enrichment-analyses/GSE250136_mast_clusterprofiler_enrichment_results.rds"
# )
# temp_enrichment_long <- temp_enrichment$enrichment_long

# readr::write_csv(
#   temp_enrichment_long,
#   "temp_GSE250136_mast_clusterprofiler.csv"
# )
