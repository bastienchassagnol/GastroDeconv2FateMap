# ============================================================================
# 0. Library and hyperparameter Setup ----
# ============================================================================

study <- "GSE250136"
technique <- "clusterprofiler_enrichment"
today <- format(Sys.Date(), "%Y-%m-%d")
input_de_file <- paste0(
  "outputs/biological-exploration/DEA-analyses/",
  "GSE250136_pseudobulk_deseq2_120h_celltype_model_results_2026-07-08.csv"
)
output_dir <- "outputs/biological-exploration/enrichment-analyses"
payload_file <- file.path(
  output_dir,
  paste0(study, "_clusterprofiler_enrichment_payload_", today, ".rds")
)
use_cached_payload <- file.exists(payload_file)

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

df_DEGs <- readr::read_csv(input_de_file, show_col_types = FALSE)

# ============================================================================
# 1. Reshape DESeq2 wide output to one row per gene and cell type ----
# ============================================================================

fc_cols <- names(df_DEGs) |>
  stringr::str_subset("__log2FoldChange$")
celltypes <- fc_cols |>
  stringr::str_remove("__log2FoldChange$")

de_long <- purrr::map_dfr(celltypes, function(celltype_level) {
  tibble::tibble(
    study = df_DEGs$study,
    de_technique = df_DEGs$technique,
    analysis_level = df_DEGs$analysis_level,
    gene_symbol = df_DEGs$gene,
    celltype = celltype_level,
    log2_fc = df_DEGs[[paste0(celltype_level, "__log2FoldChange")]],
    pvalue = df_DEGs[[paste0(celltype_level, "__pvalue")]],
    padj = df_DEGs[[paste0(celltype_level, "__padj_celltype_BH")]]
  )
}) |>
  dplyr::mutate(
    pvalue = tidyr::replace_na(pvalue, 1),
    padj = tidyr::replace_na(padj, 1),
    log2_fc = tidyr::replace_na(log2_fc, 0),
    p_for_rank = pmax(pvalue, .Machine$double.xmin),
    rank_score = log2_fc * (-log10(p_for_rank)),
    direction = dplyr::case_when(
      log2_fc > 0 ~ "up_in_neural_bias",
      log2_fc < 0 ~ "down_in_neural_bias",
      TRUE ~ "zero"
    ),
    is_de_for_ora = padj < fdr_cutoff & abs(log2_fc) >= lfc_cutoff_ora
  )

# ============================================================================
# 2. Map mouse gene symbols to Entrez IDs ----
# ============================================================================

# clusterProfiler, ReactomePA, and MSigDB Reactome sets are most stable here
# when keyed by Entrez IDs. The DESeq2 table stores mouse gene symbols, so we
# map SYMBOL -> ENTREZID with the Bioconductor mouse OrgDb (`org.Mm.eg.db`).
# Unmapped symbols are dropped from enrichment because they cannot be placed in
# the tested-gene universe. One-to-many mappings are retained at this stage and
# later collapsed for GSEA by keeping the largest absolute ranking score per
# Entrez ID, while ORA uses distinct Entrez IDs.
gene_map <- clusterProfiler::bitr(
  unique(de_long$gene_symbol),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db::org.Mm.eg.db
) |>
  tibble::as_tibble() |>
  dplyr::distinct(SYMBOL, ENTREZID) |>
  dplyr::mutate(ENTREZID = as.character(ENTREZID))

de_annot <- de_long |>
  dplyr::inner_join(gene_map, by = c("gene_symbol" = "SYMBOL")) |>
  dplyr::filter(!is.na(ENTREZID), is.finite(rank_score))

message(
  "Mapped genes: ",
  dplyr::n_distinct(de_annot$ENTREZID),
  " / ",
  dplyr::n_distinct(de_long$gene_symbol)
)

# ============================================================================
# 3. Prepare MSigDB Reactome TERM2GENE and TERM2NAME tables ----
# ============================================================================

# Retrieve the mouse-native MSigDB curated Reactome collection only.
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
    dplyr::filter(!is.na(gene), gene != "") |>
    dplyr::distinct(),
  term2name = msig_raw |>
    dplyr::transmute(
      term = .data$gs_name,
      name = dplyr::coalesce(.data$gs_description, .data$gs_name)
    ) |>
    dplyr::distinct()
)

# ============================================================================
# 4. Helpers for ORA, GSEA, result tidying, and plotting ----
# ============================================================================

# Return the tested, mappable Entrez universe for one cell type.
make_universe <- function(df_celltype) {
  df_celltype |>
    dplyr::distinct(ENTREZID) |>
    dplyr::pull(ENTREZID) |>
    as.character()
}

# Return ORA input genes using FDR and log2FC thresholds.
make_ora_genes <- function(df_celltype) {
  df_celltype |>
    dplyr::filter(is_de_for_ora) |>
    dplyr::distinct(ENTREZID) |>
    dplyr::pull(ENTREZID) |>
    as.character()
}

# Build the decreasing named ranking vector required by clusterProfiler GSEA.
make_gsea_gene_list <- function(df_celltype) {
  rank_tbl <- df_celltype |>
    dplyr::filter(is.finite(rank_score), !is.na(ENTREZID)) |>
    dplyr::group_by(ENTREZID) |>
    dplyr::slice_max(
      order_by = abs(rank_score),
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      rank_score = rank_score + dplyr::dense_rank(log2_fc) * 1e-12
    ) |>
    dplyr::arrange(dplyr::desc(rank_score))

  gene_list <- rank_tbl$rank_score
  names(gene_list) <- rank_tbl$ENTREZID
  sort(gene_list, decreasing = TRUE)
}

# Check whether an enrichment object contains at least one term.
has_terms <- function(x) {
  !is.null(x) &&
    inherits(x, c("enrichResult", "gseaResult")) &&
    nrow(as.data.frame(x)) > 0
}

# Convert one enrichment object to a tidy table with provenance columns.
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
      neglog10_padj = -log10(pmax(p.adjust, .Machine$double.xmin)),
      sig_fdr = p.adjust < fdr_cutoff,
      .before = 1
    )
}

# Add term similarities for network/tree visualisations when supported.
add_termsim <- function(obj, database) {
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

# ============================================================================
# 5. Run ORA and GSEA per cell type ----
# ============================================================================

run_enrichment_one_celltype <- function(celltype_level) {
  message("Running enrichment for: ", celltype_level)

  df_celltype <- de_annot |>
    dplyr::filter(celltype == celltype_level)

  universe <- make_universe(df_celltype)
  ora_genes <- make_ora_genes(df_celltype)
  gene_list <- make_gsea_gene_list(df_celltype)

  res <- list(ORA = list(), GSEA = list())

  for (ont in go_ontologies) {
    db_name <- paste0("GO_", ont)

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

if (use_cached_payload) {
  message("Loading existing enrichment payload: ", payload_file)
  payload <- readRDS(payload_file)
  de_long <- payload$de_long
  de_annot <- payload$de_annot
  celltype_results <- payload$celltype_results
  enrichment_long <- payload$enrichment_long
  sig_enrichment <- payload$sig_enrichment
  summary_tbl <- payload$summary_tbl
} else {
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

  sig_enrichment <- enrichment_long |>
    dplyr::filter(p.adjust < fdr_cutoff) |>
    dplyr::arrange(celltype, method, database, p.adjust)

  summary_tbl <- enrichment_long |>
    dplyr::group_by(celltype, method, database) |>
    dplyr::summarise(
      n_terms_tested = dplyr::n(),
      n_terms_FDR = sum(p.adjust < fdr_cutoff, na.rm = TRUE),
      top_term = Description[which.min(p.adjust)][1],
      top_FDR = min(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(celltype, method, database)
}

# ============================================================================
# 6. Native enrichplot visualisations saved with ggsave() ----
# ============================================================================

# Build a faceted-by-source page for one cell type and one visualisation type.
build_native_page <- function(
  celltype_level,
  plot_fun,
  title,
  needs_termsim = FALSE,
  methods = c("ORA", "GSEA")
) {
  res_ct <- celltype_results[[celltype_level]]
  db_order <- unique(c(names(res_ct$ORA), names(res_ct$GSEA)))
  plot_list <- list()

  for (method in methods) {
    for (db in db_order) {
      obj <- res_ct[[method]][[db]]

      if (!has_terms(obj)) {
        message(
          "No terms for ",
          title,
          " / ",
          celltype_level,
          " / ",
          method,
          " / ",
          db
        )
        next
      }

      if (needs_termsim) {
        obj <- add_termsim(obj, db)
      }

      plot_list[[paste(method, db, sep = "_")]] <- tryCatch(
        plot_fun(obj, method, db) +
          ggplot2::ggtitle(paste(method, db, sep = " | ")) +
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
    }
  }

  plot_list <- purrr::compact(plot_list)
  if (length(plot_list) == 0L) {
    message(
      "Skipping ",
      title,
      " page for ",
      celltype_level,
      ": no plottable terms."
    )
    return(NULL)
  }

  patchwork::wrap_plots(plot_list, ncol = length(db_order)) +
    patchwork::plot_annotation(
      title = paste0(title, " - ", celltype_level),
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

# Save one multi-page PDF using ggsave() and gridExtra::marrangeGrob().
save_native_grid_pdf <- function(
  filename,
  plot_fun,
  title,
  needs_termsim = FALSE,
  methods = c("ORA", "GSEA"),
  width = 26,
  height = 16
) {
  pages <- lapply(
    names(celltype_results),
    build_native_page,
    plot_fun = plot_fun,
    title = title,
    needs_termsim = needs_termsim,
    methods = methods
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

save_native_grid_pdf(
  filename = paste0(study, "_01_dotplot_ORA_GSEA_", today, ".pdf"),
  title = "Dotplot",
  plot_fun = function(x, method, database) {
    p <- if ("ONTOLOGY" %in% colnames(as.data.frame(x))) {
      enrichplot::dotplot(
        x,
        showCategory = show_category,
        label_format = 50,
        split = "ONTOLOGY"
      ) +
        ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free")
    } else {
      enrichplot::dotplot(x, showCategory = show_category, label_format = 50)
    }
    p +
      ggplot2::scale_size_continuous(range = c(5, 13)) +
      ggplot2::scale_colour_continuous(
        labels = scales::label_number(accuracy = 0.001)
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

save_native_grid_pdf(
  filename = paste0(study, "_02_manhattanplot_ORA_GSEA_", today, ".pdf"),
  title = "Manhattan plot",
  plot_fun = function(x, method, database) {
    enrichplot::manhattanplot(x, showCategory = 10) +
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
        label = paste0("BH padj = ", fdr_cutoff),
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

save_native_grid_pdf(
  filename = paste0(study, "_03_treeplot_ORA_GSEA_", today, ".pdf"),
  title = "Tree plot",
  plot_fun = function(x, method, database) {
    enrichplot::treeplot(
      x,
      nCluster = 6,
      label_format = 45,
      fontsize = 4,
      fontsize_cladelab = 5
    )
  },
  needs_termsim = TRUE,
  width = 32,
  height = 22
)

save_native_grid_pdf(
  filename = paste0(study, "_04_semantic_space_ORA_GSEA_", today, ".pdf"),
  title = "Semantic-space plot",
  plot_fun = function(x, method, database) {
    if ("ssplot" %in% getNamespaceExports("enrichplot")) {
      enrichplot::ssplot(x, nCluster = 6)
    } else {
      enrichplot::emapplot(x, node_label = "group", layout = "kk")
    }
  },
  needs_termsim = TRUE,
  width = 34,
  height = 24
)

save_native_grid_pdf(
  filename = paste0(study, "_05_upsetplot_ORA_only_", today, ".pdf"),
  title = "UpSet plot - ORA only",
  plot_fun = function(x, method, database) {
    enrichplot::upsetplot(x)
  },
  methods = "ORA",
  width = 28,
  height = 18
)

# ============================================================================
# 7. Volcano-style enrichment plot ----
# ============================================================================

# Split the slash-separated contributing gene IDs stored in enrichment tables.
parse_gene_ids <- function(gene_id_string) {
  if (is.na(gene_id_string) || gene_id_string == "") {
    return(character())
  }
  unlist(strsplit(gene_id_string, "/", fixed = TRUE))
}

# Estimate the ORA x-axis as signed mean DESeq2 log2FC for term genes.
compute_ora_signed_mean_lfc <- function(gene_id_string, celltype_level) {
  ids <- parse_gene_ids(gene_id_string)
  df_celltype <- de_annot |>
    dplyr::filter(celltype == celltype_level)

  if (length(intersect(ids, df_celltype$gene_symbol)) > 0) {
    mean(df_celltype$log2_fc[df_celltype$gene_symbol %in% ids], na.rm = TRUE)
  } else {
    mean(df_celltype$log2_fc[df_celltype$ENTREZID %in% ids], na.rm = TRUE)
  }
}

volcano_tbl <- enrichment_long |>
  dplyr::group_by(celltype, method, database) |>
  dplyr::arrange(p.adjust, .by_group = TRUE) |>
  dplyr::slice_head(n = 150) |>
  dplyr::ungroup() |>
  dplyr::rowwise() |>
  dplyr::mutate(
    volcano_x = dplyr::case_when(
      method == "GSEA" ~ NES,
      method == "ORA" ~ compute_ora_signed_mean_lfc(geneID, celltype),
      TRUE ~ NA_real_
    ),
    direction_label = dplyr::case_when(
      volcano_x > 0 ~ "higher in neural_bias",
      volcano_x < 0 ~ "higher in TLS",
      TRUE ~ "neutral"
    ),
    size_for_plot = dplyr::case_when(
      method == "GSEA" ~ as.numeric(setSize),
      method == "ORA" ~ as.numeric(Count),
      TRUE ~ NA_real_
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::filter(is.finite(volcano_x), is.finite(neglog10_padj))

# Build one volcano page with a shared legend and explicit FDR 0.05 guide line.
plot_volcano_one_celltype <- function(celltype_level, show_legend = FALSE) {
  dat <- volcano_tbl |>
    dplyr::filter(celltype == celltype_level) |>
    dplyr::group_by(method, database) |>
    dplyr::arrange(p.adjust, .by_group = TRUE) |>
    dplyr::mutate(
      label_me = dplyr::row_number() <= 3 |
        p.adjust < volcano_fdr_cutoff,
      label_text = paste0(
        stringr::str_trunc(Description, 55),
        "\nFDR = ",
        signif(p.adjust, 3)
      )
    ) |>
    dplyr::ungroup()

  if (nrow(dat) == 0) {
    message("No enrichment results for volcano-style plot: ", celltype_level)
    return(NULL)
  }

  y_limit <- max(dat$neglog10_padj, na.rm = TRUE)
  y_limit <- max(y_limit, -log10(volcano_fdr_cutoff) + 0.5)
  fdr_line <- -log10(volcano_fdr_cutoff)

  ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = volcano_x,
      y = neglog10_padj,
      colour = direction_label,
      size = size_for_plot
    )
  ) +
    ggplot2::geom_hline(
      yintercept = fdr_line,
      colour = "blue",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::annotate(
      "text",
      x = -Inf,
      y = fdr_line + 0.05 * y_limit,
      label = "BH padj = 0.05",
      colour = "blue",
      hjust = -0.05,
      size = 4
    ) +
    ggplot2::geom_point(alpha = 0.8) +
    ggrepel::geom_label_repel(
      data = dat |> dplyr::filter(label_me),
      ggplot2::aes(label = label_text),
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
      title = paste0("Enrichment volcano-style summary - ", celltype_level),
      subtitle = "GSEA uses NES; ORA uses signed mean log2FC of term genes",
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
    cowplot::plot_grid(page, volcano_legend, ncol = 1, rel_heights = c(1, 0.08))
  })

  ggplot2::ggsave(
    filename = file.path(
      output_dir,
      paste0(study, "_06_volcano_ORA_GSEA_", today, ".pdf")
    ),
    plot = gridExtra::marrangeGrob(volcano_pages, nrow = 1, ncol = 1),
    width = 28,
    height = 18,
    units = "in",
    limitsize = FALSE
  )
}

# ============================================================================
# 8. Save a single serialised enrichment payload ----
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
      go_ontologies = go_ontologies
    ),
    de_long = de_long,
    de_annot = de_annot,
    celltype_results = celltype_results,
    enrichment_long = enrichment_long,
    sig_enrichment = sig_enrichment,
    summary_tbl = summary_tbl,
    volcano_tbl = volcano_tbl
  ),
  file.path(
    output_dir,
    paste0(study, "_clusterprofiler_enrichment_payload_", today, ".rds")
  )
)

message("clusterProfiler enrichment workflow completed.")
