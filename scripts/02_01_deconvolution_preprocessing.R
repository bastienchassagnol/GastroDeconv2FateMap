# ==========================================================================
# 0. Libraries and Filename settings ----
# ==========================================================================

library(Seurat)
library(patchwork)
library(ggplot2)

# useful scripts
source("./R/01_03_gene_feature_annotation.R")

study <- "suppinger"
today <- format(Sys.Date(), "%Y-%m-%d")
output_dir <- "data/intermediate"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ==========================================================================
# 1. Load single-cell and bulk objects ----
# ==========================================================================

suppinger_single_cell_seurat <- readRDS(
  "./data/raw/GSE229513_gastruloidsobject.rds"
)

bulk_rnaseq_htseq <- read.table(
  file = "./data/raw/GSE229386_AllHTSeqCountsWithGeneNames.txt.gz",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

colnames(bulk_rnaseq_htseq)
dim(bulk_rnaseq_htseq)


# ==========================================================================
# 2. Compare gene names between single-cell and bulk RNA-seq data
# ==========================================================================

# --- 2a. Get gene names from single-cell data ----
single_cell_genes <- rownames(
  suppinger_single_cell_seurat@assays$RNA@meta.features
)

bulk_rnaseq_genes <- bulk_rnaseq_htseq$gene_name


# --- 2b. venn diagram of gene names --------------------------------------

venn_file <- file.path(
  output_dir,
  paste0(study, "_venn_gene_names_single_vs_bulk_", today, ".pdf")
)

gene_sets <- list(
  "Single-cell RNA-seq" = single_cell_genes,
  "Bulk RNA-seq" = bulk_rnaseq_genes
)

venn_plot <- ggVennDiagram::ggVennDiagram(
  gene_sets,
  set_size = 4
) +
  ggplot2::scale_fill_gradient(low = "grey90", high = "red") +
  ggplot2::ggtitle("Gene name overlap: single-cell vs bulk") +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0.22, 0.08))
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    plot.margin = ggplot2::margin(t = 10, r = 20, b = 10, l = 35)
  )

ggplot2::ggsave(
  filename = venn_file,
  plot = venn_plot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)

# --- 2c. understand what went wrong with the gene names --------------------------------------

shared_genes <- intersect(single_cell_genes, bulk_rnaseq_genes)

rnaseq_features_filtered <- tibble::tibble(
  gene_name = bulk_rnaseq_genes,
  ensembl_id = rownames(bulk_rnaseq_htseq)
) |>
  dplyr::filter(gene_name %in% shared_genes) |>
  # Optional but often useful: remove Ensembl version suffix (e.g. ".16")
  dplyr::mutate(ensembl_id = stringr::str_remove(ensembl_id, "\\..*$"))

print(paste(
  "Number of genes before filtering: ",
  nrow(rnaseq_features_filtered)
))

ambiguous_genes <- rnaseq_features_filtered |>
  dplyr::distinct(gene_name, ensembl_id) |>
  dplyr::count(gene_name, name = "n_ensembl_ids") |>
  dplyr::filter(n_ensembl_ids > 1) |>
  dplyr::pull(gene_name)
# All gene_name / ensembl_id pairs involved in non-1:1 mappings
ambiguous_gene_pairs <- rnaseq_features_filtered |>
  dplyr::filter(gene_name %in% ambiguous_genes) |>
  dplyr::arrange(gene_name, ensembl_id)
print(ambiguous_gene_pairs)

rnaseq_features_filtered <- rnaseq_features_filtered |>
  dplyr::filter(!gene_name %in% ambiguous_genes)

print(paste(
  "Number of genes after filtering: ",
  nrow(rnaseq_features_filtered)
))


test_annotation <- annotate_genes(
  genes = c("Gnai3", "Pbsn", "Cdc45", "H19", "Scml2", "Apoh"),
  dataset = "mmusculus_gene_ensembl",
  mirror = "www",
  verbose = TRUE
)



  mart <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = "mmusculus_gene_ensembl"
  )

  # Attributes to retrieve
  attributes <- c(
    "external_gene_name", # gene symbol
    "ensembl_gene_id", # Ensembl ID
    "description", # biological function
    "gene_biotype", # gene type
    "chromosome_name",
    "start_position",
    "end_position"
  )

  # Query biomart
  if (verbose) {
    message("Querying biomart...")
  }
  bm <- biomaRt::getBM(
    attributes = attributes,
    filters = "external_gene_name",
    values = c("Gnai3", "Pbsn"),
    mart = mart
  )


# ==========================================================================
# 3. Build SummarizedExperiment for the bulk RNA-seq data ----
# ==========================================================================

# --- 3a. Count matrix ----------------------------------------------------
# gene_name is a metadata column; keep Ensembl IDs (row.names) as row identifiers
rownames(bulk_rnaseq_htseq) <- NULL
count_matrix <- bulk_rnaseq_htseq |>
  dplyr::filter(gene_name %in% rnaseq_features_filtered$gene_name) |>
  tibble::column_to_rownames(var = "gene_name") |>
  as.matrix()
storage.mode(count_matrix) <- "integer"

# --- 3b. Row metadata (feature data) -------------------------------------
row_data <- S4Vectors::DataFrame(
  ensembl_id = rnaseq_features_filtered$ensembl_id,
  row.names = rnaseq_features_filtered$gene_name
)

# --- 3c. Column metadata (sample phenotype data) -------------------------
sample_names <- colnames(count_matrix)

time_point_id <- stringr::str_extract(sample_names, "\\d+h") |>
  as.factor()
treatment_status <- stringr::str_extract(
  sample_names,
  "(?<=_)(CTL|Early)(?=_)"
) |>
  forcats::fct_recode(control = "CTL", early_treatment = "Early") |>
  forcats::fct_relevel("control", "early_treatment")

batch_id <- stringr::str_extract(sample_names, "rep\\d+") |>
  as.factor()

col_data <- S4Vectors::DataFrame(
  time_point_id = time_point_id,
  batch_id = batch_id,
  treatment_status = treatment_status,
  row.names = sample_names
)

tinytable::tt(
  col_data |>
    as.data.frame(),
  caption = "Sample phenotype data (bulk RNA-seq)"
)

# --- 3d. Global metadata from GEO: GSE229386 -----------------------------
geo_metadata <- list(
  geo_accession = "GSE229386",
  title = "Multimodal characterization of murine gastruloid development",
  organism = "Mus musculus",
  experiment_type = "Expression profiling by high throughput sequencing",
  overall_design = "RNAseq from different timing of Wnt Pulse",
  platform = "GPL21103 — Illumina HiSeq 4000 (Mus musculus)",
  submission_date = "Apr 11, 2023",
  last_update = "Sep 07, 2023",
  contact = "Alexandre Mayran <alexandre.mayran@epfl.ch>, EPFL, Lausanne, Switzerland",
  contributors = c("Mayran A", "Suppinger S", "Aizarani N", "Liberali P"),
  citation = paste(
    "Suppinger S et al. Multimodal characterization of murine gastruloid development.",
    "Cell Stem Cell 2023 Jun 1;30(6):867-884.e11. PMID: 37209681"
  ),
  bioproject = "PRJNA954338",
  geo_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229386",
  summary = paste(
    "Gastruloids are highly scalable, three-dimensional assemblies generated from",
    "pluripotent stem cells that recapitulate fundamental principles of embryonic",
    "pattern formation in vitro. Using single cell RNA and multiome sequencing,",
    "this study provides a comprehensive resource mapping cellular states and cell",
    "types found during gastruloid development and compares them to the in vivo embryo."
  )
)

# --- 3e. Assemble SummarizedExperiment -----------------------------------
bulk_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = count_matrix),
  rowData = row_data,
  colData = col_data,
  metadata = geo_metadata
)

saveRDS(
  bulk_se,
  file = paste0(
    "./data/intermediate/",
    study,
    "_bulk_summarized_experiment_",
    today,
    ".rds"
  )
)

# ==========================================================================
# 4. Build SingleCellExperiment for the single-cell data ----
# ==========================================================================

# --- 4a. Slim the Seurat object to shared genes and minimal slots ---------
# Set RNA as active assay before any subsetting so row operations target counts
SeuratObject::DefaultAssay(suppinger_single_cell_seurat) <- "RNA"

sc_filtered <- subset(
  suppinger_single_cell_seurat,
  features = rnaseq_features_filtered$gene_name
)

sc_filtered@meta.data$timepoints <- factor(
  sc_filtered@meta.data$timepoints,
  levels = c(
    "0h",
    "24h",
    "36h",
    "48h",
    "52h",
    "56h",
    "60h",
    "72h",
    "84h",
    "96h",
    "108h",
    "120h"
  ),
  ordered = TRUE
)

tinytable::tt(
  sc_filtered@meta.data |>
    dplyr::select(
      "celltypeannotation",
      "timepoints",
      "batch",
      "nCount_RNA",
      "nFeature_RNA",
    ) |>
    head(),
  caption = "Cell-level metadata (single-cell RNA-seq)"
)

skimr::skim(sc_filtered@meta.data)

tinytable::tt(
  sc_filtered@meta.data |>
    head(),
  caption = "Cell-level metadata (single-cell RNA-seq)"
)


tinytable::tt(
  sc_filtered@meta.features |>
    head(),
  caption = "Feature metadata (single-cell RNA-seq)"
)

tinytable::tt(
  sc_filtered@meta.data |>
    dplyr::distinct(
      celltypeannotation,
      timepoints
    ) |>
    dplyr::arrange(timepoints, celltypeannotation),
  caption = "Cell types per time point"
)

# --- 4b. Convert to SingleCellExperiment ----------------------------------
suppinger_single_cell <- Seurat::as.SingleCellExperiment(
  sc_filtered,
  assay = "RNA"
)

suppinger_single_cell

# --- 4c. Save -------------------------------------------------------------
saveRDS(
  suppinger_single_cell,
  file = paste0(
    "./data/intermediate/",
    study,
    "_single_cell_",
    today,
    ".rds"
  )
)
