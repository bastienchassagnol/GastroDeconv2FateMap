##' Annotate mouse genes with Ensembl metadata
##'
##' Query Ensembl BioMart for mouse genes and return a collapsed, one-row-per-gene
##' summary including Ensembl IDs, description(s), estimated gene length, and
##' genomic position. If multiple Ensembl records map to one gene symbol, the
##' function concatenates text fields and averages gene lengths.
##'
##' @param genes Character vector of gene symbols (for example, `c("Gnai3", "H19")`).
##' @param dataset Ensembl BioMart dataset. Defaults to `"mmusculus_gene_ensembl"`.
##' @param mirror Ensembl mirror passed to [biomaRt::useEnsembl()].
##' @param verbose Logical; whether to print progress messages.
##'
##' @return A `data.frame` with one row per input gene and columns:
##'   `gene_name`, `ensembl_id`, `biological_function`, `gene_length`,
##'   `genome_position`, and `is_one_to_one`.
##'
##' @examples
##' \dontrun{
##' annotate_genes(c("Gnai3", "H19", "Scml2"))
##' }
annotate_genes <- function(
  genes,
  dataset = "mmusculus_gene_ensembl",
  mirror = "www",
  verbose = TRUE
) {
  # Required packages
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "Package 'biomaRt' is required. Install with install.packages('biomaRt')"
    )
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install with install.packages('dplyr')")
  }
  # Connect to Ensembl
  if (verbose) {
    message("Connecting to Ensembl biomart...")
  }
  mart <- biomaRt::useEnsembl(
    biomart = "genes",
    dataset = dataset,
    mirror = mirror
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
    values = genes,
    mart = mart
  )
  summarized_bm <- bm |> # Compute gene length from genomic coordinates
    dplyr::mutate(
      gene_length = end_position - start_position + 1
    ) |>
    dplyr::group_by(external_gene_name) |>
    dplyr::summarise(
      ensembl_id = paste(unique(ensembl_gene_id), collapse = ";"),
      chromosome_name = paste(unique(chromosome_name), collapse = ";"),
      biological_function = paste(unique(description), collapse = ";"),
      gene_length = mean(gene_length, na.rm = TRUE),
      n_ids = dplyr::n_distinct(ensembl_gene_id),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      ambiguous_genes = dplyr::if_else(n_ids > 1, TRUE, FALSE)
    ) |>
    dplyr::rename(gene_name = external_gene_name) |>
    dplyr::select(-n_ids)

  # Handle missing genes
  missing_genes <- setdiff(genes, summarized_bm$gene_name)
  if (length(missing_genes) > 0) {
    missing_bm <- data.frame(
      gene_name = missing_genes,
      ensembl_id = NA,
      biological_function = NA,
      gene_length = NA,
      genome_position = NA,
      is_one_to_one = NA,
      stringsAsFactors = FALSE
    )
    summarized_bm <- dplyr::bind_rows(summarized_bm, missing_df)
  }

  return(summarized_bm)
}
