##' Annotate mouse genes using local Bioconductor databases
##'
##' Offline-first annotation for mouse gene symbols by combining:
##' - `org.Mm.eg.db` for ID mapping (SYMBOL, ENSEMBL, ENTREZID, GENENAME, GO)
##' - `TxDb.Mmusculus.UCSC.mm10.knownGene` for genomic coordinates and gene length
##' - `GO.db` for GO term labels
##'
##' The function returns one row per input gene symbol, collapsing multi-mapping
##' records and averaging gene length when multiple loci/transcripts are linked.
##'
##' @param genes Character vector of mouse gene symbols
##'   (for example, `c("Gnai3", "H19")`).
##' @param verbose Logical; whether to print progress messages.
##'
##' @return A `data.frame` with one row per input gene and columns:
##'   `gene_name`, `ensembl_id`, `biological_function`, `gene_length`,
##'   `chromosome_name`, `gene_biotype`, and `ambiguous_genes`.
##'
##' @examples
##' \dontrun{
##' annotate_genes(c("Gnai3", "Pbsn", "Cdc45"))
##' }
annotate_genes_mouse <- function(genes, verbose = TRUE) {
  # deal with dependencies
  required_pkgs <- c(
    "AnnotationDbi",
    "GenomicFeatures",
    "org.Mm.eg.db",
    "TxDb.Mmusculus.UCSC.mm10.knownGene"
  )

  missing_pkgs <- required_pkgs[
    !vapply(
      required_pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing_pkgs, collapse = ", "),
      ". Install with BiocManager::install()."
    )
  }

  # check input
  if (!is.character(genes) || length(genes) == 0) {
    stop("`genes` must be a non-empty character vector.")
  }

  genes <- unique(genes)

  if (verbose) {
    message("Querying local mouse annotation databases...")
  }

  # get correspondance between gene symbols and other gene identifiers
  id_map <- AnnotationDbi::select(
    x = org.Mm.eg.db::org.Mm.eg.db,
    keys = genes,
    columns = c("SYMBOL", "ENSEMBL", "ENTREZID", "GENENAME", "GENETYPE"),
    keytype = "SYMBOL"
  ) |>
    dplyr::rename(
      gene_name = SYMBOL,
      ensembl_id = ENSEMBL,
      entrez_id = ENTREZID,
      gene_name_long = GENENAME,
      gene_type_raw = GENETYPE
    )
  # Get genomic coordinates and gene length
  gene_ranges <- GenomicFeatures::genes(
    TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
    single.strand.genes.only = TRUE
  )
  chromosome_names <- data.frame(
    entrez_id = names(gene_ranges),
    chromosome_name = as.character(GenomicRanges::seqnames(gene_ranges)),
    gene_length = as.numeric(BiocGenerics::width(gene_ranges)),
    stringsAsFactors = FALSE
  )
  # Annotate genes with genomic coordinates and gene length
  annotated_genes <- id_map |>
    dplyr::distinct(
      gene_name,
      ensembl_id,
      entrez_id,
      gene_name_long,
      .keep_all = TRUE
    ) |>
    dplyr::left_join(chromosome_names, by = "entrez_id") |>
    dplyr::mutate(
      biological_function = gene_name_long,
      gene_biotype = dplyr::case_when(
        is.na(gene_type_raw) ~ NA_character_,
        gene_type_raw %in%
          c("protein-coding", "protein_coding") ~ "protein_coding",
        grepl("lnc|long non", gene_type_raw, ignore.case = TRUE) ~ "lncRNA",
        grepl("miRNA|microRNA", gene_type_raw, ignore.case = TRUE) ~ "miRNA",
        grepl("snRNA", gene_type_raw, ignore.case = TRUE) ~ "snRNA",
        grepl("snoRNA", gene_type_raw, ignore.case = TRUE) ~ "snoRNA",
        grepl("rRNA", gene_type_raw, ignore.case = TRUE) ~ "rRNA",
        grepl("tRNA", gene_type_raw, ignore.case = TRUE) ~ "tRNA",
        grepl(
          "ncRNA|non.?coding",
          gene_type_raw,
          ignore.case = TRUE
        ) ~ "other_non_coding_RNA",
        TRUE ~ "other"
      )
    ) |>
    dplyr::group_by(gene_name) |>
    dplyr::summarise(
      entrez_id = paste(unique(stats::na.omit(entrez_id)), collapse = ";"),
      ensembl_id = paste(unique(stats::na.omit(ensembl_id)), collapse = ";"),
      biological_function = paste(
        unique(stats::na.omit(biological_function)),
        collapse = ";"
      ),
      gene_length = mean(gene_length, na.rm = TRUE),
      chromosome_name = paste(
        unique(stats::na.omit(chromosome_name)),
        collapse = ";"
      ),
      gene_biotype = paste(
        unique(stats::na.omit(gene_biotype)),
        collapse = ";"
      ),
      n_ensembl_id = dplyr::n_distinct(stats::na.omit(ensembl_id)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      ambiguous_genes = n_ensembl_id > 1
    ) |>
    dplyr::select(
      gene_name,
      ensembl_id,
      biological_function,
      gene_length,
      chromosome_name,
      gene_biotype,
      ambiguous_genes
    )

  missing_genes <- setdiff(genes, annotated_genes$gene_name)
  if (length(missing_genes) > 0) {
    missing_df <- tibble::tibble(
      gene_name = missing_genes,
      ensembl_id = NA_character_,
      biological_function = NA_character_,
      gene_length = NA_real_,
      chromosome_name = NA_character_,
      gene_biotype = NA_character_,
      ambiguous_genes = FALSE
    )
    annotated_genes <- dplyr::bind_rows(annotated_genes, missing_df)
  }

  return(annotated_genes)
}
