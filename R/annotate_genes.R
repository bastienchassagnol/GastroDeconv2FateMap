#' Annotate mouse genes using local Bioconductor databases
#'
#' Offline-first annotation for mouse gene symbols by combining:
#' - `org.Mm.eg.db` for ID mapping (SYMBOL, ENSEMBL, ENTREZID, GENENAME, GO)
#' - `TxDb.Mmusculus.UCSC.mm10.knownGene` for genomic coordinates and gene length
#' - `GO.db` for GO term labels
#'
#' The function returns one row per input gene symbol, collapsing multi-mapping
#' records and averaging gene length when multiple loci/transcripts are linked.
#'
#' @param genes Character vector of mouse gene symbols
#'   (for example, `c("Gnai3", "H19")`).
#' @param drop_missing_genes Logical; whether to drop genes
#'   that are not found in the annotation database. Default is `TRUE`.
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A `data.frame` with one row per input gene and columns:
#'   `gene_name`, `ensembl_id`, `biological_function`, `gene_length`,
#'   `chromosome_name`, `gene_biotype`, and `ambiguous_genes`.
#'
#' @examples
#' \dontrun{
#' annotate_genes_local(c("Gnai3", "Pbsn", "Cdc45"))
#' }
#' @importFrom rlang .data
#' @export
annotate_genes_local <- function(
  genes,
  drop_missing_genes = TRUE,
  verbose = TRUE
) {
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
      gene_name = "SYMBOL",
      ensembl_id = "ENSEMBL",
      entrez_id = "ENTREZID",
      gene_name_long = "GENENAME",
      gene_type_raw = "GENETYPE"
    )
  # Get genomic coordinates and gene length
  gene_ranges <- suppressMessages(
    GenomicFeatures::genes(
      TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene,
    )
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
      "gene_name",
      "ensembl_id",
      "entrez_id",
      "gene_name_long",
      .keep_all = TRUE
    ) |>
    dplyr::left_join(chromosome_names, by = "entrez_id") |>
    dplyr::mutate(
      biological_function = .data[["gene_name_long"]],
      gene_biotype = dplyr::case_when(
        is.na(.data[["gene_type_raw"]]) ~ NA_character_,
        .data[["gene_type_raw"]] %in%
          c("protein-coding", "protein_coding") ~ "protein_coding",
        grepl("lnc|long non", .data[["gene_type_raw"]], ignore.case = TRUE) ~
          "lncRNA",
        grepl("miRNA|microRNA", .data[["gene_type_raw"]], ignore.case = TRUE) ~
          "miRNA",
        grepl("snRNA", .data[["gene_type_raw"]], ignore.case = TRUE) ~ "snRNA",
        grepl("snoRNA", .data[["gene_type_raw"]], ignore.case = TRUE) ~ "snoRNA",
        grepl("rRNA", .data[["gene_type_raw"]], ignore.case = TRUE) ~ "rRNA",
        grepl("tRNA", .data[["gene_type_raw"]], ignore.case = TRUE) ~ "tRNA",
        grepl(
          "ncRNA|non.?coding",
          .data[["gene_type_raw"]],
          ignore.case = TRUE
        ) ~ "other_non_coding_RNA",
        TRUE ~ "other"
      )
    ) |>
    dplyr::group_by(.data[["gene_name"]]) |>
    dplyr::summarise(
      entrez_id = paste(
        unique(stats::na.omit(.data[["entrez_id"]])),
        collapse = ";"
      ),
      ensembl_id = paste(
        unique(stats::na.omit(.data[["ensembl_id"]])),
        collapse = ";"
      ),
      biological_function = paste(
        unique(stats::na.omit(.data[["biological_function"]])),
        collapse = ";"
      ),
      gene_length = mean(.data[["gene_length"]], na.rm = TRUE),
      chromosome_name = paste(
        unique(stats::na.omit(.data[["chromosome_name"]])),
        collapse = ";"
      ),
      gene_biotype = paste(
        unique(stats::na.omit(.data[["gene_biotype"]])),
        collapse = ";"
      ),
      n_ensembl_id = dplyr::n_distinct(stats::na.omit(.data[["ensembl_id"]])),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      ambiguous_genes = .data[["n_ensembl_id"]] > 1
    ) |>
    dplyr::select(
      "gene_name",
      "ensembl_id",
      "biological_function",
      "gene_length",
      "chromosome_name",
      "gene_biotype",
      "ambiguous_genes"
    )

  # Convert empty strings, and inputs to NA (homogenize missing values)
  annotated_genes <- annotated_genes |>
    dplyr::mutate(
      dplyr::across(dplyr::where(is.character), \(x) dplyr::na_if(x, "")),
      dplyr::across(dplyr::where(is.numeric), \(x) ifelse(is.nan(x), NA_real_, x)),
    )

  # Drop genes in which not all biological annotation is available
  if (drop_missing_genes) {
    annotated_genes <- annotated_genes |>
      tidyr::drop_na()
  }

  return(annotated_genes)
}


#' Connect to Ensembl BioMart for mouse gene annotation
#'
#' Tries several connection strategies with retries. Regional mirrors
#' (`useast`, `asia`) often fail SSL checks; `host` and `version` routes
#' are preferred fallbacks when the default www mirror probe times out.
#'
#' @noRd
.connect_ensembl_biomart <- function(
  verbose = TRUE,
  max_attempts = 12L,
  retry_pause_sec = 5
) {
  biomart <- "genes"
  dataset <- "mmusculus_gene_ensembl"
  ensembl_host <- "https://www.ensembl.org"
  last_error <- NULL

  connect_strategies <- list(
    version_115 = function() {
      biomaRt::useEnsembl(
        biomart = biomart,
        dataset = dataset,
        version = 115,
        verbose = FALSE
      )
    },
    host = function() {
      biomaRt::useEnsembl(
        biomart = biomart,
        dataset = dataset,
        host = ensembl_host,
        verbose = FALSE
      )
    },
    mirror_www = function() {
      biomaRt::useEnsembl(
        biomart = biomart,
        dataset = dataset,
        mirror = "www",
        verbose = FALSE
      )
    }
  )

  for (attempt in seq_len(max_attempts)) {
    for (strategy_name in names(connect_strategies)) {
      if (verbose) {
        message(
          "Connecting to Ensembl BioMart (",
          strategy_name,
          ", attempt ",
          attempt,
          "/",
          max_attempts,
          ")..."
        )
      }

      mart <- tryCatch(
        connect_strategies[[strategy_name]](),
        error = function(err) {
          last_error <<- conditionMessage(err)
          NULL
        }
      )

      if (!is.null(mart)) {
        return(mart)
      }
    }

    if (attempt < max_attempts) {
      Sys.sleep(retry_pause_sec * attempt)
    }
  }

  stop(
    "Could not connect to Ensembl BioMart after ",
    max_attempts,
    " attempt(s). Last error: ",
    last_error,
    call. = FALSE
  )
}


#' Query Ensembl BioMart in batches
#'
#' @noRd
.get_biomart_batch <- function(
  mart,
  genes,
  attributes,
  max_attempts = 3L
) {
  last_error <- NULL
  for (attempt in seq_len(max_attempts)) {
    result <- tryCatch(
      biomaRt::getBM(
        attributes = attributes,
        filters = "external_gene_name",
        values = genes,
        mart = mart
      ),
      error = function(err) {
        last_error <<- conditionMessage(err)
        NULL
      }
    )
    if (!is.null(result)) {
      return(result)
    }
    Sys.sleep(2 * attempt)
  }

  if (length(genes) > 1L) {
    split_at <- ceiling(length(genes) / 2)
    return(dplyr::bind_rows(
      .get_biomart_batch(
        mart,
        genes[seq_len(split_at)],
        attributes,
        max_attempts
      ),
      .get_biomart_batch(
        mart,
        genes[seq(split_at + 1L, length(genes))],
        attributes,
        max_attempts
      )
    ))
  }

  stop(
    "BioMart query failed after ",
    max_attempts,
    " attempt(s). Last error: ",
    last_error,
    call. = FALSE
  )
}


.query_biomart_by_gene_batches <- function(
  mart,
  genes,
  attributes,
  batch_size = 200L,
  verbose = TRUE
) {
  if (length(genes) <= batch_size) {
    return(.get_biomart_batch(mart, genes, attributes))
  }

  batch_ids <- ceiling(seq_along(genes) / batch_size)
  batches <- split(genes, batch_ids)
  n_batches <- length(batches)

  bm_list <- lapply(seq_along(batches), function(i) {
    if (verbose) {
      message(
        "BioMart batch ",
        i,
        "/",
        n_batches,
        " (",
        length(batches[[i]]),
        " genes)..."
      )
    }
    .get_biomart_batch(mart, batches[[i]], attributes)
  })

  dplyr::bind_rows(bm_list)
}


#' Annotate mouse genes with Ensembl metadata
#'
#' Query Ensembl BioMart for mouse genes and return a collapsed, one-row-per-gene
#' summary including Ensembl IDs, description(s), estimated gene length, and
#' genomic coordinates. If multiple Ensembl records map to one gene symbol, the
#' function concatenates text fields and averages gene lengths.
#'
#' Output structure matches [annotate_genes_local()].
#'
#' @param genes Character vector of mouse gene symbols
#'   (for example, `c("Gnai3", "H19")`).
#' @param drop_missing_genes Logical; whether to drop genes
#'   that are not found in the annotation database. Default is `TRUE`.
#' @param verbose Logical; whether to print progress messages.
#' @param batch_size Maximum number of gene symbols per BioMart request.
#'   Large inputs are split automatically to reduce timeouts.
#'
#' @return A `data.frame` with one row per input gene and columns:
#'   `gene_name`, `ensembl_id`, `biological_function`, `gene_length`,
#'   `chromosome_name`, `gene_biotype`, and `ambiguous_genes`.
#'
#' @examples
#' \dontrun{
#' annotate_genes_biomart(c("Gnai3", "H19", "Scml2"))
#' }
#' @importFrom rlang .data
#' @export
annotate_genes_biomart <- function(
  genes,
  drop_missing_genes = TRUE,
  verbose = TRUE,
  batch_size = 200L
) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(
      "Package 'biomaRt' is required. Install with install.packages('biomaRt')"
    )
  }
  if (!is.character(genes) || length(genes) == 0) {
    stop("`genes` must be a non-empty character vector.")
  }
  genes <- unique(genes)

  mart <- .connect_ensembl_biomart(verbose = verbose)

  attributes <- c(
    "external_gene_name",
    "ensembl_gene_id",
    "description",
    "gene_biotype",
    "chromosome_name",
    "start_position",
    "end_position"
  )

  if (verbose) {
    message("Querying biomart...")
  }
  bm <- .query_biomart_by_gene_batches(
    mart = mart,
    genes = genes,
    attributes = attributes,
    batch_size = batch_size,
    verbose = verbose
  )

  annotated_genes <- bm |>
    dplyr::mutate(
      gene_length = .data[["end_position"]] - .data[["start_position"]] + 1
    ) |>
    dplyr::group_by(.data[["external_gene_name"]]) |>
    dplyr::summarise(
      ensembl_id = paste(
        unique(stats::na.omit(.data[["ensembl_gene_id"]])),
        collapse = ";"
      ),
      biological_function = paste(
        unique(stats::na.omit(.data[["description"]])),
        collapse = ";"
      ),
      gene_length = mean(.data[["gene_length"]], na.rm = TRUE),
      chromosome_name = paste(
        unique(stats::na.omit(.data[["chromosome_name"]])),
        collapse = ";"
      ),
      gene_biotype = paste(
        unique(stats::na.omit(.data[["gene_biotype"]])),
        collapse = ";"
      ),
      n_ensembl_id = dplyr::n_distinct(stats::na.omit(.data[["ensembl_gene_id"]])),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      ambiguous_genes = .data[["n_ensembl_id"]] > 1
    ) |>
    dplyr::rename(gene_name = "external_gene_name") |>
    dplyr::select(
      "gene_name",
      "ensembl_id",
      "biological_function",
      "gene_length",
      "chromosome_name",
      "gene_biotype",
      "ambiguous_genes"
    )

  # Output cleaning: Remove genes with partial annotation ----
  # This is to ensure that the output is a one-to-one mapping of gene symbols to annotations.
  annotated_genes <- annotated_genes |>
    dplyr::mutate(
      dplyr::across(dplyr::where(is.character), \(x) dplyr::na_if(x, "")),
      dplyr::across(dplyr::where(is.numeric), \(x) ifelse(is.nan(x), NA_real_, x)),
    )

  if (drop_missing_genes) {
    annotated_genes <- annotated_genes |>
      tidyr::drop_na()
  }

  return(annotated_genes)
}
