# Annotate mouse genes using local Bioconductor databases

Offline-first annotation for mouse gene symbols by combining:

- `org.Mm.eg.db` for ID mapping (SYMBOL, ENSEMBL, ENTREZID, GENENAME,
  GO)

- `TxDb.Mmusculus.UCSC.mm10.knownGene` for genomic coordinates and gene
  length

- `GO.db` for GO term labels

## Usage

``` r
annotate_genes_local(genes, drop_missing_genes = TRUE, verbose = TRUE)
```

## Arguments

- genes:

  Character vector of mouse gene symbols (for example,
  `c("Gnai3", "H19")`).

- drop_missing_genes:

  Logical; whether to drop genes that are not found in the annotation
  database. Default is `TRUE`.

- verbose:

  Logical; whether to print progress messages.

## Value

A `data.frame` with one row per input gene and columns: `gene_name`,
`ensembl_id`, `biological_function`, `gene_length`, `chromosome_name`,
`gene_biotype`, and `ambiguous_genes`.

## Details

The function returns one row per input gene symbol, collapsing
multi-mapping records and averaging gene length when multiple
loci/transcripts are linked.

## Examples

``` r
if (FALSE) { # \dontrun{
annotate_genes_local(c("Gnai3", "Pbsn", "Cdc45"))
} # }
```
