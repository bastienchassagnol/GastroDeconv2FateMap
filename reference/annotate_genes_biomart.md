# Annotate mouse genes with Ensembl metadata

Query Ensembl BioMart for mouse genes and return a collapsed,
one-row-per-gene summary including Ensembl IDs, description(s),
estimated gene length, and genomic coordinates. If multiple Ensembl
records map to one gene symbol, the function concatenates text fields
and averages gene lengths.

## Usage

``` r
annotate_genes_biomart(
  genes,
  drop_missing_genes = TRUE,
  verbose = TRUE,
  batch_size = 200L
)
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

- batch_size:

  Maximum number of gene symbols per BioMart request. Large inputs are
  split automatically to reduce timeouts.

## Value

A `data.frame` with one row per input gene and columns: `gene_name`,
`ensembl_id`, `biological_function`, `gene_length`, `chromosome_name`,
`gene_biotype`, and `ambiguous_genes`.

## Details

Output structure matches
[`annotate_genes_local()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/annotate_genes_local.md).

## Examples

``` r
if (FALSE) { # \dontrun{
annotate_genes_biomart(c("Gnai3", "H19", "Scml2"))
} # }
```
