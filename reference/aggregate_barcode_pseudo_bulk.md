# Aggregate pseudo-bulk samples by barcode (organoid)

Strategy 1 from `vignettes/pseudo_bulk_generation.qmd`: treat each
barcode as one bulk sample by summing raw single-cell counts within
barcodes, analogous to
[`AggregateExpression`](https://satijalab.org/seurat/reference/AggregateExpression.html)
with `slot = "counts"` and no normalisation.

## Usage

``` r
aggregate_barcode_pseudo_bulk(
  seurat_obj,
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  assay = "RNA",
  cell_mask = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object with raw counts in the RNA assay.

- phenotype_col:

  Metadata column for organoid-level phenotype labels (e.g.
  `"Morphotype"`).

- barcode_col:

  Metadata column for barcode / organoid identifiers (e.g.
  `"Sample.barcode"`).

- assay:

  Assay name from which raw counts are read. Default `"RNA"`.

- cell_mask:

  Optional logical vector aligned with `seurat_obj` cells; when
  provided, only these cells are aggregated.

## Value

A SummarizedExperiment with one column per distinct barcode. `colData`
contains `barcode_id`, the phenotype label, `library_depth` (column sum
of counts), and `simulation_method = "barcode_aggregation"`.

## See also

[`simulate_bootstrap_samples()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/simulate_bootstrap_samples.md),
[`simulate_generative_models()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/simulate_generative_models.md),
`vignettes/pseudo_bulk_generation.qmd`
