# Aggregate pseudo-bulk samples by cell type and barcode

Extension of
[`aggregate_barcode_pseudo_bulk()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/aggregate_barcode_pseudo_bulk.md):
sum raw single-cell counts within each cell type x barcode combination,
yielding one pseudo-bulk column per organoid within each annotated cell
type. Uses the same raw-count aggregation principle as Strategy 1 in
`vignettes/pseudo_bulk_generation.qmd`.

## Usage

``` r
aggregate_celltype_barcode_pseudo_bulk(
  seurat_obj,
  celltype_col = "luque_cluster_annotation",
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  assay = "RNA",
  cell_mask = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object with raw counts in the RNA assay.

- celltype_col:

  Metadata column for cell-type labels (e.g.
  `"luque_cluster_annotation"`).

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

A SummarizedExperiment with one column per distinct cell type x barcode
pair. `colData` contains `celltype`, `barcode_id`, the phenotype label,
`n_cells`, `library_depth`, and
`simulation_method = "celltype_barcode_aggregation"`.

## See also

[`aggregate_barcode_pseudo_bulk()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/aggregate_barcode_pseudo_bulk.md),
[`simulate_bootstrap_samples()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/simulate_bootstrap_samples.md),
[`simulate_generative_models()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/simulate_generative_models.md),
`vignettes/pseudo_bulk_generation.qmd`
