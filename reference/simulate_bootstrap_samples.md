# Simulate pseudo-bulk samples by cell resampling (bootstrap)

Strategies 2 and 3 from `vignettes/pseudo_bulk_generation.qmd`.

- `replicate_type = "biological"`:

  Phenotype-stratified bootstrap: resample cells from all organoids
  sharing the target phenotype (strategy 2).

- `replicate_type = "technical"`:

  Organoid-preserving bootstrap: resample cells within a single barcode
  before aggregation (strategy 3).

## Usage

``` r
simulate_bootstrap_samples(
  seurat_obj,
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  n_samples = 100L,
  replicate_type = c("biological", "technical"),
  cells_per_sample = NULL,
  assay = "RNA",
  seed = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object with raw counts.

- phenotype_col:

  Metadata column for phenotype labels.

- barcode_col:

  Metadata column for barcode identifiers.

- n_samples:

  Total number of pseudo-bulk samples to generate. Phenotypes are
  allocated as evenly as possible across samples.

- replicate_type:

  `"biological"` (phenotype-level pooling) or `"technical"`
  (within-barcode resampling).

- cells_per_sample:

  Number of cells aggregated per pseudo-bulk sample. Defaults to the
  median barcode cell count in the object.

- assay:

  Assay from which raw counts are read.

- seed:

  Optional random seed for reproducibility.

## Value

A SummarizedExperiment with raw summed counts and balanced phenotypes in
`colData`.

## See also

[`aggregate_barcode_pseudo_bulk()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/aggregate_barcode_pseudo_bulk.md),
[`simulate_generative_models()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/simulate_generative_models.md)
