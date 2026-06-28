# Simulate pseudo-bulk samples from generative count models

Splatter- and HADACA3-inspired generative simulation (strategies 7 and
10 in `vignettes/pseudo_bulk_generation.qmd`). Parameters are inferred
at the biological (phenotype) level; technical library scaling is
applied only for the log-normal model. Simulated cells are summed to
sample-level raw bulk counts.

## Usage

``` r
simulate_generative_models(
  seurat_obj,
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  n_samples = 100L,
  model = c("lognormal", "negative_binomial"),
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

  Total pseudo-bulk samples (balanced across phenotypes).

- model:

  `"lognormal"` for log-normal gene counts, or `"negative_binomial"` for
  Splatter-style negative-binomial counts.

- cells_per_sample:

  Cells simulated and summed per pseudo-bulk sample.

- assay:

  Assay from which raw counts are read.

- seed:

  Optional random seed.

## Value

A SummarizedExperiment with raw summed counts.

## Details

**Log-normal path**: per-gene \\(\hat{\mu}\_g, \hat{\sigma}\_g)\\ on
\\\log(1+x)\\ scale per phenotype, with barcode library scaling \\s_b\\
at simulation time.

**Negative-binomial path**: per-gene \\(\hat{\mu}\_g, \hat{\theta}\_g)\\
via method-of-moments (HADACA3-style) and per-cell library factors from
[`splatter::splatEstimate`](http://oshlacklab.com/splatter/reference/splatEstimate.md).
No extra barcode scaling is applied, because Splatter already models
library-size variation.

Bulk aggregation: \$\$Y\_{gi} = \sum\_{c \in \mathcal{C}\_i}
x\_{gc}^{\mathrm{sim}}\$\$

## See also

[`aggregate_barcode_pseudo_bulk()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/aggregate_barcode_pseudo_bulk.md),
[`simulate_bootstrap_samples()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/simulate_bootstrap_samples.md),
<https://oshlacklab.com/splatter/articles/splatter.html>
