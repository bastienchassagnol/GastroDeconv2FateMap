# Simulate log-normal single-cell counts

Draws `n_cells` count vectors from per-gene log-normal models,
optionally scaled by a barcode library factor.

## Usage

``` r
.simulate_lognormal_cells(bio_params, tech_scale, n_cells, seed = NULL)
```

## Arguments

- bio_params:

  List from
  [`.infer_lognormal_gene_params()`](https://bastienchassagnol.github.io/GastroDeconv2FateMap/reference/dot-infer_lognormal_gene_params.md).

- tech_scale:

  Barcode library-size scaling factor.

- n_cells:

  Number of cells to simulate.

- seed:

  Optional random seed.

## Value

Numeric matrix (genes x cells).

## Details

For each cell \\c\\ and gene \\g\\: \$\$z\_{gc} \sim
\mathcal{N}(\hat{\mu}\_g, \hat{\sigma}\_g^2), \quad x\_{gc} =
\max\left\\0, \left\lfloor e^{z\_{gc}} \cdot s_b - 1 \right\rfloor
\right\\\$\$
