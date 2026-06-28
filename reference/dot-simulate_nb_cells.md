# Simulate negative-binomial single-cell counts

Draws `n_cells` count vectors from per-gene NB models with cell-specific
library depth from Splatter estimates.

## Usage

``` r
.simulate_nb_cells(bio_params, n_cells, seed = NULL)
```

## Arguments

- bio_params:

  List with `gene` and `splatter` components.

- n_cells:

  Number of cells to simulate.

- seed:

  Optional random seed.

## Value

Numeric matrix (genes x cells).

## Details

For cell \\c\\ with library factor \\\ell_c\\ (from Splatter) and gene
\\g\\: \$\$x\_{gc} \sim \mathrm{NB}\\\left(\mu = \hat{\mu}\_g \cdot
\ell_c,\\ \mathrm{size} = \hat{\theta}\_g\right)\$\$
