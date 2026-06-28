# Draw per-cell library-size factors from Splatter estimates

Samples cell-level library-depth multipliers from the log-normal library
model inferred by
[`splatter::splatEstimate`](http://oshlacklab.com/splatter/reference/splatEstimate.md).

## Usage

``` r
.draw_splatter_library_factors(n_cells, splatter_params)
```

## Arguments

- n_cells:

  Number of cells to simulate.

- splatter_params:

  A `SplatParams` object.

## Value

Numeric vector of length `n_cells`.

## Details

With Splatter library parameters \\(\mathrm{loc}\_L,
\mathrm{scale}\_L)\\: \$\$\log L_c \sim \mathcal{N}(\mathrm{loc}\_L,
\mathrm{scale}\_L^2)\$\$ Factors are normalised to mean 1 within the
simulated cell batch.
