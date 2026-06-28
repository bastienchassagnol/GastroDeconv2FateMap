# Validate two compositional vectors

Validate two compositional vectors

## Usage

``` r
.validate_compositions(p_obs, p_estimated, trim_shared_zeros = TRUE)
```

## Arguments

- p_obs:

  Numeric vector of observed cellular ratios.

- p_estimated:

  Numeric vector of estimated cellular ratios.

- trim_shared_zeros:

  Logical; if `TRUE`, entries where both vectors are exactly zero are
  removed before normalisation (e.g. jointly missing cell types).

## Value

A named list with normalised vectors `p_obs` and `p_estimated`.
