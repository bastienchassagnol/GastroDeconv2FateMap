# Pearson correlation between two compositional vectors

`p_obs` and `p_estimated` must be numeric vectors of the same length,
representing observed and estimated cellular ratios for the same set of
cell types.

## Usage

``` r
eval_Pearson(p_obs, p_estimated, trim_shared_zeros = TRUE)
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

A numeric scalar score.

## Details

The Pearson correlation is: \$\$ r = \mathrm{cor}(p^{obs}, \hat{p}).
\$\$

If all estimated proportions have zero variance, the worst score `-1` is
returned.
