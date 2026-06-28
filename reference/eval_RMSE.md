# Root mean squared error for cellular proportion vectors

The root mean squared error is: \$\$ \mathrm{RMSE} = \sqrt{ \frac{1}{J}
\sum\_{j = 1}^{J} \left(p^{obs}\_{j} - \hat{p}\_{j}\right)^2 }. \$\$

## Usage

``` r
eval_RMSE(p_obs, p_estimated, trim_shared_zeros = TRUE)
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
