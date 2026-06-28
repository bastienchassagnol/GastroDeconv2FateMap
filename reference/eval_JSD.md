# Jensen-Shannon divergence between two compositional vectors

Define the mixture composition: \$\$ m\_{j} = \frac{1}{2}(p^{obs}\_{j} +
\hat{p}\_{j}). \$\$

## Usage

``` r
eval_JSD(p_obs, p_estimated, min_ratio = 1e-09, trim_shared_zeros = TRUE)
```

## Arguments

- p_obs:

  Numeric vector of observed cellular ratios.

- p_estimated:

  Numeric vector of estimated cellular ratios.

- min_ratio:

  Numeric pseudo-count replacing ratios below this value before
  computing logarithms. Jensen-Shannon divergence is undefined for zero
  components.

- trim_shared_zeros:

  Logical; if `TRUE`, entries where both vectors are exactly zero are
  removed before normalisation (e.g. jointly missing cell types).

## Value

A numeric scalar score.

## Details

The Jensen-Shannon divergence is: \$\$ \mathrm{JSD} = \frac{1}{2}
\sum\_{j = 1}^{J} p^{obs}\_{j}
\log\left(\frac{p^{obs}\_{j}}{m\_{j}}\right) + \frac{1}{2} \sum\_{j =
1}^{J} \hat{p}\_{j} \log\left(\frac{\hat{p}\_{j}}{m\_{j}}\right). \$\$
