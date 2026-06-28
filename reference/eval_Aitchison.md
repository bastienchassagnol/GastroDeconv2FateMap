# Aitchison distance between two compositional vectors

For a composition \\p\\, the centred log-ratio transform is: \$\$
\mathrm{clr}(p)\_j = \log(p\_{j}) - \frac{1}{J}\sum\_{k =
1}^{J}\log(p\_{k}). \$\$

## Usage

``` r
eval_Aitchison(p_obs, p_estimated, min_ratio = 1e-09, trim_shared_zeros = TRUE)
```

## Arguments

- p_obs:

  Numeric vector of observed cellular ratios.

- p_estimated:

  Numeric vector of estimated cellular ratios.

- min_ratio:

  Numeric pseudo-count replacing ratios below this value before the
  centred log-ratio transform. Aitchison distance is undefined for zero
  components.

- trim_shared_zeros:

  Logical; if `TRUE`, entries where both vectors are exactly zero are
  removed before normalisation (e.g. jointly missing cell types).

## Value

A numeric scalar score.

## Details

The Aitchison distance is the Euclidean distance between the two centred
log-ratio vectors: \$\$ d_A(p^{obs}, \hat{p}) = \sqrt{ \sum\_{j = 1}^{J}
\left\[ \mathrm{clr}(p^{obs})\_j - \mathrm{clr}(\hat{p})\_j \right\]^2
}. \$\$
