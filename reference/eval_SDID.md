# Square-root sine distance between two compositional vectors

This metric first computes the cosine similarity between `p_obs` and
`p_estimated`: \$\$ \cos(\theta) = \frac{ \sum\_{j = 1}^{J}
p^{obs}\_{j}\hat{p}\_{j} }{ \lVert p^{obs} \rVert_2 \lVert \hat{p}
\rVert_2 }. \$\$

## Usage

``` r
eval_SDID(p_obs, p_estimated, trim_shared_zeros = TRUE)
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

The returned score is: \$\$ \mathrm{SDID} = \sqrt{ \sqrt{ 1 - \min(1,
\cos^2(\theta)) } }. \$\$

This keeps the positive distance scale rather than changing the sign as
in the MPRA formulation.

## References

<https://mpra.ub.uni-muenchen.de/84387/>
