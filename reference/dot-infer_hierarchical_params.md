# Hierarchical parameter estimation (phenotype + barcode)

Estimates generative-model parameters at the biological level
(phenotype) and, for the log-normal model only, at the technical level
(barcode library-size scaling).

## Usage

``` r
.infer_hierarchical_params(
  counts,
  meta,
  phenotype_col,
  barcode_col,
  model = c("lognormal", "negative_binomial")
)
```

## Arguments

- counts:

  Count matrix (genes x cells).

- meta:

  Cell metadata.

- phenotype_col:

  Phenotype column name.

- barcode_col:

  Barcode column name.

- model:

  `"lognormal"` or `"negative_binomial"`.

## Value

List with `biological` and `technical` components.

## Details

Infer biological and technical simulation parameters

**Biological level** (per phenotype \\p\\):

- Log-normal: \\(\hat{\mu}\_g, \hat{\sigma}\_g)\\ on \\\log(1+x)\\
  scale.

- Negative binomial: per-gene \\(\hat{\mu}\_g, \hat{\theta}\_g)\\ plus
  [`splatter::splatEstimate`](http://oshlacklab.com/splatter/reference/splatEstimate.md)
  for library-size parameters.

**Technical level** (log-normal only): per barcode \\b\\, \$\$s_b =
\frac{\sum\_{c: b(c)=b} L_c}{\overline{\sum\_{c} L_c}}\$\$ where \\L_c =
\sum_g x\_{gc}.\\ This factor is applied when simulating log-normal
cells anchored to barcode \\b\\.

For the negative-binomial path, barcode library scaling is *not* applied
separately: Splatter already estimates library-size variation via
`splatEstimate`, and an extra scaling factor would be redundant.
