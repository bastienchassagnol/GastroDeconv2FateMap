# Extract variance components from an hrRMSE mixed model

Retrieves the three variance components estimated by
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) for the hrRMSE
hierarchical model.

## Usage

``` r
.extract_hrrmse_variances(fit)
```

## Arguments

- fit:

  A fitted `lmerMod` object from
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html).

## Value

A named numeric vector with elements `population`, `sample`, and
`residual`.

## Details

For a fitted model of the form \$\$ y\_{ijm} = \mu + \beta\\ I(m =
\text{estimated}) + u_j + v\_{ij} + \varepsilon\_{ijm}, \$\$ the
random-effect standard deviations returned by
[`lme4::VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) are
squared to obtain \$\$\hat\sigma^2\_{\mathrm{population}} =
\mathrm{Var}(u_j),\$\$ \$\$\hat\sigma^2\_{\mathrm{sample}} =
\mathrm{Var}(v\_{ij}),\$\$ from the `(1 | cell_type)` and
`(1 | cell_type:sample)` terms, respectively. The residual variance is
\$\$\hat\sigma_e^2 = \sigma(\mathrm{fit})^2,\$\$ where
`\sigma(\mathrm{fit})` is given by
[`stats::sigma()`](https://rdrr.io/r/stats/sigma.html).
