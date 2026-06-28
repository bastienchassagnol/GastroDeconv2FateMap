# Hierarchical relative RMSE against a shared compositional reference

Computes a conditional, common-reference hierarchical relative RMSE
(hrRMSE) from multiple deconvolution replicates compared with one pooled
gold-standard composition. This is appropriate for technical replicates
of the same mixture, but is a strong assumption for biological
replicates.

## Usage

``` r
eval_hrRMSE(
  p_obs,
  p_estimated,
  trim_shared_zeros = TRUE,
  method = c("REML", "ML")
)
```

## Arguments

- p_obs:

  Named numeric vector of the shared reference composition (\\p_j\\),
  repeated for every sample.

- p_estimated:

  Sample-level estimated compositions. A matrix or data frame with
  samples in rows and cell types in columns (matching `p_obs` names), or
  a named list of named numeric vectors.

- trim_shared_zeros:

  Logical; if `TRUE`, cell types that are zero in both the reference and
  all estimates are removed before fitting.

- method:

  Character; `"REML"` (default) or `"ML"` for variance-component
  estimation.

## Value

A numeric scalar hrRMSE. Attributes `variance_components` (named vector
of estimated variances) and `method` record the fitted components and
estimation method.

## Details

After arranging the data, for sample \\i\\, cell type \\j\\, and
measurement source \\m\\: \$\$ y\_{ijm} = \begin{cases} p_j, & m =
\text{gold standard}, \\ \hat{p}\_{ij}, & m = \text{estimated}.
\end{cases} \$\$

The hierarchical mixed-effects model is: \$\$ y\_{ijm} = \mu + \beta\\
I(m = \text{estimated}) + u_j + v\_{ij} + \varepsilon\_{ijm}, \$\$ with
\$\$u_j \sim \mathcal{N}(0, \sigma^2\_{\mathrm{population}}),\$\$
\$\$v\_{ij} \sim \mathcal{N}(0, \sigma^2\_{\mathrm{sample}}),\$\$
\$\$\varepsilon\_{ijm} \sim \mathcal{N}(0, \sigma_e^2).\$\$

Biological signal variance is \$\$ V_T =
\sigma^2\_{\mathrm{population}} + \sigma^2\_{\mathrm{sample}}, \$\$ and
\$\$ \mathrm{hrRMSE} = \sqrt{\frac{\sigma_e^2}{V_T}}. \$\$

Variance components are estimated by REML or ML. Values below 1 indicate
residual disagreement smaller than the variance separating cell-type
abundances.
