# Per-gene NB mean and dispersion estimation

Estimates per-gene mean count and NB dispersion within a cell subset,
following the HADACA3 / Splatter count-modelling rationale.

## Usage

``` r
.infer_nb_gene_params(counts, cell_idx, max_cells = 2000L)
```

## Arguments

- counts:

  Sparse or dense count matrix (genes x cells).

- cell_idx:

  Integer indices of cells in the estimation subset.

- max_cells:

  Maximum cells used for estimation (random subsample).

## Value

List with `mu` and `size` numeric vectors.

## Details

Infer negative-binomial gene parameters from single cells

For gene \\g\\ and cells \\c \in \mathcal{C}\\: \$\$\hat{\mu}\_g =
\frac{1}{\|\mathcal{C}\|}\sum\_{c \in \mathcal{C}} x\_{gc}\$\$
\$\$\hat{\theta}\_g = \frac{\hat{\mu}\_g^2}{
\widehat{\mathrm{Var}}(x\_{g\cdot}) - \hat{\mu}\_g}\$\$

Simulated counts use
[`stats::rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)
parameterised by \\(\mu_g, \mathrm{size}=\theta_g)\\.

## References

Zappia et al. (2017) Splatter: simulation of single-cell RNA sequencing
data. *Genome Biology* 18, 174.
[doi:10.1186/s13059-017-1303-1](https://doi.org/10.1186/s13059-017-1303-1)
