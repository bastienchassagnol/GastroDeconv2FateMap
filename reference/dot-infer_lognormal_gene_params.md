# Log-normal parameter estimation on `log1p` counts

Estimates per-gene mean and standard deviation on the \\\log(1+x)\\
scale within a cell subset. Used for the generative log-normal
pseudo-bulk path (strategy 7).

## Usage

``` r
.infer_lognormal_gene_params(counts, cell_idx, max_cells = 2000L)
```

## Arguments

- counts:

  Sparse or dense count matrix (genes x cells).

- cell_idx:

  Integer indices of cells in the estimation subset.

- max_cells:

  Maximum cells used for estimation (random subsample).

## Value

List with `mean_log` and `sd_log` numeric vectors.

## Details

Infer log-normal gene parameters from single cells

For gene \\g\\ and cells \\c \in \mathcal{C}\\, with subsampled counts
\\x\_{gc}\\: \$\$\hat{\mu}\_g = \frac{1}{\|\mathcal{C}\|}\sum\_{c \in
\mathcal{C}} \log(1 + x\_{gc})\$\$ \$\$\hat{\sigma}\_g = \mathrm{sd}\_{c
\in \mathcal{C}} \left(\log(1 + x\_{gc})\right)\$\$

Simulated counts are drawn as \\x\_{gc}^{\mathrm{sim}} = \max\\0,
\lfloor e^{\tilde{z}\_{gc}} - 1 \rfloor\\\\ with \\\tilde{z}\_{gc} \sim
\mathcal{N}(\hat{\mu}\_g, \hat{\sigma}\_g^2)\\.

This stabilises zeros before log transformation, following common
pseudo-bulk pipelines (strategy 7 in
`vignettes/pseudo_bulk_generation.qmd`) and mean-pseudo-bulk approaches
using `log2(counts + epsilon)`.

## References

<https://oshlacklab.com/splatter/articles/splatter.html> (log-normal
library-size modelling motivation).
