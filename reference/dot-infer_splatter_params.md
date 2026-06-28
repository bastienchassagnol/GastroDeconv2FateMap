# Splatter `splatEstimate` on a phenotype cell subset

Wraps
[`splatter::splatEstimate`](http://oshlacklab.com/splatter/reference/splatEstimate.md)
to infer global scRNA-seq simulation parameters, including library-size
location and scale.

## Usage

``` r
.infer_splatter_params(counts, cell_idx, max_cells = 2000L)
```

## Arguments

- counts:

  Sparse or dense count matrix (genes x cells).

- cell_idx:

  Integer indices of cells in the estimation subset.

- max_cells:

  Maximum cells passed to `splatEstimate`.

## Value

A
[`splatter::SplatParams`](http://oshlacklab.com/splatter/reference/SplatParams.md)
object.

## Details

Estimate Splatter simulation parameters from single cells

Splatter models total library size \\L_c\\ per cell (often log-normal)
and gene-wise means as functions of \\L_c\\. `splatEstimate` returns a
`SplatParams` object with, among others, library-size
\\(\mathrm{loc}\_L, \mathrm{scale}\_L)\\ used to draw technical depth
variation at simulation time.

## References

<https://oshlacklab.com/splatter/reference/splatEstimate.html>; Zappia
et al. (2017)
[doi:10.1186/s13059-017-1303-1](https://doi.org/10.1186/s13059-017-1303-1)
.
