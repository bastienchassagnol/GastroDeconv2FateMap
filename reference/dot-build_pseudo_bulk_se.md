# Build a SummarizedExperiment from pseudo-bulk counts

Assembles raw count matrix, sample metadata, and gene annotations into a
SummarizedExperiment object.

## Usage

``` r
.build_pseudo_bulk_se(count_matrix, col_data, feature_data)
```

## Arguments

- count_matrix:

  Numeric matrix (genes x samples).

- col_data:

  Sample metadata `data.frame`.

- feature_data:

  Gene metadata `DataFrame`.

## Value

A SummarizedExperiment with assay `"counts"`.
