# Extract gene feature metadata from a Seurat assay

Returns `assay@meta.features`, or a minimal
[`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
of gene names when feature metadata is empty.

## Usage

``` r
.get_feature_metadata(seurat_obj, assay = "RNA")
```

## Arguments

- seurat_obj:

  A Seurat object.

- assay:

  Assay name.

## Value

A `DataFrame` with one row per gene.
