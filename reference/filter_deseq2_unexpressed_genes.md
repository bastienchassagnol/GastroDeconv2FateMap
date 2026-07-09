# Filter pseudo-bulk genes for DESeq2 by cell type

Aggregates raw counts by cell type and biological replicate, analogous
to
[`Seurat::AggregateExpression()`](https://satijalab.org/seurat/reference/AggregateExpression.html),
then applies
[`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html)
within that cell type. This keeps the tested gene universe aligned with
cell-type-specific pseudo-bulk `Morphotype` contrasts. Data-adaptive
methods such as HTSFilter or IHW can be added downstream, but are not
applied here.

## Usage

``` r
filter_deseq2_unexpressed_genes(
  seurat_obj,
  celltype_col = "luque_cluster_annotation",
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  assay = "RNA",
  cell_mask = NULL,
  min_cells_per_sample = 10L,
  min_count = 5L,
  min_total_count = 15L,
  large_n = 10L,
  min_prop = 0.7
)
```

## Arguments

- seurat_obj:

  A Seurat object.

- celltype_col, phenotype_col, barcode_col:

  Metadata columns containing cell-type labels, condition labels, and
  biological replicate identifiers.

- assay:

  Assay used for raw counts. Default `"RNA"`.

- cell_mask:

  Optional logical vector aligned with cells in `seurat_obj`; only
  selected cells are used.

- min_cells_per_sample:

  Minimum cells required for a cell-type x barcode pseudo-bulk sample to
  be retained.

- min_count, min_total_count, large_n, min_prop:

  Parameters passed to
  [`edgeR::filterByExpr()`](https://rdrr.io/pkg/edgeR/man/filterByExpr.html).
  The default `min_count = 5` is slightly less stringent than bulk
  RNA-seq defaults because cell-type pseudo-bulk libraries are often
  smaller.

## Value

Named list. Names are cell types and values are selected gene names.
