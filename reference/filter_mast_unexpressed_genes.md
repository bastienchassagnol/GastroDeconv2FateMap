# Filter single-cell genes for MAST by cell type

Implements a sample-aware detection filter before MAST: a gene is
retained for a cell type when it is detected in at least
`min_cell_detection` of cells for enough biological replicates in at
least one `Morphotype`. This avoids strong global filters that mostly
keep housekeeping genes.

## Usage

``` r
filter_mast_unexpressed_genes(
  seurat_obj,
  celltype_col = "luque_cluster_annotation",
  phenotype_col = "Morphotype",
  barcode_col = "Sample.barcode",
  assay = "RNA",
  layer = "data",
  cell_mask = NULL,
  min_cell_detection = 0.05,
  min_sample_fraction = 0.5,
  min_total_positive_cells = 20L,
  min_cells_per_group = 10L
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

- layer:

  Assay layer used for detection. Default `"data"`; use `"counts"` for
  raw-count detection.

- cell_mask:

  Optional logical vector aligned with cells in `seurat_obj`; only
  selected cells are used.

- min_cell_detection:

  Minimum within-sample cell detection fraction.

- min_sample_fraction:

  Fraction of biological replicates per phenotype that must pass
  `min_cell_detection`.

- min_total_positive_cells:

  Minimum number of positive cells in the tested cell type.

- min_cells_per_group:

  Minimum cells required in a barcode x phenotype x cell-type group
  before that group contributes to the replicate count.

## Value

Named list. Names are cell types and values are selected gene names.
