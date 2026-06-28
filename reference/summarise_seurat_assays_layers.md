# Print a hierarchical summary of assays and layers

**Why RNA can be ~27k genes while `dim(seurat)` is ~2944 features**:

- `dim(seurat)` uses `DefaultAssay()` (here often `"integrated"`).

- The integrated / SCT assays keep the integration anchor feature set.

- Assay `"RNA"` still holds the full transcriptome (all detected genes).

**LayerData** (see
[`?SeuratObject::LayerData`](https://satijalab.github.io/seurat-object/reference/Layers.html)):
main arguments:

- `layer`: e.g. `"counts"`, `"data"`, `"scale.data"`, or `NULL` for the
  default layer(s).

- `features`, `cells`: optional row/column subsets (Assay / Assay5).

- `assay`: when `object` is a Seurat, which assay to query.

- `search` in
  [`SeuratObject::Layers`](https://satijalab.github.io/seurat-object/reference/Layers.html):
  `NA` = all layer names; `NULL` = default layer(s).

**Default layer**:

- Read:
  [`SeuratObject::DefaultLayer`](https://satijalab.github.io/seurat-object/reference/DefaultLayer.html)
  on an assay.

- Set: `DefaultLayer(assay) <- "counts"` (or another existing layer).

- Reading defaults: `LayerData(seurat, assay = "RNA", layer = NULL)`
  uses the default layer for that assay (often `"data"`).

## Usage

``` r
summarise_seurat_assays_layers(seurat_obj)
```

## Arguments

- seurat_obj:

  A Seurat object (v4/v5), with one or more assays.

## Value

`NULL`, invisibly. The summary is printed to the standard output
connection via [`cat`](https://rdrr.io/r/base/cat.html).

## See also

[LayerData](https://satijalab.github.io/seurat-object/reference/Layers.html),
[Layers](https://satijalab.github.io/seurat-object/reference/Layers.html),
[DefaultLayer](https://satijalab.github.io/seurat-object/reference/DefaultLayer.html)
