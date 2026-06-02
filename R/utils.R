
# ============================================================================
# Seurat v5: assay and layer diagnostics ----
# ============================================================================

#' Print a hierarchical summary of assays and layers
#'
#' @description
#' **Why RNA can be ~27k genes while `dim(seurat)` is ~2944 features**:
#' \itemize{
#'   \item `dim(seurat)` uses `DefaultAssay()` (here often `"integrated"`).
#'   \item The integrated / SCT assays keep the integration anchor feature set.
#'   \item Assay `"RNA"` still holds the full transcriptome (all detected
#'     genes).
#' }
#' **LayerData** (see \code{?SeuratObject::LayerData}): main arguments:
#' \itemize{
#'   \item `layer`: e.g. `"counts"`, `"data"`, `"scale.data"`, or `NULL` for
#'     the default layer(s).
#'   \item `features`, `cells`: optional row/column subsets (Assay / Assay5).
#'   \item `assay`: when `object` is a Seurat, which assay to query.
#'   \item \code{search} in \code{SeuratObject::Layers}: \code{NA} = all layer
#'     names; \code{NULL} = default layer(s).
#' }
#' **Default layer**:
#' \itemize{
#'   \item Read: \code{SeuratObject::DefaultLayer} on an assay.
#'   \item Set: `DefaultLayer(assay) <- "counts"` (or another existing layer).
#'   \item Reading defaults: `LayerData(seurat, assay = "RNA", layer = NULL)`
#'     uses the default layer for that assay (often `"data"`).
#' }
#'
#' @param seurat_obj A \pkg{Seurat} object (v4/v5), with one or more assays.
#'
#' @return \code{NULL}, invisibly. The summary is printed to the standard
#'   output connection via \code{\link[base]{cat}}.
#'
#' @seealso \link[SeuratObject]{LayerData}, \link[SeuratObject]{Layers},
#'   \link[SeuratObject]{DefaultLayer}
#'
#' @keywords internal
summarise_seurat_assays_layers <- function(seurat_obj) {
  # ======================================================================
  # Shared setup ----
  # ======================================================================
  assay_nms <- names(methods::slot(seurat_obj, "assays"))
  def_assay <- SeuratObject::DefaultAssay(seurat_obj)
  cat(
    "Seurat object: ",
    length(assay_nms),
    " assay(s); DefaultAssay = \"",
    def_assay,
    "\"\n",
    sep = ""
  )

  # ======================================================================
  # Per-assay and per-layer tree ----
  # ======================================================================
  for (nm in assay_nms) {
    ast <- seurat_obj[[nm]]
    layers_all <- SeuratObject::Layers(ast, search = NA)
    layers_default <- SeuratObject::Layers(ast, search = NULL)
    def_one <- tryCatch(
      SeuratObject::DefaultLayer(ast),
      error = function(e) NA_character_
    )
    suffix <- if (identical(nm, def_assay)) " [DefaultAssay]" else ""
    cat("+- assay: ", nm, suffix, "\n", sep = "")
    for (ly in layers_all) {
      mat <- SeuratObject::LayerData(ast, layer = ly)
      if (is.null(mat)) {
        cat("|  +- layer: ", ly, "  (missing)\n", sep = "")
        next
      }
      d <- dim(mat)
      is_def <- ly %in%
        layers_default ||
        (!is.na(def_one) && identical(ly, def_one))
      tag <- if (is_def) " [default layer]" else ""
      cat(
        "|  +- layer: ",
        ly,
        "  ",
        d[[1L]],
        " x ",
        d[[2L]],
        tag,
        "\n",
        sep = ""
      )
    }
  }
  invisible(NULL)
}
