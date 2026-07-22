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
#' @importFrom glue glue
summarise_seurat_assays_layers <- function(seurat_obj) {
  # ======================================================================
  # Shared setup ----
  # ======================================================================
  assay_nms <- names(methods::slot(seurat_obj, "assays"))
  def_assay <- SeuratObject::DefaultAssay(seurat_obj)
  cat(
    glue::glue(
      "Seurat object: {length(assay_nms)} assay(s); ",
      "DefaultAssay = \"{def_assay}\"\n"
    )
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


#' Read a doubly gzip-compressed RDS file
#'
#' @description
#' Loads an R object from an \code{.rds.gz} file that has been compressed
#' twice with gzip. A single \code{\link[base]{readRDS}} call on
#' \code{\link[base]{gzfile}} only removes the outer compression layer; the
#' remaining bytes are still gzip-encoded, which triggers an *unknown input
#' format* error. This helper decompresses both layers in memory, then calls
#' \code{\link[base]{readRDS}}.
#'
#' Some GEO supplementary archives (e.g. per-sample files inside
#' \code{GSE250136_RAW.tar}) are stored in this double-gzip form.
#'
#' @param path Character scalar. Path to the \code{.rds.gz} file.
#'
#' @return The R object stored in the RDS file (class depends on the file;
#'   often a \pkg{Seurat} object for single-cell supplementary data).
#'
#' @details
#' Peak memory use is roughly the size of the first decompressed gzip layer
#' plus the second. For large objects, ensure sufficient RAM is available.
#'
#' @seealso \link[base]{readRDS}, \link[base]{gzfile}, \link[base]{gzcon},
#'   \link[base]{rawConnection}, \link[base]{readBin}
#'
#' @examples
#' \dontrun{
#' obj <- read_double_gz_rds(
#'   "./data/raw/GSE250136/GSM7974412_df48_final.rds.gz"
#' )
#' }
#'
#' @export
read_double_gz_rds <- function(path) {
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  inner_gz <- readBin(con, "raw", n = 1e9)
  readRDS(gzcon(rawConnection(inner_gz)))
}


#' Realign Luque merged Seurat colnames to RNA assay barcodes
#'
#' Merged Luque objects store numeric colnames while the RNA assay uses cell
#' barcodes. Recreates the object from RNA counts and realigned metadata, and
#' re-attaches the SCT assay when present.
#'
#' @param obj A \pkg{Seurat} object from the Luque merged time-course.
#'
#' @return A \pkg{Seurat} object with colnames matching the RNA assay.
#'
#' @keywords internal
realign_luque_seurat <- function(obj) {
  counts <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "counts")
  meta <- obj@meta.data
  rownames(meta) <- colnames(counts)

  sct_assay <- NULL
  if ("SCT" %in% SeuratObject::Assays(obj)) {
    sct_assay <- obj[["SCT"]]
  }

  out <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)

  if (!is.null(sct_assay)) {
    sct_cells <- intersect(colnames(out), colnames(sct_assay))
    out[["SCT"]] <- SeuratObject::CreateAssayObject(
      counts = SeuratObject::GetAssayData(sct_assay, layer = "counts")[, sct_cells]
    )
    sct_data <- SeuratObject::GetAssayData(sct_assay, layer = "data")
    if (ncol(sct_data) > 0L) {
      out[["SCT"]] <- SeuratObject::SetAssayData(
        out[["SCT"]],
        layer = "data",
        new.data = sct_data[, sct_cells]
      )
    }
  }

  out
}
