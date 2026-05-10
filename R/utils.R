#' Row-wise normalisation to the unit simplex
#'
#' @description
#' Rescales each **row** of a numeric matrix-like object so that non-degenerate
#' rows sum to 1 (`\\sum_j x_{ij} = 1`). This enforces mixture proportions along
#' the simplex commonly used after bulk deconvolution when each row is one bulk
#' sample and columns are cell types (or analogous factors).
#'
#' @details
#' **Logic**
#' \enumerate{
#'   \item Coerce to a numeric matrix (`as.matrix`).
#'   \item For each row, compute the row sum (`rowSums`, with `na.rm = TRUE` so a
#'     row dominated by missing values still yields a usable sum).
#' }
#'
#' If every row fails the positivity check, **`x`** is returned unchanged via the
#' early `return`; otherwise the altered rows are merged back into **`m`** and
#' **`m`** is returned (possibly with mixed normalised / untouched rows).
#'
#' @param x A matrix-like object containing non-negative cell estimates (e.g. from omnideconv).
#'
#' @return A numeric matrix of the same size as `as.matrix(x)` (after coercion).
#'

#'
#' @keywords internal
normalise_cell_estimates <- function(x) {
  # Coerce first so downstream logic always sees consistent dimnames/type.
  m <- as.matrix(x)
  rs <- rowSums(m, na.rm = TRUE)
  idx <- rs > 0 & is.finite(rs)
  if (!any(idx)) {
    return(x)
  }
  # Row-major rescaling so each usable row divides by sum_j x_ij (unit simplex).
  m[idx, ] <- sweep(m[idx, , drop = FALSE], 1L, rs[idx], FUN = "/")
  m
}
