# All metrics are based on observed and estimated cell type proportions.
# They are adapted from the HADACA3 framework:
# https://github.com/bioinfo-LIG/hadaca3_framework/blob/main/wrapper/06_scoring.R

#' Pearson correlation between two compositional vectors
#'
#' `p_obs` and `p_estimated` must be numeric vectors of the same length,
#' representing observed and estimated cellular ratios for the same set of cell
#' types.
#'
#' The Pearson correlation is:
#' \deqn{
#' r = \mathrm{cor}(p^{obs}, \hat{p}).
#' }
#'
#' If all estimated proportions have zero variance, the worst score `-1` is
#' returned.
#'
#' @inheritParams .validate_compositions
#'
#' @return A numeric scalar score.
#' @export
eval_Pearson <- function(p_obs, p_estimated, trim_shared_zeros = TRUE) {
  validated <- .validate_compositions(
    p_obs = p_obs,
    p_estimated = p_estimated,
    trim_shared_zeros = trim_shared_zeros
  )
  p_obs <- validated$p_obs
  p_estimated <- validated$p_estimated

  if (stats::var(p_estimated) == 0) {
    return(-1)
  }

  stats::cor(p_obs, p_estimated, method = "pearson")
}

#' Square-root sine distance between two compositional vectors
#'
#' @inheritParams .validate_compositions
#'
#' This metric first computes the cosine similarity between `p_obs` and
#' `p_estimated`:
#' \deqn{
#' \cos(\theta) =
#' \frac{
#'   \sum_{j = 1}^{J} p^{obs}_{j}\hat{p}_{j}
#' }{
#'   \lVert p^{obs} \rVert_2
#'   \lVert \hat{p} \rVert_2
#' }.
#' }
#'
#' The returned score is:
#' \deqn{
#' \mathrm{SDID} =
#' \sqrt{
#'   \sqrt{
#'     1 - \min(1, \cos^2(\theta))
#'   }
#' }.
#' }
#'
#' This keeps the positive distance scale rather than changing the sign as in
#' the MPRA formulation.
#'
#' @return A numeric scalar score.
#' @references \url{https://mpra.ub.uni-muenchen.de/84387/}
#' @export
eval_SDID <- function(p_obs, p_estimated, trim_shared_zeros = TRUE) {
  validated <- .validate_compositions(
    p_obs = p_obs,
    p_estimated = p_estimated,
    trim_shared_zeros = trim_shared_zeros
  )
  p_obs <- validated$p_obs
  p_estimated <- validated$p_estimated

  denom <- base::norm(matrix(p_obs, ncol = 1), type = "2") *
    base::norm(matrix(p_estimated, ncol = 1), type = "2")

  if (denom == 0) {
    return(NA_real_)
  }

  cos_angle <- as.numeric(sum(p_obs * p_estimated) / denom)
  sin_angle <- sqrt(1 - base::min(1, cos_angle^2))
  sqrt(sin_angle)
}

#' Aitchison distance between two compositional vectors
#'
#' @inheritParams .validate_compositions
#' @param min_ratio Numeric pseudo-count replacing ratios below this value
#'   before the centred log-ratio transform. Aitchison distance is undefined for
#'   zero components.
#'
#' For a composition \eqn{p}, the centred log-ratio transform is:
#' \deqn{
#' \mathrm{clr}(p)_j =
#' \log(p_{j}) -
#' \frac{1}{J}\sum_{k = 1}^{J}\log(p_{k}).
#' }
#'
#' The Aitchison distance is the Euclidean distance between the two centred
#' log-ratio vectors:
#' \deqn{
#' d_A(p^{obs}, \hat{p}) =
#' \sqrt{
#'   \sum_{j = 1}^{J}
#'   \left[
#'     \mathrm{clr}(p^{obs})_j -
#'     \mathrm{clr}(\hat{p})_j
#'   \right]^2
#' }.
#' }
#'
#' @return A numeric scalar score.
#' @export
eval_Aitchison <- function(
    p_obs,
    p_estimated,
    min_ratio = 1e-9,
    trim_shared_zeros = TRUE
) {
  validated <- .validate_compositions(
    p_obs = p_obs,
    p_estimated = p_estimated,
    trim_shared_zeros = trim_shared_zeros
  )
  p_obs <- validated$p_obs
  p_estimated <- validated$p_estimated

  p_obs[p_obs < min_ratio] <- min_ratio
  p_estimated[p_estimated < min_ratio] <- min_ratio

  clr_obs <- .calculate_clr(p_obs)
  clr_estimated <- .calculate_clr(p_estimated)

  sqrt(sum((clr_obs - clr_estimated)^2))
}

#' Jensen-Shannon divergence between two compositional vectors
#'
#' @inheritParams .validate_compositions
#' @inheritParams eval_Aitchison
#'
#' Define the mixture composition:
#' \deqn{
#' m_{j} = \frac{1}{2}(p^{obs}_{j} + \hat{p}_{j}).
#' }
#'
#' The Jensen-Shannon divergence is:
#' \deqn{
#' \mathrm{JSD} =
#' \frac{1}{2}
#' \sum_{j = 1}^{J}
#' p^{obs}_{j}
#' \log\left(\frac{p^{obs}_{j}}{m_{j}}\right)
#' +
#' \frac{1}{2}
#' \sum_{j = 1}^{J}
#' \hat{p}_{j}
#' \log\left(\frac{\hat{p}_{j}}{m_{j}}\right).
#' }
#'
#' @return A numeric scalar score.
#' @export
eval_JSD <- function(
    p_obs,
    p_estimated,
    min_ratio = 1e-9,
    trim_shared_zeros = TRUE
) {
  validated <- .validate_compositions(
    p_obs = p_obs,
    p_estimated = p_estimated,
    trim_shared_zeros = trim_shared_zeros
  )
  p_obs <- validated$p_obs
  p_estimated <- validated$p_estimated

  obs <- p_obs
  estimated <- p_estimated
  obs[obs < min_ratio] <- min_ratio
  estimated[estimated < min_ratio] <- min_ratio

  mixture <- 0.5 * (obs + estimated)
  0.5 * (sum(obs * log(obs / mixture)) + sum(estimated * log(estimated / mixture)))
}

#' Root mean squared error for cellular proportion vectors
#'
#' @inheritParams .validate_compositions
#'
#' The root mean squared error is:
#' \deqn{
#' \mathrm{RMSE} =
#' \sqrt{
#'   \frac{1}{J}
#'   \sum_{j = 1}^{J}
#'   \left(p^{obs}_{j} - \hat{p}_{j}\right)^2
#' }.
#' }
#'
#' @return A numeric scalar score.
#' @export
eval_RMSE <- function(p_obs, p_estimated, trim_shared_zeros = TRUE) {
  validated <- .validate_compositions(
    p_obs = p_obs,
    p_estimated = p_estimated,
    trim_shared_zeros = trim_shared_zeros
  )
  p_obs <- validated$p_obs
  p_estimated <- validated$p_estimated
  sqrt(mean((p_obs - p_estimated)^2, na.rm = TRUE))
}

#' Mean absolute error for cellular proportion vectors
#'
#' @inheritParams .validate_compositions
#'
#' The mean absolute error is:
#' \deqn{
#' \mathrm{MAE} =
#' \frac{1}{J}
#' \sum_{j = 1}^{J}
#' \left|p^{obs}_{j} - \hat{p}_{j}\right|.
#' }
#'
#' @return A numeric scalar score.
#' @export
eval_MAE <- function(p_obs, p_estimated, trim_shared_zeros = TRUE) {
  validated <- .validate_compositions(
    p_obs = p_obs,
    p_estimated = p_estimated,
    trim_shared_zeros = trim_shared_zeros
  )
  p_obs <- validated$p_obs
  p_estimated <- validated$p_estimated
  mean(abs(p_obs - p_estimated), na.rm = TRUE)
}

#' Validate two compositional vectors
#'
#' @param p_obs Numeric vector of observed cellular ratios.
#' @param p_estimated Numeric vector of estimated cellular ratios.
#' @param trim_shared_zeros Logical; if `TRUE`, entries where both vectors are
#'   exactly zero are removed before normalisation (e.g. jointly missing cell
#'   types).
#'
#' @return A named list with normalised vectors `p_obs` and `p_estimated`.
#' @keywords internal
.validate_compositions <- function(
    p_obs,
    p_estimated,
    trim_shared_zeros = TRUE
) {
  if (!is.logical(trim_shared_zeros) || length(trim_shared_zeros) != 1L ||
    is.na(trim_shared_zeros)) {
    stop("`trim_shared_zeros` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(p_obs) || !is.numeric(p_estimated)) {
    stop("`p_obs` and `p_estimated` must be numeric vectors.", call. = FALSE)
  }
  if (!is.vector(p_obs) || !is.vector(p_estimated)) {
    stop("`p_obs` and `p_estimated` must be vectors, not matrices.", call. = FALSE)
  }
  if (length(p_obs) != length(p_estimated)) {
    stop("`p_obs` and `p_estimated` must have the same length.", call. = FALSE)
  }
  if (length(p_obs) == 0) {
    stop("`p_obs` and `p_estimated` must be non-empty vectors.", call. = FALSE)
  }
  if (is.null(names(p_obs)) || is.null(names(p_estimated))) {
    stop("`p_obs` and `p_estimated` must be named vectors.", call. = FALSE)
  }
  if (anyNA(names(p_obs)) || anyNA(names(p_estimated)) ||
    any(names(p_obs) == "") || any(names(p_estimated) == "")) {
    stop("`p_obs` and `p_estimated` names must be non-missing and non-empty.", call. = FALSE)
  }
  if (!identical(names(p_obs), names(p_estimated))) {
    stop("`p_obs` and `p_estimated` must have exactly matching names in the same order.", call. = FALSE)
  }
  if (trim_shared_zeros) {
    keep_idx <- !(p_obs == 0 & p_estimated == 0)
    p_obs <- p_obs[keep_idx]
    p_estimated <- p_estimated[keep_idx]
  }
  if (length(p_obs) == 0) {
    stop("No entries remain after trimming shared zeros.", call. = FALSE)
  }

  list(
    p_obs = normalise_cell_estimates(p_obs),
    p_estimated = normalise_cell_estimates(p_estimated)
  )
}

#' Normalise a composition vector to the unit simplex
#'
#' @description
#' Rescales a non-negative numeric vector so that its entries sum to 1
#' (`\\sum_j x_j = 1`), i.e. the unit simplex constraint for compositions.
#'
#' @details
#' If the input already satisfies the simplex constraint (within
#' `sqrt(.Machine$double.eps)`), it is returned unchanged. Otherwise, it is
#' divided by its sum.
#'
#' @param x Numeric vector containing non-negative cell estimates.
#'
#' @return A numeric vector on the unit simplex with the same names as `x`.
#'
#' @keywords internal
normalise_cell_estimates <- function(x) {

  if (length(x) == 0) {
    stop("`x` must be a non-empty vector.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`x` must contain only finite, non-missing values.", call. = FALSE)
  }

  total <- sum(x)
  if (!is.finite(total) || total <= 0) {
    stop("`x` must have a strictly positive sum.", call. = FALSE)
  }

  tolerance <- sqrt(.Machine$double.eps)
  if (abs(total - 1) <= tolerance) {
    return(x)
  }

  x / total
}

# Helper function to calculate the centred log-ratio transform
.calculate_clr <- function(x) {
  log_x <- log(x)
  log_x - mean(log_x)
}
