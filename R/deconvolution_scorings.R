# All metrics are based on observed and estimated cell type proportions.
# They are adapted from the HADACA3 framework:
# https://github.com/bioinfo-LIG/hadaca3_framework/blob/main/wrapper/06_scoring.R

#' Global Pearson correlation for cellular proportion matrices
#'
#' `p_obs` and `p_estimated` must have the same dimensions. Each column
#' \eqn{i \in \{1, \ldots, N\}} is one cellular profile across
#' \eqn{J} cell types, and each row \eqn{j \in \{1, \ldots, J\}} is a cell
#' type. In other words, the input matrices are \eqn{J \times N}, with
#' \eqn{p^{obs}_{ji}} and \eqn{\hat{p}_{ji}} denoting the observed and
#' estimated ratios for cell type \eqn{j} in sample \eqn{i}.
#'
#' The global Pearson correlation is computed after vectorising both
#' matrices:
#' \deqn{
#' r =
#' \mathrm{cor}\left(
#'   \mathrm{vec}(p^{obs}),
#'   \mathrm{vec}(\hat{p})
#' \right).
#' }
#'
#' If all estimated proportions have zero variance, the worst score `-1` is
#' returned.
#'
#' @param p_obs Numeric matrix of observed cellular ratios. Rows are cell
#'   types and columns are samples.
#' @param p_estimated Numeric matrix of estimated cellular ratios with the
#'   same dimensions and ordering as `p_obs`.
#'
#' @return A numeric scalar score.
#' @export
correlationP_tot <- function(p_obs, p_estimated) {
  if (all(stats::var(c(p_estimated)) == 0)) {
    return(-1)
  }

  stats::cor(c(p_obs), c(p_estimated), method = "pearson")
}

#' Global Spearman correlation for cellular proportion matrices
#'
#' @inheritParams correlationP_tot
#'
#' The global Spearman correlation is computed after vectorising both
#' matrices and ranking their entries:
#' \deqn{
#' \rho =
#' \mathrm{cor}\left(
#'   \mathrm{rank}(\mathrm{vec}(p^{obs})),
#'   \mathrm{rank}(\mathrm{vec}(\hat{p}))
#' \right).
#' }
#'
#' If all estimated proportions have zero variance, the worst score `-1` is
#' returned.
#'
#' @return A numeric scalar score.
#' @export
correlationS_tot <- function(p_obs, p_estimated) {
  if (all(stats::var(c(p_estimated)) == 0)) {
    return(-1)
  }

  stats::cor(c(p_obs), c(p_estimated), method = "spearman")
}

#' Mean square-root sine distance between observed and estimated profiles
#'
#' @inheritParams correlationP_tot
#'
#' For each sample \eqn{i}, this metric first computes the cosine similarity
#' between the observed and estimated cellular profiles:
#' \deqn{
#' \cos(\theta_i) =
#' \frac{
#'   \sum_{j = 1}^{J} p^{obs}_{ji}\hat{p}_{ji}
#' }{
#'   \lVert p^{obs}_{\cdot i} \rVert_2
#'   \lVert \hat{p}_{\cdot i} \rVert_2
#' }.
#' }
#'
#' It then reports the mean sample-wise square-root sine distance:
#' \deqn{
#' \mathrm{SDID} =
#' \frac{1}{N}
#' \sum_{i = 1}^{N}
#' \sqrt{
#'   \sqrt{
#'     1 - \min(1, \cos^2(\theta_i))
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
eval_SDID <- function(p_obs, p_estimated) {
  cos_angle <- vapply(
    seq_len(ncol(p_obs)),
    \(i) {
      as.numeric(
        p_obs[, i] %*% p_estimated[, i] /
          base::norm(p_obs[, i, drop = FALSE], type = "2") /
          base::norm(p_estimated[, i, drop = FALSE], type = "2")
      )
    },
    numeric(1)
  )

  sin_angle <- vapply(
    cos_angle,
    \(x) sqrt(1 - base::min(1, x^2)),
    numeric(1)
  )

  mean(sqrt(sin_angle), na.rm = TRUE)
}

#' Mean Aitchison distance between observed and estimated profiles
#'
#' @inheritParams correlationP_tot
#' @param min_ratio Numeric pseudo-count replacing ratios below this value
#'   before the centred log-ratio transform. Aitchison distance is undefined for
#'   zero components.
#'
#' For a composition \eqn{p_{\cdot i}}, the centred log-ratio transform is:
#' \deqn{
#' \mathrm{clr}(p_{\cdot i})_j =
#' \log(p_{ji}) -
#' \frac{1}{J}\sum_{k = 1}^{J}\log(p_{ki}).
#' }
#'
#' The per-sample Aitchison distance is the Euclidean distance between the
#' centred log-ratio profiles:
#' \deqn{
#' d_A(p^{obs}_{\cdot i}, \hat{p}_{\cdot i}) =
#' \sqrt{
#'   \sum_{j = 1}^{J}
#'   \left[
#'     \mathrm{clr}(p^{obs}_{\cdot i})_j -
#'     \mathrm{clr}(\hat{p}_{\cdot i})_j
#'   \right]^2
#' }.
#' }
#'
#' The returned score is the mean across samples:
#' \deqn{
#' \frac{1}{N}\sum_{i = 1}^{N}
#' d_A(p^{obs}_{\cdot i}, \hat{p}_{\cdot i}).
#' }
#'
#' @return A numeric scalar score.
#' @export
eval_Aitchison <- function(p_obs, p_estimated, min_ratio = 1e-9) {
  p_obs[p_obs < min_ratio] <- min_ratio
  p_estimated[p_estimated < min_ratio] <- min_ratio

  distances <- vapply(
    seq_len(ncol(p_obs)),
    \(i) {
      clr_obs <- .calculate_clr(p_obs[, i])
      clr_estimated <- .calculate_clr(p_estimated[, i])

      sqrt(sum((clr_obs - clr_estimated)^2))
    },
    numeric(1)
  )

  mean(distances, na.rm = TRUE)
}

#' Mean Jensen-Shannon divergence between observed and estimated profiles
#'
#' @inheritParams correlationP_tot
#' @inheritParams eval_Aitchison
#'
#' For each sample \eqn{i}, define the mixture composition:
#' \deqn{
#' m_{ji} = \frac{1}{2}(p^{obs}_{ji} + \hat{p}_{ji}).
#' }
#'
#' The Jensen-Shannon divergence is:
#' \deqn{
#' \mathrm{JSD}_i =
#' \frac{1}{2}
#' \sum_{j = 1}^{J}
#' p^{obs}_{ji}
#' \log\left(\frac{p^{obs}_{ji}}{m_{ji}}\right)
#' +
#' \frac{1}{2}
#' \sum_{j = 1}^{J}
#' \hat{p}_{ji}
#' \log\left(\frac{\hat{p}_{ji}}{m_{ji}}\right).
#' }
#'
#' The returned score is \eqn{N^{-1}\sum_{i = 1}^{N}\mathrm{JSD}_i}.
#'
#' @return A numeric scalar score.
#' @export
eval_JSD <- function(p_obs, p_estimated, min_ratio = 1e-9) {
  jsd <- vapply(
    seq_len(ncol(p_obs)),
    \(i) {
      obs <- p_obs[, i]
      estimated <- p_estimated[, i]

      obs[obs < min_ratio] <- min_ratio
      estimated[estimated < min_ratio] <- min_ratio

      mixture <- 0.5 * (obs + estimated)
      0.5 * (
        sum(obs * log(obs / mixture)) +
          sum(estimated * log(estimated / mixture))
      )
    },
    numeric(1)
  )

  mean(jsd, na.rm = TRUE)
}

#' Root mean squared error for cellular proportion matrices
#'
#' @inheritParams correlationP_tot
#'
#' The root mean squared error is:
#' \deqn{
#' \mathrm{RMSE} =
#' \sqrt{
#'   \frac{1}{JN}
#'   \sum_{j = 1}^{J}\sum_{i = 1}^{N}
#'   \left(p^{obs}_{ji} - \hat{p}_{ji}\right)^2
#' }.
#' }
#'
#' @return A numeric scalar score.
#' @export
eval_RMSE <- function(p_obs, p_estimated) {
  sqrt(mean((p_obs - p_estimated)^2, na.rm = TRUE))
}

#' Centre-scale a numeric score vector to `[0, 1]`
#'
#' @param x Numeric vector of scores.
#'
#' @return A numeric vector with the same length as `x`.
#' @export
CenterScaleNorm <- function(x) {
  centred <- x - base::min(x, na.rm = TRUE)
  centred / base::max(centred, na.rm = TRUE)
}

.calculate_clr <- function(x) {
  log_x <- log(x)
  log_x - mean(log_x)
}
