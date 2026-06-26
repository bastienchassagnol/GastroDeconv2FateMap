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

#' Hierarchical relative RMSE against a shared compositional reference
#'
#' @description
#' Computes a conditional, common-reference hierarchical relative RMSE
#' (hrRMSE) from multiple deconvolution replicates compared with one pooled
#' gold-standard composition. This is appropriate for technical replicates of
#' the same mixture, but is a strong assumption for biological replicates.
#'
#' @param p_obs Named numeric vector of the shared reference composition
#'   (\eqn{p_j}), repeated for every sample.
#' @param p_estimated Sample-level estimated compositions. A matrix or data
#'   frame with samples in rows and cell types in columns (matching `p_obs`
#'   names), or a named list of named numeric vectors.
#' @param trim_shared_zeros Logical; if `TRUE`, cell types that are zero in
#'   both the reference and all estimates are removed before fitting.
#' @param method Character; `"REML"` (default) or `"ML"` for variance-component
#'   estimation.
#'
#' @details
#' After arranging the data, for sample \eqn{i}, cell type \eqn{j}, and
#' measurement source \eqn{m}:
#' \deqn{
#' y_{ijm} =
#' \begin{cases}
#' p_j, & m = \text{gold standard}, \\
#' \hat{p}_{ij}, & m = \text{estimated}.
#' \end{cases}
#' }
#'
#' The hierarchical mixed-effects model is:
#' \deqn{
#' y_{ijm} = \mu + \beta\, I(m = \text{estimated}) + u_j + v_{ij} + \varepsilon_{ijm},
#' }
#' with
#' \deqn{u_j \sim \mathcal{N}(0, \sigma^2_{\mathrm{population}}),}
#' \deqn{v_{ij} \sim \mathcal{N}(0, \sigma^2_{\mathrm{sample}}),}
#' \deqn{\varepsilon_{ijm} \sim \mathcal{N}(0, \sigma_e^2).}
#'
#' Biological signal variance is
#' \deqn{
#' V_T = \sigma^2_{\mathrm{population}} + \sigma^2_{\mathrm{sample}},
#' }
#' and
#' \deqn{
#' \mathrm{hrRMSE} = \sqrt{\frac{\sigma_e^2}{V_T}}.
#' }
#'
#' Variance components are estimated by REML or ML. Values below 1 indicate
#' residual disagreement smaller than the variance separating cell-type
#' abundances.
#'
#' @return A numeric scalar hrRMSE. Attributes `variance_components` (named
#'   vector of estimated variances) and `method` record the fitted components
#'   and estimation method.
#' @export
eval_hrRMSE <- function(
    p_obs,
    p_estimated,
    trim_shared_zeros = TRUE,
    method = c("REML", "ML")
) {
  method <- match.arg(method)
  prepared <- .prepare_hrrmse_compositions(
    p_obs = p_obs,
    p_estimated = p_estimated,
    trim_shared_zeros = trim_shared_zeros
  )
  p_obs <- prepared$p_obs
  p_estimated <- prepared$p_estimated

  n_samples <- nrow(p_estimated)
  n_cell_types <- ncol(p_estimated)
  if (n_samples < 2L) {
    stop(
      "At least two samples are required to estimate nested sample variance.",
      call. = FALSE
    )
  }
  if (n_cell_types < 2L) {
    stop(
      "At least two cell types are required to fit the hierarchical model.",
      call. = FALSE
    )
  }

  cell_types <- colnames(p_estimated)
  samples <- rownames(p_estimated)
  ref_grid <- expand.grid(
    sample = samples,
    cell_type = cell_types,
    stringsAsFactors = FALSE
  )
  ref_grid$source <- "reference"
  ref_grid$y <- p_obs[ref_grid$cell_type]

  est_idx <- cbind(
    match(ref_grid$sample, samples),
    match(ref_grid$cell_type, cell_types)
  )
  est_grid <- ref_grid
  est_grid$source <- "estimated"
  est_grid$y <- p_estimated[est_idx]

  model_data <- rbind(ref_grid, est_grid)
  model_data$source <- factor(
    model_data$source,
    levels = c("reference", "estimated")
  )
  model_data$cell_type <- factor(model_data$cell_type, levels = cell_types)
  model_data$sample <- factor(model_data$sample, levels = samples)
  model_data$cell_type_sample <- interaction(
    model_data$cell_type,
    model_data$sample,
    drop = TRUE
  )

  fit <- lme4::lmer(
    y ~ source + (1 | cell_type) + (1 | cell_type_sample),
    data = model_data,
    REML = identical(method, "REML"),
    control = lme4::lmerControl(check.conv.singular = "ignore")
  )

  variance_components <- .extract_hrrmse_variances(fit)
  biological_variance <- variance_components["population"] +
    variance_components["sample"]

  if (!is.finite(biological_variance) || biological_variance <= 0) {
    stop(
      "Biological signal variance must be positive to compute hrRMSE.",
      call. = FALSE
    )
  }

  hrrmse <- sqrt(variance_components["residual"] / biological_variance)
  attr(hrrmse, "variance_components") <- variance_components
  attr(hrrmse, "method") <- method
  hrrmse
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

.prepare_hrrmse_compositions <- function(
    p_obs,
    p_estimated,
    trim_shared_zeros = TRUE
) {
  if (!is.numeric(p_obs) || !is.vector(p_obs)) {
    stop("`p_obs` must be a named numeric vector.", call. = FALSE)
  }
  if (is.null(names(p_obs)) || any(names(p_obs) == "")) {
    stop("`p_obs` must be a named vector.", call. = FALSE)
  }

  if (is.list(p_estimated) && !is.data.frame(p_estimated)) {
    sample_names <- names(p_estimated)
    if (is.null(sample_names) || any(sample_names == "")) {
      stop(
        "List `p_estimated` entries must be named by sample.",
        call. = FALSE
      )
    }
    cell_types <- names(p_obs)
    p_estimated <- do.call(
      rbind,
      lapply(p_estimated, function(sample_comp) {
        if (!identical(names(sample_comp), cell_types)) {
          stop(
            "Each estimated composition must match `p_obs` names.",
            call. = FALSE
          )
        }
        sample_comp
      })
    )
    rownames(p_estimated) <- sample_names
  } else if (is.data.frame(p_estimated)) {
    p_estimated <- as.matrix(p_estimated)
  } else if (!is.matrix(p_estimated)) {
    stop(
      "`p_estimated` must be a matrix, data frame, or named list.",
      call. = FALSE
    )
  }

  if (is.null(colnames(p_estimated)) || any(colnames(p_estimated) == "")) {
    stop("`p_estimated` columns must be named cell types.", call. = FALSE)
  }
  if (!identical(names(p_obs), colnames(p_estimated))) {
    stop(
      "`p_obs` names must match `p_estimated` column names in the same order.",
      call. = FALSE
    )
  }
  if (nrow(p_estimated) < 2L) {
    stop("`p_estimated` must contain at least two samples.", call. = FALSE)
  }

  if (trim_shared_zeros) {
    zero_in_all_estimates <- apply(p_estimated, 2, function(col) {
      all(col == 0)
    })
    keep_idx <- !(p_obs == 0 & zero_in_all_estimates)
    p_obs <- p_obs[keep_idx]
    p_estimated <- p_estimated[, keep_idx, drop = FALSE]
  }
  if (length(p_obs) < 2L) {
    stop(
      "At least two cell types must remain after trimming shared zeros.",
      call. = FALSE
    )
  }

  p_obs <- normalise_cell_estimates(p_obs)
  p_estimated <- apply(
    p_estimated,
    1,
    normalise_cell_estimates,
    simplify = FALSE
  )
  p_estimated <- do.call(rbind, p_estimated)
  colnames(p_estimated) <- names(p_obs)

  if (is.null(rownames(p_estimated))) {
    rownames(p_estimated) <- paste0("sample_", seq_len(nrow(p_estimated)))
  }

  list(
    p_obs = p_obs,
    p_estimated = p_estimated
  )
}

#' Extract variance components from an hrRMSE mixed model
#'
#' @description
#' Retrieves the three variance components estimated by
#' [lme4::lmer()] for the hrRMSE hierarchical model.
#'
#' @details
#' For a fitted model of the form
#' \deqn{
#' y_{ijm} = \mu + \beta\, I(m = \text{estimated}) + u_j + v_{ij} + \varepsilon_{ijm},
#' }
#' the random-effect standard deviations returned by [lme4::VarCorr()] are
#' squared to obtain
#' \deqn{\hat\sigma^2_{\mathrm{population}} = \mathrm{Var}(u_j),}
#' \deqn{\hat\sigma^2_{\mathrm{sample}} = \mathrm{Var}(v_{ij}),}
#' from the `(1 | cell_type)` and `(1 | cell_type:sample)` terms,
#' respectively. The residual variance is
#' \deqn{\hat\sigma_e^2 = \sigma(\mathrm{fit})^2,}
#' where `\sigma(\mathrm{fit})` is given by [stats::sigma()].
#'
#' @param fit A fitted `lmerMod` object from [lme4::lmer()].
#'
#' @return A named numeric vector with elements `population`, `sample`, and
#'   `residual`.
#'
#' @keywords internal
.extract_hrrmse_variances <- function(fit) {
  variance_cor <- lme4::VarCorr(fit)
  population_var <- attr(variance_cor$cell_type, "stddev")^2
  sample_var <- attr(variance_cor$cell_type_sample, "stddev")^2
  residual_var <- stats::sigma(fit)^2

  if (any(!is.finite(c(population_var, sample_var, residual_var)))) {
    stop("Could not extract finite variance components from the model.", call. = FALSE)
  }

  stats::setNames(
    c(population_var, sample_var, residual_var),
    c("population", "sample", "residual")
  )
}
