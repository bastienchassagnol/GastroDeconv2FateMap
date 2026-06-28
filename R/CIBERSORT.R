#' Main functions
#'
#' The Main function of CIBERSORT
#' @param sig_matrix  sig_matrix file path to gene expression from isolated cells, or a matrix of expression profile of cells.
#'
#' @param mixture_file mixture_file file path to heterogenous mixed expression file, or a matrix of heterogenous mixed expression
#'
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @param maxSize maximum size for the computation, to be passed to the future.global.maxSize. Default to 500 (in MB)
#' @import utils
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom stats sd
#' @export
#' @examples
#' \dontrun{
#'   ## example 1
#'   sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#'   mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
#'   results <- cibersort(sig_matrix, mixture_file)
#'   ## example 2
#'   data(LM22)
#'   data(mixed_expr)
#'   results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr)
#' }
cibersort <- function(sig_matrix, mixture_file, maxSize = 500, QN = TRUE) {
  #read in data
  if (is.character(sig_matrix)) {
    X <- read.delim(
      sig_matrix,
      header = T,
      sep = "\t",
      row.names = 1,
      check.names = F
    )
    X <- data.matrix(X)
  } else {
    X <- sig_matrix
  }

  if (is.character(mixture_file)) {
    Y <- read.delim(
      mixture_file,
      header = T,
      sep = "\t",
      row.names = 1,
      check.names = F
    )
    Y <- data.matrix(Y)
  } else {
    Y <- mixture_file
  }

  #order
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]

  #anti-log if max < 50 in mixture file
  if (max(Y) < 50) {
    Y <- 2^Y
  }

  #quantile normalization of mixture file
  if (QN == TRUE) {
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX, ]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY, ]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  #print(nulldist)

  header <- c('Mixture', colnames(X), "P-value", "Correlation", "RMSE")
  #print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- NA_real_

  #iterate through mixtures
  while (itor <= mixtures) {
    y <- Y[, itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y, maxSize)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #print output
    out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
    if (itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }

    itor <- itor + 1
  }

  #return matrix object containing all results
  obj <- rbind(header, output)
  obj <- obj[, -1]
  obj <- obj[-1, ]
  obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X), "P-value", "Correlation", "RMSE")
  obj
}


#' Configure a furrr parallel backend for CIBERSORT
#'
#' @param nThreads Number of parallel workers. Defaults to `availableCores() - 2`.
#' @param maxSize Maximum global object size (MB) for future exports.
#'
#' @keywords internal
.enable_parallel_cibersort <- function(nThreads = NULL, maxSize = 500) {
  if (is.null(nThreads)) {
    nThreads <- max(1L, future::availableCores() - 2L)
  }
  future::plan(future::multisession, workers = nThreads)
  options(future.globals.maxSize = maxSize * 1024^2)
}

#' Core algorithm
#'
#' The core algorithm of CIBERSORT which uses nu-SVR (`e1071::svm`).
#'
#' @param X Cell-specific gene expression matrix.
#' @param y Mixed expression vector for one sample.
#' @param maxSize Maximum global object size (MB) for parallel workers.
#'
#' @return A list with `w` (weights), `mix_rmse`, and `mix_r`.
#'
#' @importFrom furrr future_map
#' @importFrom future availableCores
#' @importFrom stats cor
#' @importFrom e1071 svm
#' @export
CoreAlg <- function(X, y, maxSize = 500) {
  #try different values of nu
  svn_itor <- 3

  res <- function(i) {
    if (i == 1) {
      nus <- 0.25
    }
    if (i == 2) {
      nus <- 0.5
    }
    if (i == 3) {
      nus <- 0.75
    }
    model <- e1071::svm(
      X,
      y,
      type = "nu-regression",
      kernel = "linear",
      nu = nus,
      scale = FALSE
    )
    model
  }

  .enable_parallel_cibersort(maxSize = maxSize)

  if (Sys.info()[["sysname"]] == "Windows") {
    out <- furrr::future_map(1:svn_itor, res)
  } else {
    if (svn_itor <= future::availableCores() - 2L) {
      .enable_parallel_cibersort(nThreads = svn_itor, maxSize = maxSize)
    } else {
      .enable_parallel_cibersort(maxSize = maxSize)
    }
    out <- furrr::future_map(1:svn_itor, res)
  }

  nusvm <- rep(0, svn_itor)
  corrv <- rep(0, svn_itor)

  #do cibersort
  t <- 1
  while (t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights < 0)] <- 0
    w <- weights / sum(weights)
    u <- sweep(X, MARGIN = 2, w, '*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- stats::cor(k, y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q < 0)] <- 0
  w <- (q / sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
}
