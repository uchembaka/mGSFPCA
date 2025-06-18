#' Compute Fitted Trajectories and Scores from mGSFPCA
#'
#' Estimates fitted trajectories and principal component scores for each subject using
#' the results from mGSFPCA().
#'
#' @param mGSFPCA_obj A list object returned by the \code{mGSFPCA} function, containing
#' the estimated principal components, mean function, and other model parameters.
#' @param eval_pts Numeric vector specifying the time points at which to evaluate the
#' fitted trajectories. If \code{NULL} (default), the evaluation grid from
#' \code{mGSFPCA_obj$pars$evalGrid} is used.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{xHat}: A matrix of dimension \code{n x length(eval_pts)}, where \code{n}
#'   is the number of subjects, containing the fitted trajectories for each subject
#'   evaluated at \code{eval_pts}.
#'   \item \code{Xi}: A matrix of dimension \code{n x p}, where \code{p} is the number
#'   of principal components, containing the principal component scores for each subject.
#' }
#'
#' @details
#' The \code{get_xHat} function computes the fitted trajectories and principal component
#' scores for each subject based on the output of the \code{mGSFPCA} function. It uses the
#' PACE (Principal Analysis by Conditional Expectation) method to estimate the scores.
#'
#' @references
#' Yao, F., Müller, H.-G., & Wang, J.-L. (2005). Functional data analysis for sparse
#' longitudinal data. Journal of the American Statistical Association, 100(470), 577–590.
#'
#' @import fda
#' @export
get_xHat <- function(mGSFPCA_obj, eval_pts = NULL) {

  if (is.null(eval_pts)) {
    eval_pts <- mGSFPCA_obj$pars$evalGrid
  }
  aT <- mGSFPCA_obj$pars$range[1]; bT <- mGSFPCA_obj$pars$range[2]
  if (all(eval_pts >= aT) && all(eval_pts <= bT)) {
    eval_pts <- (eval_pts - aT) / (bT - aT)
  } else {
    stop("eval_pts outside of observation point range")
  }

  k <- mGSFPCA_obj$pars$k
  p <- mGSFPCA_obj$pars$p
  data <- mGSFPCA_obj$pars$binData
  IDs <- unique(data[,1])
  n <- length(IDs)
  eval_mu <- fda::eval.basis(eval_pts, mGSFPCA_obj$pars$mu_basis) %*%
    mGSFPCA_obj$pars$mu_fdobj$fd$coefs

  eval_U <- as.matrix(eval_mGSFPCA(mGSFPCA_obj, eval_pts))
  est_U <- eval_mGSFPCA(mGSFPCA_obj, mGSFPCA_obj$pars$workGrid)
  D <- diag(mGSFPCA_obj$Lambda)
  est_Cov <- (est_U %*% D %*% t(est_U)) + mGSFPCA_obj$sig2 * diag(nrow(est_U))

  xHat_pace <- matrix(0, nrow = n, ncol = length(eval_pts))
  xi <- matrix(NA, nrow = n, ncol = p)

  for (i in 1:n) {
    iind <- which(data[,1] == IDs[i])
    ti <- data[iind,5]
    mi <- length(ti)
    if(length(ti) > 0) {
      Cyy <- topdm_cpp(est_Cov[ti, ti, drop = FALSE])
      # PACE estimator
      xii <- (D %*% t(est_U[ti, , drop = FALSE])) %*% solve(Cyy) %*% data[iind, 4]
      xhati <- eval_mu + eval_U %*% xii
      xHat_pace[i,] <- xhati
      xi[i,] <- xii
    }
  }

  list(xHat = xHat_pace, Xi = xi)
}
