#' Evaluate Principal Component Functions from mGSFPCA
#'
#' Evaluates the principal component functions (eigenfunctions) from mGSFPCA()
#' object at specified time points.
#'
#' @param mGSFPCA_obj A list object returned by the \code{mGSFPCA} function,
#' containing the estimated principal components, basis functions, and other
#' model parameters.
#' @param eval_pts Numeric vector specifying the time points at which to
#' evaluate the principal component functions.
#'
#' @return A matrix of dimension \code{length(eval_pts) x p}, where \code{p} is
#' the number of principal components, containing the evaluated principal
#' component functions at the specified \code{eval_pts}.
#'
#' @details
#' The \code{eval_mGSFPCA} function evaluates the principal component functions
#' (eigenfunctions) from an \code{mGSFPCA} object at the provided \code{eval_pts}.
#'
#' @import fda
#' @export
eval_mGSFPCA <- function(mGSFPCA_obj, eval_pts = NULL) {

  k <- mGSFPCA_obj$pars$k
  p <- mGSFPCA_obj$pars$p

  if (is.null(eval_pts)) {
    eval_pts <- mGSFPCA_obj$pars$evalGrid
  } else if (setequal(eval_pts, mGSFPCA_obj$pars$workGrid)) {
    eval_pts <- mGSFPCA_obj$pars$workGrid
  } else { # for external input
    eval_pts <- sort(eval_pts)
    aT <- mGSFPCA_obj$pars$range[1]; bT <- mGSFPCA_obj$pars$range[2]
    if (all(eval_pts >= aT) && all(eval_pts <= bT)) {
      eval_pts <- (eval_pts - aT) / (bT - aT)
    } else {
      stop("eval_pts outside of observation point range")
    }
  }

  cvec <- mGSFPCA_obj$pars$coeffs
  basis <- mGSFPCA_obj$pars$eigBasis
  npts <- length(eval_pts)

  if (npts > p) {
    Phi <- get_UDUT(cvec, basis, eval_pts, p, k)$U*sqrt(npts)
  } else {
    D_grid <- seq(0,1,0.001) #Evaluate at fine grid
    U <- get_UDUT(cvec, basis, D_grid, p, k)$U*sqrt(length(D_grid))
    B <- fda::eval.basis(D_grid, basis)
    BtB <- (t(B)%*%(B))*(D_grid[2]-D_grid[1])
    R.inv <- solve(chol(BtB))
    B_orth <- (B) %*% R.inv
    beta <- solve((t(B_orth)%*%(B_orth))) %*% t(B_orth) %*% U

    eval_B <- fda::eval.basis(eval_pts, basis)
    eval_B_orth <- eval_B %*% R.inv
    Phi <- eval_B_orth %*% beta

  }

  return(Phi)
}
