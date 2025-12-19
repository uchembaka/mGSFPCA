#' Evaluate Principal Component Functions from mGSFPCA
#'
#' Evaluates the principal component functions (eigenfunctions) from mGSFPCA()
#' object at specified time points.
#'
#' @param mGSFPCA_obj A list object returned by the \code{mGSFPCA} function,
#' containing the estimated principal components, basis functions, and other
#' model parameters.
#' @param eval_pts Numeric vector specifying the time points at which to
#' evaluate the principal component functions. The time points should either be
#' in the original data range or \eqn{\in [0, 1]}, in which case the output is
#' equivalent to calling fda::eval.fd on mGSFPCA_obj$eigenfunction.
#' @param matrix_orth Logical indicating whether to enforce strict orthogonality
#' of the evaluated eigenfunctions \eqn{\Phi^\top \Phi = I_p} under standard
#' matrix multiplication (discrete inner product). Diagnostic comparison plots
#' are displayed. Requires `length(eval_pts) > p`.
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
eval_mGSFPCA <- function(mGSFPCA_obj, eval_pts = NULL, matrix_orth = FALSE) {

  k <- mGSFPCA_obj$pars$k
  p <- mGSFPCA_obj$pars$p

  if (is.null(eval_pts)) {
    return(mGSFPCA_obj$Phi)
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
  evpts <- binData(eval_pts, mGSFPCA_obj$pars$evalGrid)

  if (matrix_orth) {
    if (npts < p) stop("The number of evaluation points (eval_pts) < p while matrix_orth set to TRUE")
    cat('See plot to compare eigenfunctions evaluation at eval_pts to original grid\n')
    Phi <- get_UDUT(cvec, basis, evpts, p, k)$U
    for (i in 1:p) {
      plot(mGSFPCA_obj$pars$evalGrid, mGSFPCA_obj$Phi[,i],
           ylim = range(c(mGSFPCA_obj$Phi[,i], Phi[,i])), lwd = 2,
           type = 'l', ylab = expression(phi(t)), xlab = 't', main = paste0("Eigenfunction ", i))
      points(evpts, Phi[,i], col = 2, pch = 19)
      legend("top", legend=c("eval_pts"), pch = 19,
             col=c("red"), cex=0.8, box.lty=0)
      readline(prompt = "Press ENTER to see the next plot...")
    }
  } else {
    Phi <- sapply(1:p, function(k) fda::eval.fd(evpts, mGSFPCA_obj$eigenfunctions[[k]]))

    # # Similarly (As in the supplementary)
    # B <- fda::eval.basis(mGSFPCA_obj$pars$evalGrid, mGSFPCA_obj$pars$eigBasis)
    # int_BBt <- fda::inprod(mGSFPCA_obj$pars$eigBasis, mGSFPCA_obj$pars$eigBasis)
    # R.inv <- solve(chol(int_BBt))
    #
    # Orth_B <- B %*% R.inv
    #
    # tilde_C <- int_Bt_phi(Orth_B, mGSFPCA_obj$Phi, mGSFPCA_obj$pars$evalGrid)
    #
    # B_eval_pts<- fda::eval.basis(binData(evpts,mGSFPCA_obj$pars$workGrid), mGSFPCA_obj$pars$eigBasis)
    # Orth_eval_pts <- B_eval_pts %*% R.inv
    # Phi <- Orth_eval_pts %*% tilde_C

  }

  return(Phi)
}

