#' Evaluate Principal Component Functions from spMultFPCA
#'
#' Evaluates the principal component functions (eigenfunctions) from
#' spMultFPCA() object at specified time points (eval_pts).
#'
#' @param spMultFPCA_obj A list object returned by the \code{spMultFPCA}
#' function, containing the estimated principal components, basis functions,
#' and other model parameters.
#' @param eval_pts Numeric vector specifying the time points at which to
#' evaluate the principal component functions.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{Ckk}: A list containing the eigenfunctions, scores, and curve
#'   prediction for each variable at eval_pts.
#'   \item \code{scores}: A matrix of multivariate principal component scores
#'   for each subject.
#'   \item \code{eigFunctions}: A list of matrices, each containing the
#'   multivariate eigenfunctions for a variable.
#'   \item \code{values}: A vector of multivariate eigenvalues.
#'   \item \code{vectors}: A matrix representing the eigenvectors associated
#'   with the combined univariate score vectors
#'   \item \code{npc}: The number of multivariate principal components retained.
#'   \item \code{Cov}: The estimated multivariate covariance matrix.
#'   \item \code{fullEigFun}: A matrix of concatenated multivariate
#'   eigenfunctions across all variables.
#' }
#'
#' @details
#' The \code{eval_spMultFPCA} function evaluates the principal component
#' functions (eigenfunctions) from an \code{spMultFPCA} object at the provided
#' \code{eval_pts}.
#'
#' @import fda
#' @export
eval_spMultFPCA <- function(spMultFPCA_obj, eval_pts = NULL) {
  p <- length(spMultFPCA_obj$Ckk)
  Ckk <- vector("list", p)
  bigXi <- NULL

  for (i in 1:p) {
    result2 <- get_xHat(spMultFPCA_obj$Ckk[[i]], eval_pts)
    CkkStr <- list(
      Phi = eval_mGSFPCA(spMultFPCA_obj$Ckk[[i]], eval_pts),
      Xi = result2$Xi,
      xHat = result2$xHat
    )

    Ckk[[i]] <- CkkStr
    bigXi <- cbind(bigXi, result2$Xi)
  }

  Z <- cov(bigXi)
  svd_result <- svd(Z)
  cm <- svd_result$u
  vm <- diag(svd_result$d)

  if (is.null(M_npc)) {
    elbow_result <- elbow(
      data = data.frame(
        x = 1:length(svd_result$d),
        y = svd_result$d/sum(svd_result$d)
      )
    )
    M_npc <- elbow_result$x_selected
  }

  psi <- vector("list", p)
  rho <- bigXi %*% cm

  fullEigFun <- NULL
  ed <- 0
  st <- 1

  for (i in 1:p) {
    Mi <- Ckk[[i]]$pars$p
    ed <- ed + Mi
    psi[[i]] <- Ckk[[i]]$Phi %*% cm[st:ed, ]
    fullEigFun <- rbind(fullEigFun, Ckk[[i]]$Phi %*% cm[st:ed, ])
    st <- ed + 1
  }

  # Prepare outputs
  outs <- list(
    Ckk = Ckk,
    scores = rho,
    eigFunctions = psi,
    values = diag(vm),
    vectors = cm,
    npc = M_npc,
    Cov = fullEigFun %*% vm %*% t(fullEigFun),
    fullEigFun = fullEigFun
  )

  return(outs)
}
