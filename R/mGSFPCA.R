## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
## usethis namespace: start
#' @useDynLib mGSFPCA, .registration = TRUE
## usethis namespace: end
NULL
#' Estimate Functional Principal Components from Sparse Data
#'
#' Functional principal component analysis with modified Gram-Schmidt
#' Orthornormalization and MLE.
#'
#' @param data A matrix or data frame with three columns: ID (subject ID),
#'  time (observation time points), and value (observed values).
#' @param p Integer vector specifying the candidate number of principal
#' components to consider. Default is 2:5.
#' @param k Integer vector specifying the candidate number of basis functions to
#'  consider. Default is c(5, 10, 15).
#' @param basis_type Character string specifying the type of basis functions to
#' use. Options are "bspline" (default) or "fourier".
#' @param maxit Integer specifying the maximum number of iterations for the
#' optimization algorithm. Default is 500.
#' @param optim_tol Numeric specifying the relative tolerance for the
#' optimization convergence (Based on optim()). Default is 1e-5.
#' @param optim_trace Integer specifying the level of tracing information from
#' the optimization algorithm (Based on optim()). Default is 0 (no tracing).
#' @param nRegGrid Integer specifying the number of points in the regular grid
#' for evaluation. Default is 101.
#' @param init_coeff Numeric vector of initial coefficients, or "LSQ".
#' If "LSQ", initialization is based on a least-squares estimate of the
#' covariance structure.
#' If NULL (default), coefficients are initialized randomly.
#' @param obs_range Numeric vector of length 2 specifying the observation range
#' (aT, bT). If NULL (default), it is set to c(0,1).
#' @param mu_nbasis Integer specifying the number of basis functions for the
#' mean function. Default is 15.
#' @param bin_size Integer specifying the number of bins for data binning.
#' Default is 101.
#' @param use_kp_grid Logical indicating whether to evaluate all combinations of
#'  p and k (TRUE, default) or use a stepwise model selection approach (FALSE).
#' @param alpha Numeric value \eqn{\in [0, 1]} controlling density-based weighting
#' of the binned observations. Default is 0.
#' @param skip_check Logical indicating whether to skip input validation checks.
#'  Default is FALSE.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item \code{Phi}: Matrix of estimated eigenfunctions.
#' \item \code{Lambda}: Vector of estimated eigenvalues.
#' \item \code{eigenfunctions}: List of eigenfunctions on \eqn{[0, 1]} as fd object.
#' \item \code{mu}: Vector of the estimated mean function evaluated on the grid.
#' \item \code{sig2}: Estimated variance of the error term.
#' \item \code{pars}: A list containing model parameters and results:
#' \itemize{
#' \item \code{p}: Optimal number of principal components.
#' \item \code{k}: Optimal number of basis functions.
#' \item \code{AIC}: Table of AIC values for different p and k combinations.
#' \item \code{coeffs}: Optimized coefficients.
#' \item \code{mu_fdobj}: Functional data object for the mean function.
#' \item \code{mu_basis}: Basis object for the mean function.
#' \item \code{workGrid}: Grid points used for estimation.
#' \item \code{evalGrid}: Evaluation grid points.
#' \item \code{eigBasis}: Basis object for the eigenfunctions.
#' \item \code{binData}: Binned data used for estimation.
#' \item \code{range}: Observation range.
#' \item \code{convergence}: Convergence status of the optimization
#' (Based on optim()). 0 indicates successful completion.
#' }
#' }
#'
#' @details
#' The mGSFPCA function implements functional principal component analysis using
#'  MLE with modified Gram-Schmidt orthornormalization.
#'
#' @examples
#' # Plot 3 random subjects
#' set.seed(111)
#' samp <- sample.int(100, 3)
#' # bspline_sim is automatically loaded
#' Y <- bspline_sim$Y; X <- bspline_sim$X
#' id <- unique(Y[,1])
#' grid <- seq(0, 1, length = 51)
#'
#' plot(Y[,2], Y[,3], type = 'n', ylab = 'Y(t)', xlab = 't')
#' for (i in 1:3) {
#'   points(Y[id == samp[i], 2], Y[id == samp[i], 3], col = i, pch = 19)
#'   lines(grid, X[samp[i], ], col = i, lwd = 2)
#' }
#'
#' # Estimate
#' fit <- mGSFPCA(Y, p = 3:5, k = c(5, 10, 15), basis_type = "bspline",
#' nRegGrid = 51, bin_size = 101)
#'
#' c(fit$pars$p, fit$pars$k)
#'
#' fit$Lambda
#'
#' plot(grid, fit$Phi[,1], type ='l', lty = 2)
#' lines(grid, bspline_sim$eigFun[1,], col = 2) # change sign if needed
#' min(Metrics::rmse(fit$Phi[,1], bspline_sim$eigFun[1,]),
#'  Metrics::rmse(-fit$Phi[,1], bspline_sim$eigFun[1,]))
#'
#' x_pred <- get_xHat(fit)
#' Metrics::rmse(x_pred$xHat, bspline_sim$X)
#'
#'
#' @references
#' Peng, J. and Paul, D. (2009). A geometric approach to maximum likelihood
#' estimation of the functional principal components from sparse longitudinal
#' data. Journal of Computational and Graphical Statistics.
#'
#' Yao, F., Müller, H.-G., & Wang, J.-L. (2005). Functional data analysis for
#' sparse longitudinal data. Journal of the American Statistical Association,
#' 100(470), 577–590.
#'
#' @import fda pracma Metrics stats
#' @export
#'


mGSFPCA <- function(data,
                    p = 2:5,
                    k = c(5,10,15),
                    basis_type = "bspline",
                    maxit = 500,
                    optim_tol = 1e-5,
                    optim_trace = 0,
                    nRegGrid = 51,
                    init_coeff = NULL,
                    obs_range = NULL,
                    mu_nbasis = 15,
                    bin_size = 101,
                    use_kp_grid = TRUE,
                    alpha = 0,
                    skip_check = TRUE){

  inputs <- input_check(data, p, k, basis_type, optim_trace,
                        maxit, optim_tol, bin_size, nRegGrid,
                        init_coeff, obs_range, mu_nbasis,
                        use_kp_grid, alpha, skip_check)

  comp_data <- data_prep(data = inputs$data, k = inputs$k, p = inputs$p,
                         obs_range = inputs$obs_range, aT = inputs$aT,
                         bT = inputs$bT, bin_size = inputs$bin_size,
                         nRegGrid = inputs$nRegGrid,
                         basis_type = inputs$basis_type, inputs$init_coeff,
                         inputs$mu_nbasis)

  opt_result <- fit_model(p = inputs$p, k = inputs$k,
                          binData = comp_data$binData,
                          gcvData = comp_data$gcvData,
                          basis_type = inputs$basis_type,
                          ini = comp_data$ini,
                          bin_size = comp_data$bin_size,
                          optim_control = inputs$optim_control,
                          all_kp = inputs$use_kp_grid,
                          alpha = inputs$alpha)


  # Extract parameters
  cvec <- opt_result$cvec
  p <- opt_result$p
  k <- opt_result$k
  basis <- opt_result$basis
  sig2 <- exp(cvec[length(cvec)])

  # Get functional components
  UD_result <- get_UDUT(cvec, basis, inputs$evalGrid, p, k)
  Phi <- UD_result$U
  Lambda <- UD_result$D/comp_data$bin_size
  mu <- fda::eval.basis(inputs$evalGrid, comp_data$mu_results$basis) %*%
    comp_data$mu_results$fdobj$fd$coefs

  # Return results
  list(
    Phi = Phi,
    Lambda = diag(Lambda),
    eigenfunctions = eigenFunctions_fd(Phi, basis),
    mu = mu,
    sig2 = sig2,
    pars = list(
      p = p,
      k = k,
      AIC = as.data.frame(opt_result$GCV_res[,1:5]),
      coeffs = cvec,
      mu_fdobj = comp_data$mu_results$fdobj,
      mu_basis = comp_data$mu_results$basis,
      workGrid = comp_data$binData$est_pts,
      evalGrid = inputs$evalGrid,
      eigBasis = basis,
      binData = comp_data$binData$newData,
      range = inputs$obs_range,
      convergence = opt_result$convergence
    )
  )
}
