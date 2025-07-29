#' Sparse Multivariate Functional Principal Component Analysis
#'
#' Performs sparse multivariate functional principal component analysis
#' (spMultFPCA) on a list of functional datasets, estimating multivariate
#' principal components using the modified Gram-Schmidt Functional Principal
#' Component Analysis (mGSFPCA) for each variable and combining results via
#' singular value decomposition (SVD).
#'
#' @param dataCell A list of matrices or data frames, each with three columns:
#' ID (subject identifier), time (observation time points), and value
#' (observed values) for each functional variable.
#' @param r A list of integer vectors specifying candidate numbers of principal
#' components for each variable. If \code{NULL}, defaults to \code{3:5} for
#' each variable.
#' @param G A list of integer vectors specifying candidate numbers of basis
#' functions for each variable. If \code{NULL}, defaults to \code{c(5, 10, 15)}
#' for each variable.
#' @param basis_type Character string or vector specifying the type of basis
#' functions (\code{"bspline"} or \code{"fourier"}) for each variable. If a
#' single string, it is applied to all variables. Defaults to \code{"bspline"}.
#' @param nRegGrid Integer or list of integers specifying the number of points
#' in the regular grid for evaluation for each variable. If a single integer,
#' it is applied to all variables. If
#' \code{NULL}, defaults to \code{50} for each variable.
#' @param grid_range A list of numeric vectors of length 2 specifying the
#' observation range \code{[aT, bT]} for each variable. If a single range is
#' provided, it is applied to all variables. If \code{NULL}, ranges are
#' determined from the data.
#' @param mu_nbasis Integer or list of integers specifying the number of basis
#' functions for the mean function of each variable. If a single integer, it is
#' applied to all variables. If \code{NULL}, defaults to \code{15} for each
#' variable.
#' @param M_npc Integer specifying the number of multivariate principal
#' components to retain. If \code{NULL}, determined automatically using the
#' elbow method.
#' @param maxit Integer specifying the maximum number of iterations for the
#' optimization algorithm. Default is \code{500}.
#' @param optim_tol Numeric specifying the tolerance for optimization
#' convergence. Default is \code{1e-5}.
#' @param optim_trace Integer specifying the level of tracing information from
#' the optimization algorithm. Default is \code{0} (no tracing).
#' @param bin_size Integer or list of integers specifying the number of bins for
#'  data binning for each variable. If a single integer, it is applied to all
#'  variables.
#' If \code{NULL}, defaults to \code{51} for each variable.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{Ckk}: A list of results from \code{mGSFPCA} for each variable.
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
#' @examples
#' sp_mult_sim
#' tmp <- spMultFPCA(sp_mult_sim$obs_data,
#'                   r = rep(list(3:5), 3),
#'                   G = rep(list(5:10), 3),
#'                   nRegGrid = 100)
#' Metrics::rmse(tmp$Cov, sp_mult_sim$True_Cov)
#' grid <- seq(0, 1, length = 300)
#' plot(grid, sp_mult_sim$True_Eigs$funcs[,1], type ='l', ylab = 'phi')
#' lines(grid, -tmp$fullEigFun[,1], col = 2)
#'
#' @details
#' The \code{spMultFPCA} function performs multivariate functional principal
#' component analysis on sparse functional data by applying \code{mGSFPCA} to
#' each variable in \code{dataCell}.
#'
#' @references
#' Happ, C., & Greven, S. (2018). Multivariate functional principal component
#' analysis for data observed on different dimensional domains. Journal of the
#' American Statistical Association, 113(522), 649–659.
#'
#' @import fda
#' @export
spMultFPCA <- function(dataCell,
                       r = NULL,
                       G = NULL,
                       basis_type = "bspline",
                       nRegGrid = 51,
                       grid_range = NULL,
                       mu_nbasis = 15,
                       M_npc = NULL,
                       maxit = 500,
                       optim_tol = 1e-5,
                       optim_trace = 0,
                       bin_size = 51) {

  p <- length(dataCell)

  # Rank
  if (is.null(r)) {
    r <- rep(list(3:5), p)
  } else if (!is.list(r)) {
    stop("r should be a list")
  } else if (length(r) != p) {
    stop("length of r list should match length of variables")
  }

  # Number of basis
  if (is.null(G)) {
    G <- rep(list(c(5, 10, 15)), p)
  } else if (!is.list(G)) {
    stop("G should be a list")
  } else if (length(G) != p) {
    stop("length of G list should match length of variables")
  }

  # basis_type
  if (is.null(basis_type)) {
    basis_type <- rep("bspline", p)
  } else if (length(basis_type) == 1) {
    basis_type <- rep(basis_type, p)
  }

  # bin_size
  if (is.null(bin_size)) {
    bin_size <- vector("list", p)
  } else if (is.numeric(bin_size) && length(bin_size) == 1) {
    bin_size <- rep(list(bin_size), p)
  }

  # nRegGrid
  if (is.null(nRegGrid)) {
    nRegGrid <- vector("list", p)
  } else if (is.numeric(nRegGrid) && length(nRegGrid) == 1) {
    nRegGrid <- rep(list(nRegGrid), p)
  }

  # grid_range
  if (is.null(grid_range)) {
    grid_range <- vector("list", p)
  } else if (!is.list(grid_range)) {
    stop("grid_range should be a list")
  } else if (length(grid_range) == 1) {
    if (length(grid_range[[1]]) != 2) {
      stop("Range should be two values")
    } else {
      grid_range <- rep(grid_range, p)
    }
  }

  # mu_nbasis
  if (is.null(mu_nbasis)) {
    mu_nbasis <- vector("list", p)
  } else if (is.numeric(mu_nbasis) && length(mu_nbasis) == 1) {
    mu_nbasis <- rep(list(mu_nbasis), p)
  }

  Ckk <- vector("list", p)
  bigXi <- NULL

  for (i in 1:p) {
    message(paste0("C_", i))
    datai <- dataCell[[i]]
    result <- mGSFPCA(data = datai, p = r[[i]], k = G[[i]],
                      basis_type = basis_type[i], bin_size = bin_size[[i]],
                      nRegGrid = nRegGrid[[i]], init_coeff = NULL,
                      obs_range = grid_range[[i]], mu_nbasis = mu_nbasis[[i]],
                      maxit = maxit, optim_tol = optim_tol,
                      optim_trace = optim_trace)
    result2 <- get_xHat(result)
    CkkStr <- list(
      Cov = result$Phi %*% diag(result$Lambda) %*% t(result$Phi),
      Phi = result$Phi,
      Lambda = result$Lambda,
      mu = result$mu,
      sig2 = result$sig2,
      pars = result$pars,
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
