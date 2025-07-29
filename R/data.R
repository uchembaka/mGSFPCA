
#' @title Simulated B-Spline Functional Data
#'
#' @description
#' A simulated dataset containing true and observed functional data based on
#' B-spline functions
#'
#' @format A list with 6 components:
#' \describe{
#'   \item{X}{A numeric matrix of dimension 100 x 51, representing the true
#'   functional data for 100 subjects evaluated at 51 time points.}
#'   \item{Y}{A numeric matrix of dimension 617 x 3, representing observed data
#'   with the following columns:
#'     \describe{
#'       \item{ID}{Integer, subject identifiers (1 to 100).}
#'       \item{time}{Numeric, observation time points in the range (0, 1).}
#'       \item{measurement}{Numeric, observed values at the corresponding time
#'       points.}
#'     }
#'   }
#'   \item{eigFun}{A numeric matrix of dimension 5 x 51, representing the true
#'   eigenfunctions (principal components) for 5 components evaluated on a
#'   51-point grid.}
#'   \item{eigVal}{A numeric vector of length 5, containing the true eigenvalues
#'    associated with the eigenfunctions, in descending order.}
#'   \item{nbasis}{An integer, the true number of B-spline basis functions used
#'   to generate the data, equal to 10.}
#'   \item{sig}{A numeric value, the true error variance (sigma^2) of the
#'   observation noise, equal to 0.167.}
#' }
#'
#' @examples
#' # example code
#'  bspline_sim
#'
#' @source
#' Simulated data based on:
#' Peng, J. and Paul, D. (2009). A geometric approach to maximum likelihood
#' estimation of the functional principal components from sparse longitudinal
#' data. Journal of Computational and Graphical Statistics.
#' \url{http://www.tandfonline.com/doi/abs/10.1198/jcgs.2009.08011}
#'
#' @name bspline_sim
NULL


#' @title Simulated Sparse Multivariate Functional Data
#'
#' @description
#' A simulated dataset for testing sparse multivariate functional principal
#' component analysis (e.g., \code{spMultFPCA}). It contains observed data for
#' three functional variables, true covariance, true trajectories,
#' eigenfunctions, eigenvalues, scores, and error variance for 100 subjects
#' evaluated over a 100-point grids.
#'
#' @format A list with 6 components:
#' \describe{
#'   \item{obs_data}{A list of 3 data frames, each representing observed data
#'   for one functional variable:
#'     \describe{
#'       \item{y1}{A data frame with 482 observations and 3 variables:
#'         \describe{
#'           \item{subj}{Integer, subject identifiers (1 to 100).}
#'           \item{argvals}{Numeric, observation time points in the range (0, 1).}
#'           \item{y}{Numeric, observed values for the first functional variable.}
#'         }
#'       }
#'       \item{y2}{A data frame with 489 observations and 3 variables:
#'         \describe{
#'           \item{subj}{Integer, subject identifiers (1 to 100).}
#'           \item{argvals}{Numeric, observation time points in the range (0, 1).}
#'           \item{y}{Numeric, observed values for the second functional variable.}
#'         }
#'       }
#'       \item{y3}{A data frame with 510 observations and 3 variables:
#'         \describe{
#'           \item{subj}{Integer, subject identifiers (1 to 100).}
#'           \item{argvals}{Numeric, observation time points in the range (0, 1).}
#'           \item{y}{Numeric, observed values for the third functional variable.}
#'         }
#'       }
#'     }
#'   }
#'   \item{True_Cov}{A numeric matrix of dimension 300 x 300, representing the
#'   true multivariate covariance matrix across all variables evaluated on a
#'   300-point grid.}
#'   \item{True_obs}{A list of 3 numeric matrices, each of dimension 100 x 100,
#'   representing the true trajectories for 100 subjects for each of the three
#'   functional variables.}
#'   \item{True_Eigs}{A list of 2 components:
#'     \describe{
#'       \item{funcs}{A numeric matrix of dimension 300 x 300, representing the
#'       true multivariate eigenfunctions evaluated on a 300-point grid.}
#'       \item{vals}{A numeric vector of length 300, containing the true
#'       multivariate eigenvalues.}
#'     }
#'   }
#'   \item{True_scores}{A numeric matrix of dimension 100 x 9, containing the
#'   true principal component scores for 100 subjects across 9 multivariate
#'   components.}
#'   \item{sigma}{A numeric value, the true error variance (sigma^2) of the
#'   observation noise.}
#' }
#'
#' @examples
#' # example code
#'  sp_mult_sim
#'
#' @source
#' Simulated data based on:
#' Li, C., Xiao, L., and Luo, S. (2020). Fast covariance estimation for multi-
#' variate sparse functional data. Stat, 9(1):e245.
#' \url{https://onlinelibrary.wiley.com/doi/10.1002/sta4.245}
#'
#' @name sp_mult_sim
NULL
