binData <- function(X, grid) {
  X <- as.vector(X)
  grid <- as.vector(grid)
  Y <- numeric(length(X))

  # For each element in X, find closest value in grid
  for (i in seq_along(X)) {
    idx <- which.min(abs(grid - X[i]))
    Y[i] <- grid[idx]
  }

  return(Y)
}


regrid <- function(data, grid_range) {
  a <- grid_range[1]
  b <- grid_range[2]
  dp <- 0
  x <- sum(unique(data))
  x <- abs(x)

  # Find decimal places needed
  while (floor(x * 10^dp) != x * 10^dp) {
    dp <- dp + 1
  }

  data_fullGrid <- seq(a, b, by = 10^-dp)
  std_grid <- seq(0, 1, length.out = length(data_fullGrid))
  inds <- sapply(data, function(i) which(data_fullGrid == i)[1])
  return(std_grid[inds])
}


averageByIndex <- function(X, t) {
  unique_indices <- unique(t)
  Y <- matrix(0, nrow = length(unique_indices), ncol = ncol(X))
  firstOccurrence <- logical(length(t))

  for (i in 1:length(unique_indices)) {
    idx <- unique_indices[i]
    mask <- t == idx
    values <- X[mask, , drop = FALSE]
    Y[i, ] <- colMeans(values)
    firstOccurrence[which(mask)[1]] <- TRUE
  }

  return(list(Y = Y, firstOccurrence = firstOccurrence))
}


get_UDUT <- function(cvec, basis, evalGrid, p, k, wq = NULL) {
  indB <- cumsum(c(1, rep(k, p)))
  M <- length(evalGrid)
  U <- matrix(0, nrow = M, ncol = p)
  grid <- evalGrid

  B <- fda::eval.basis(grid, basis)

  if(is.null(wq)) {
    wq <- quad_weights(evalGrid)
  }

  CU <- cvec[1:(k * p)]
  CDCSig2 <- cvec[indB[p + 1]:length(cvec)]
  CD <- CDCSig2[1:(length(CDCSig2) - 1)]
  Csig2 <- cvec[length(cvec)]

  C <- matrix(CU, nrow = k, ncol = p)
  U <- B %*% C
  D <- exp(CD)

  ord <- order(D, decreasing = TRUE)
  D <- diag(D[ord], nrow = length(D))
  U_ortho <- mwGS(U, wq)
  U <- matrix(unlist(U_ortho), nrow = nrow(U))
  U <- U[, ord, drop = FALSE]

  list(U = U, D = D)
}



# get_mean function - Functional mean estimation
get_mean <- function(data, estGrid, type = "bspline",
                     knots = NULL, lambdaVec = NULL) {

  if (is.null(knots)) {
    knots <- length(estGrid)
  }

  if (is.null(lambdaVec)) {
    lambdaVec <- seq(1e-4, 1, length.out = 20)
  }

  grid <- estGrid
  raw_data <- data[,3]
  wk_ind <- data[,2]

  # Create basis
  if (type == "bspline") {
    basis <- fda::create.bspline.basis(rangeval = c(0,1),
                                       nbasis = knots)
  } else {
    basis <- fda::create.fourier.basis(rangeval = c(0,1),
                                       nbasis = knots)
  }

  gcv_vec <- numeric(length(lambdaVec))
  fdPar_list <- vector("list", length(lambdaVec))

  for (i in seq_along(lambdaVec)) {
    fdParobj <- fda::fdPar(fdobj = basis, Lfdobj = 2,
                           lambda = lambdaVec[i])
    smooth <- fda::smooth.basis(argvals = wk_ind, y = raw_data,
                                fdParobj = fdParobj)
    gcv_vec[i] <- smooth$gcv
    fdPar_list[[i]] <- fdParobj
  }

  best_idx <- which.min(gcv_vec)
  lam <- lambdaVec[best_idx]
  fdParobj <- fdPar_list[[best_idx]]
  fdobj <- fda::smooth.basis(argvals = wk_ind, y = raw_data, fdParobj = fdParobj)
  mu <- fda::eval.basis(grid, basis) %*% fdobj$fd$coefs
  fits <- fda::eval.basis(wk_ind, basis) %*% fdobj$fd$coefs

  list(mu = mu, fits = fits, basis = basis, fdobj = fdobj)
}

# Function to check if regridding is needed
chngRange <- function(X, r) {
  a <- r[1]; b <- r[2]
  min(X) > b || max(X) > b || min(X) < a || max(X) < a
}


quad_weights <- function(grid_points) {

  n_points <- length(grid_points)
  weights <- numeric(n_points)

  if (n_points == 1) {
    return(1)
  }

  # Interior points
  for (i in 2:(n_points - 1)) {
    weights[i] <- (grid_points[i + 1] - grid_points[i - 1]) / 2
  }

  # Endpoints
  weights[1] <- (grid_points[2] - grid_points[1]) / 2
  weights[n_points] <- (grid_points[n_points] - grid_points[n_points - 1]) / 2

  weights <- weights / sum(weights)
  return(weights)
}


density_weights <- function(observation_times, grid_points, bandwidth = "nrd0") {
  kde <- density(observation_times, bw = bandwidth,
                 from = min(grid_points), to = max(grid_points),
                 na.rm = TRUE)

  # Interpolate density at grid points
  densities <- approx(kde$x, kde$y, xout = grid_points)$y
  densities[is.na(densities)] <- min(densities, na.rm = TRUE)

  weights <- densities/ sum(densities)

  return(weights)
}

get_weights <- function(observation_times, grid_points, alpha = 0) {

  density_weights <- density_weights(observation_times, grid_points)

  quadrature_weights <- quad_weights(grid_points)

  combined_weights <- alpha * density_weights + (1 - alpha) * quadrature_weights
  combinde_weights <- combined_weights / sum(combined_weights)

  return(combined_weights)
}


eigenFunctions_fd <- function(Phi, basis) {
  grid <- seq(0,1, length = nrow(Phi))
  B <- fda::eval.basis(grid, basis)
  int_BBT <- fda::inprod(basis, basis)
  Gram_mat <- solve(chol(int_BBT))
  orthB <- B %*% Gram_mat
  CU <- Gram_mat %*% int_Bt_phi(orthB, Phi, grid)
  p <- ncol(Phi)
  return(lapply(1:p, function(i) fda::fd(CU[,i], basis)))
}


int_Bt_phi <- function(B, phi, x) {

  n1 <- nrow(B)
  n2 <- nrow(phi)

  k <- ncol(B)
  p <- ncol(phi)
  n <- n1

  mats <- array(0, dim = c(k, p, n))

  for (i in 1:n) {
    Bi   <- B[i, , drop = FALSE]
    phii <- phi[i, , drop = FALSE]

    mats[,,i] <- t(Bi) %*% phii
  }

  result <- matrix(0, k, p)

  for (i in 1:k) {
    for (j in 1:p) {
      result[i, j] <- pracma::trapz(x, mats[i, j, ])
    }
  }

  return(result)
}



int_phi_t_phi <- function(phi, x) {
  n <- nrow(phi)
  p <- ncol(phi)
  mats <- array(0, dim = c(p, p, n))

  for (i in 1:n) {
    phi_i <- phi[i, , drop = FALSE]
    mats[,,i] <- t(phi_i) %*% phi_i
  }
  result <- matrix(0, p, p)

  for (i in 1:p) {
    for (j in 1:p) {
      result[i, j] <- pracma::trapz(x, mats[i, j, ])
    }
  }

  result
}

