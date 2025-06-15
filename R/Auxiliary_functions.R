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


get_UDUT <- function(cvec, basis, evalGrid, p, k) {
  indB <- cumsum(c(1, rep(k, p)))
  M <- length(evalGrid)
  U <- matrix(0, nrow = M, ncol = p)
  grid <- evalGrid
  
  B <- fda::eval.basis(grid, basis)
  
  CU <- cvec[1:(k * p)]
  CDCSig2 <- cvec[indB[p + 1]:length(cvec)]
  CD <- CDCSig2[1:(length(CDCSig2) - 1)]
  Csig2 <- cvec[length(cvec)]
  
  C <- matrix(CU, nrow = k, ncol = p)  
  U <- B %*% C  
  D <- exp(CD)
  
  ord <- order(D, decreasing = TRUE)
  D <- diag(D[ord], nrow = length(D))
  U_ortho <- pracma::gramSchmidt(U)$Q
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



