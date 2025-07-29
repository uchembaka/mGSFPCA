#Auxiliary Functions
init_coeffs <- function(p, k, ini = NULL) {

  if (is.null(ini)) {
    ncoefs <- k*p + p + 1
    set.seed(1)
    ini <- runif(ncoefs, min = -1, max = 1)
    # set.seed(1)
    # ini <- rnorm(ncoefs)
    ini[(k*p+1) : (ncoefs-1)] <- log(rev(seq(1, 100, length.out = p)))
    ini[ncoefs] <- log(0.1)
  }

  indB <- cumsum(c(1, rep(k, p)))
  cvec0 <- ini[1:(k * p + p + 1)]
  cvec0[indB[p + 1]:(length(cvec0) - 1)] <- ini[(k*p+1):(length(ini)-1)]
  cvec0[length(cvec0)] <- ini[length(ini)]

  return(cvec0)
}

run_optimization <- function(p_val, k_val, ini, SiCell, newData, est_pts, basis_type, bin_size, optim_control) {
  if (identical(basis_type, "bspline")) {
    basis <- fda::create.bspline.basis(c(0, 1), k_val)
  } else {
    basis <- fda::create.fourier.basis(c(0, 1), k_val)
  }
  B <- fda::eval.basis(est_pts, basis)

  # Initialize cvec0
  if (is.null(ini)) {
    cvec0 <- init_coeffs(p_val, k_val, ini)

    opt_init <- optim(
      par = cvec0,
      fn = NLLKf,
      gr = NLLKg,
      SiCell = SiCell,
      data = newData,
      p = p_val,
      k = k_val,
      B = B,
      estGrid = est_pts,
      method = "BFGS",
      control = list(maxit = 3)
    )
    cvec0 <- opt_init$par
  } else {
    cvec0 <- init_coeffs(p_val, k_val, ini)
  }

  n = length(SiCell)
  ncoefs <- p_val * k_val + p_val + 1
  nobs <- nrow(newData)
  maxPars <- p_val * length(est_pts) * p_val + 1
  M <- length(est_pts)

  opt_result <- tryCatch(
    {
      optim_result <- optim(
        par = cvec0,
        fn = NLLKf,
        gr = NLLKg,
        SiCell = SiCell,
        data = newData,
        p = p_val,
        k = k_val,
        B = B,
        estGrid = est_pts,
        method = "BFGS",
        control = optim_control
      )

      # Compute AICs
      gcv_val <- n*optim_result$value + (k_val*p_val*p_val) + 1
      aic <- optim_result$par[ncoefs] + p_val*((n+M)/(n*M)) + log((n*M)/(n+M))

      list(par = optim_result$par, value = optim_result$value,
           gcv = gcv_val, aic = aic, convergence = optim_result$convergence)

    },
    error = function(e) {
      list(par = rep(-99, ncoefs), value = NaN, gcv = NaN, aic = NaN, convergence = 99)
    }
  )

  return(list(
    cvec = opt_result$par,
    gcv = gcv_val,
    basis = basis,
    convergence = opt_result$convergence,
    LLK = opt_result$value,
    aic = aic
  ))
}

extract_results <- function(results, indices, GCVkp, cvecCell, basisCell, convergenceCell, LLKCell, aicCell) {
  for (i in indices) {
    if (!is.null(results[[i]])) {
      cvecCell[[i]] <- results[[i]]$cvec
      GCVkp[i] <- results[[i]]$gcv
      basisCell[[i]] <- results[[i]]$basis
      convergenceCell[i] <- results[[i]]$convergence
      LLKCell[i] <- results[[i]]$LLK
      aicCell[i] <- results[[i]]$aic
    }
  }
  return(list(GCVkp = GCVkp, cvecCell = cvecCell, basisCell = basisCell, convergenceCell = convergenceCell, LLKCell = LLKCell, aicCell = aicCell))
}


# Main function
doGCV <- function(p, k, binData, gcvData, basis_type, ini, bin_size, optim_control, all_kp) {

  kp1 <- length(k) == 1 & length(p) == 1
  if (kp1) {
    est_pts <- binData$est_pts
    SiCell <- binData$SiCell
    newData <- binData$newData
  } else {
    est_pts <- gcvData$est_pts
    SiCell <- gcvData$SiCell
    newData <- gcvData$newData
  }


  # Initialize results list
  results <- vector("list", 0)

  if (all_kp) {

    if (length(p) > 1 && length(k) > 1) {
      cat("Selecting k & p\n")
    } else if (length(p) > 1) {
      cat("Selecting p\n")
    } else if (length(k) > 1) {
      cat("Selecting k\n")
    }
    # Evaluate all combinations of p and k
    grid <- expand.grid(p = p, k = k)
    p <- grid$p
    k <- grid$k
    use <- k >= p
    p <- p[use]
    k <- k[use]

    GCVkp <- numeric(length(p))
    cvecCell <- vector("list", length(p))
    basisCell <- vector("list", length(p))
    convergenceCell <- numeric(length(p))
    LLKCell <- numeric(length(p))
    aicCell <- numeric(length(p))

    # Loop over valid p, k combinations
    for (i in 1:length(p)) {
      results[[i]] <- run_optimization(
        p_val = p[i],
        k_val = k[i],
        ini = ini,
        SiCell = SiCell,
        newData = newData,
        est_pts = est_pts,
        basis_type = basis_type,
        bin_size = bin_size,
        optim_control = optim_control
      )
    }

    # Extract results
    extracted <- extract_results(results, 1:length(p), GCVkp, cvecCell, basisCell, convergenceCell, LLKCell, aicCell)
    GCVkp <- extracted$GCVkp
    cvecCell <- extracted$cvecCell
    basisCell <- extracted$basisCell
    convergenceCell <- extracted$convergenceCell
    LLKCell <- extracted$LLKCell
    aicCell <- extracted$aicCell

  } else {
    # Stepwise optimization: first optimize k with pmax, then p with optimal k
    if (length(k) > 1) {
      pmax <- p[length(p)]
      GCVk <- numeric(length(k))
      cvecCell <- vector("list", length(k))
      basisCell <- vector("list", length(k))
      convergenceCell <- numeric(length(k))
      LLKCell <- numeric(length(k))
      aicCell <- numeric(length(k))

      cat("Selecting k with p =", pmax, "\n")

      for (i in 1:length(k)) {
        if (k[i] < pmax) next  # Skip if k < pmax (this line should never run)
        results[[i]] <- run_optimization(
          p_val = pmax,
          k_val = k[i],
          ini = ini,
          SiCell = SiCell,
          newData = newData,
          est_pts = est_pts,
          basis_type = basis_type,
          bin_size = bin_size,
          optim_control = optim_control
        )
      }

      # Extract results for step 1
      extracted <- extract_results(results, 1:length(k), GCVk, cvecCell, basisCell, convergenceCell, LLKCell, aicCell)
      GCVk <- extracted$GCVkp
      cvecCell <- extracted$cvecCell
      basisCell <- extracted$basisCell
      convergenceCell <- extracted$convergenceCell
      LLKCell <- extracted$LLKCell
      aicCell <- extracted$aicCell

      # Find optimal k
      valid_GCV <- GCVk >= 0 & !is.na(GCVk)
      if (all(is.na(GCVk))) stop("No valid GCV values found for k with p =", pmax)
      minGCV <- which.min(GCVk[valid_GCV])
      best_k <- k[valid_GCV][minGCV]
      GCVk_res <- cbind(p = pmax, k = k[k >= pmax], AIC = GCVk, converge = convergenceCell, LLK = LLKCell, AIC_p = aicCell)
    } else {
      GCVk <- c()
      GCVk_res <- c()
      best_k <- k
    }

    if (length(p) > 1) {
      # Step 2: Optimize p with fixed k = best_k
      results <- vector("list", 0)  # Reset results
      GCVp <- numeric(length(p))
      cvecCell_p <- vector("list", length(p))
      basisCell_p <- vector("list", length(p))
      convergenceCell_p <- numeric(length(p))
      LLKCell_p <- numeric(length(p))
      aicCell_p <- numeric(length(p))

      cat("Selecting p with k =", best_k, "\n")

      for (i in 1:length(p)) {
        if (p[i] > best_k) next  # Skip if p > best_k (this line should never run)
        results[[i]] <- run_optimization(
          p_val = p[i],
          k_val = best_k,
          ini = ini,
          SiCell = SiCell,
          newData = newData,
          est_pts = est_pts,
          basis_type = basis_type,
          bin_size = bin_size,
          optim_control = optim_control
        )
      }

      # Extract results for step 2
      extracted <- extract_results(results, 1:length(p), GCVp, cvecCell_p, basisCell_p, convergenceCell_p, LLKCell_p, aicCell_p)
      GCVp <- extracted$GCVkp
      cvecCell_p <- extracted$cvecCell
      basisCell_p <- extracted$basisCell
      convergenceCell_p <- extracted$convergenceCell
      LLKCell_p <- extracted$LLKCell
      aicCell_p <- extracted$aicCell

      # Find optimal p
      valid_GCV <- GCVp >= 0 & !is.na(GCVp)
      if (all(is.na(GCVp))) stop("No valid GCV values found for p with k =", best_k)
      minGCV <- which.min(GCVp[valid_GCV])
      best_p <- p[valid_GCV][minGCV]
      GCVp_res <- cbind(p = p[p <= best_k], k = best_k, AIC = GCVp, converge = convergenceCell_p, LLK = LLKCell_p, AIC_p = aicCell_p)
    } else {
      GCVp <- c()
      GCVp_res <- c()
      best_p <- p
    }

    # Combine results for GCV_res
    GCVkp <- c(GCVk, GCVp)
    cvecCell <- c(cvecCell, cvecCell_p)
    basisCell <- c(basisCell, basisCell_p)
    convergenceCell <- c(convergenceCell, convergenceCell_p)
    LLKCell <- c(LLKCell, LLKCell_p)
    aicCell <- c(aicCell, aicCell_p)
  }

  # Compute GCV results
  GCV_res <- if (all_kp) {
    cbind(p = p, k = k, AIC = GCVkp, convergence = convergenceCell, LLK = LLKCell, AIC_p = aicCell)
  } else {
    rbind(GCVk_res, GCVp_res)
  }

  # Find best p and k
  valid_GCV <- GCVkp >= 0 & !is.na(GCVkp)
  if (all(is.na(GCVkp))) stop("No valid GCV values found")
  minGCV <- which.min(GCVkp[valid_GCV])
  best_k <- if (all_kp) k[valid_GCV][minGCV] else best_k
  best_p <- if (all_kp) p[valid_GCV][minGCV] else best_p
  best_cvec <- cvecCell[valid_GCV][[minGCV]]
  basis <- basisCell[valid_GCV][[minGCV]]
  convergence <- convergenceCell[valid_GCV][[minGCV]]

  result <- run_optimization(
    p_val = best_p,
    k_val = best_k,
    ini = NULL,
    SiCell = binData$SiCell,
    newData = binData$newData,
    est_pts = binData$est_pts,
    basis_type = basis_type,
    bin_size = bin_size,
    optim_control = list(abstol = 1e-7, maxit = 1000)
  )
  best_cvec <- result$cvec

  return(list(
    k = best_k,
    p = best_p,
    cvec = best_cvec,
    GCV_res = GCV_res,
    basis = basis,
    convergence = convergence
  ))
}
