# Auxiliary Functions
init_coeffs <- function(p, k, ini = NULL) {

  ncoefs <- k * p + p + 1
  if (is.null(ini)) {
    set.seed(1)
    ini <- runif(ncoefs, min = -1, max = 1)
    ini[(k*p+1) : (ncoefs-1)] <- log(rev(seq(1, 100, length.out = p)))
    ini[ncoefs] <- log(0.1)
  }
  indB <- cumsum(c(1, rep(k, p)))
  cvec0 <- ini[1:ncoefs]
  cvec0[indB[p + 1]:(ncoefs - 1)] <- ini[(k*p+1):(ncoefs-1)]
  cvec0[ncoefs] <- ini[ncoefs]

  return(cvec0)
}

run_optimization <- function(p_val, k_val, ini, SiCell, newData, est_pts,
                             basis_type, bin_size, optim_control, alpha) {

  # Create basis
  basis <- if (identical(basis_type, "bspline")) {
    fda::create.bspline.basis(c(0, 1), k_val)
  } else {
    fda::create.fourier.basis(c(0, 1), k_val)
  }

  B <- fda::eval.basis(est_pts, basis)
  wghts <- get_weights(newData[, 2], est_pts, alpha)

  # Initialize coefficients
  cvec0 <- if (is.null(ini)) {
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
      weights = wghts,
      method = "BFGS",
      control = list(maxit = 3)
    )
    opt_init$par
  } else if (identical(tolower("LSQ"), tolower(ini))) {
    cvec0 <- init_coeffs(p_val, k_val, NULL)
    opt_init <- optim(
      par = cvec0,
      fn = LSQf,
      gr = LSQg,
      SiCell = SiCell,
      data = newData,
      p = p_val,
      k = k_val,
      B = B,
      estGrid = est_pts,
      weights = wghts,
      method = "BFGS",
      control = list(trace = 0, maxit = 10, fnscale = 1, reltol = 1e-5)
    )
    opt_init$par
  } else {
    init_coeffs(p_val, k_val, ini)
  }

  # Perform main optimization
  n <- length(SiCell)
  ncoefs <- p_val * k_val + p_val + 1

  opt_result <- tryCatch({
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
      weights = wghts,
      method = "BFGS",
      control = optim_control
    )

    # Compute AIC
    aic_value <- n * optim_result$value + (k_val * p_val * p_val) + 1

    list(
      par = optim_result$par,
      value = optim_result$value,
      aic = aic_value,
      convergence = optim_result$convergence
    )
  }, error = function(e) {
    list(
      par = rep(NA, ncoefs),
      value = NA,
      aic = NA,
      convergence = 99
    )
  })

  return(list(
    cvec = opt_result$par,
    aic = opt_result$aic,
    basis = basis,
    convergence = opt_result$convergence,
    LLK = opt_result$value
  ))
}

extract_results <- function(results, indices, AICkp, cvec_cell, basis_cell,
                            convergence_cell, LLK_cell) {
  for (i in indices) {
    if (!is.null(results[[i]])) {
      AICkp[i] <- results[[i]]$aic
      cvec_cell[[i]] <- results[[i]]$cvec
      basis_cell[[i]] <- results[[i]]$basis
      convergence_cell[i] <- results[[i]]$convergence
      LLK_cell[i] <- results[[i]]$LLK
    }
  }

  return(list(
    AICkp = AICkp,
    cvec_cell = cvec_cell,
    basis_cell = basis_cell,
    convergence_cell = convergence_cell,
    LLK_cell = LLK_cell
  ))
}

# Main function
fit_model <- function(p, k, binData, gcvData, basis_type, ini, bin_size,
                      optim_control, all_kp, alpha) {

  # Select appropriate dataset
  data_to_use <- if (length(k) == 1 && length(p) == 1) binData else gcvData
  est_pts <- data_to_use$est_pts
  SiCell <- data_to_use$SiCell
  newData <- data_to_use$newData

  # Initialize storage
  results <- list()

  if (all_kp) {
    # All combinations of p and k
    cat(if (length(p) > 1 && length(k) > 1) "Selecting k & p\n"
        else if (length(p) > 1) "Selecting p\n"
        else if (length(k) > 1) "Selecting k\n")

    grid <- expand.grid(p = p, k = k)
    grid <- grid[grid$k >= grid$p, ]

    n_models <- nrow(grid)
    AICkp <- numeric(n_models)
    cvec_cell <- vector("list", n_models)
    basis_cell <- vector("list", n_models)
    convergence_cell <- numeric(n_models)
    LLK_cell <- numeric(n_models)

    # Grid search
    for (i in seq_len(n_models)) {
      results[[i]] <- run_optimization(
        p_val = grid$p[i],
        k_val = grid$k[i],
        ini = ini,
        SiCell = SiCell,
        newData = newData,
        est_pts = est_pts,
        basis_type = basis_type,
        bin_size = bin_size,
        optim_control = optim_control,
        alpha = alpha
      )
    }

    extracted <- extract_results(results, seq_len(n_models), AICkp, cvec_cell,
                                 basis_cell, convergence_cell, LLK_cell)
    AICkp <- extracted$AICkp
    cvec_cell <- extracted$cvec_cell
    basis_cell <- extracted$basis_cell
    convergence_cell <- extracted$convergence_cell
    LLK_cell <- extracted$LLK_cell

  } else {
    # Stepwise optimization
    if (length(k) > 1) {
      pmax <- max(p)
      cat("Selecting k with p =", pmax, "\n")

      # Step 1: Optimize k with p = pmax
      valid_k <- k[k >= pmax]
      n_k <- length(valid_k)

      AICk <- numeric(n_k)
      cvec_cell_k <- vector("list", n_k)
      basis_cell_k <- vector("list", n_k)
      convergence_cell_k <- numeric(n_k)
      LLK_cell_k <- numeric(n_k)

      for (i in seq_len(n_k)) {
        results[[i]] <- run_optimization(
          p_val = pmax,
          k_val = valid_k[i],
          ini = ini,
          SiCell = SiCell,
          newData = newData,
          est_pts = est_pts,
          basis_type = basis_type,
          bin_size = bin_size,
          optim_control = optim_control,
          alpha = alpha
        )
      }

      extracted_k <- extract_results(results, seq_len(n_k), AICk,
                                     cvec_cell_k,basis_cell_k,
                                     convergence_cell_k, LLK_cell_k)

      # Find best k
      best_k_idx <- which.min(extracted_k$AICkp)
      best_k <- valid_k[best_k_idx]
      AICk_res <- cbind(p = pmax, k = valid_k, AIC = extracted_k$AICkp,
                        converge = extracted_k$convergence_cell,
                        LLK = extracted_k$LLK_cell)
    } else {
      best_k <- k
      AICk_res <- NULL
      extracted_k <- NULL
    }

    if (length(p) > 1) {
      # Step 2: Optimize p with fixed k = best_k
      cat("Selecting p with k =", best_k, "\n")

      valid_p <- p[p <= best_k]
      n_p <- length(valid_p)

      AICp <- numeric(n_p)
      cvec_cell_p <- vector("list", n_p)
      basis_cell_p <- vector("list", n_p)
      convergence_cell_p <- numeric(n_p)
      LLK_cell_p <- numeric(n_p)

      results <- list()  # Reset for p optimization

      for (i in seq_len(n_p)) {
        results[[i]] <- run_optimization(
          p_val = valid_p[i],
          k_val = best_k,
          ini = ini,
          SiCell = SiCell,
          newData = newData,
          est_pts = est_pts,
          basis_type = basis_type,
          bin_size = bin_size,
          optim_control = optim_control,
          alpha = alpha
        )
      }

      extracted_p <- extract_results(results, seq_len(n_p), AICp,
                                     cvec_cell_p, basis_cell_p,
                                     convergence_cell_p, LLK_cell_p)

      # Find best p
      best_p_idx <- which.min(extracted_p$AICkp)
      best_p <- valid_p[best_p_idx]
      AICp_res <- cbind(p = valid_p, k = best_k, AIC = extracted_p$AICkp,
                        converge = extracted_p$convergence_cell,
                        LLK = extracted_p$LLK_cell)
    } else {
      best_p <- p
      AICp_res <- NULL
      extracted_p <- NULL
    }

    # Combine results
    AICkp <- c(if (!is.null(extracted_k)) extracted_k$AICkp else NULL,
               if (!is.null(extracted_p)) extracted_p$AICkp else NULL)
    cvec_cell <- c(if (!is.null(extracted_k)) extracted_k$cvec_cell else NULL,
                   if (!is.null(extracted_p)) extracted_p$cvec_cell else NULL)
    basis_cell <- c(if (!is.null(extracted_k)) extracted_k$basis_cell else NULL,
                    if (!is.null(extracted_p)) extracted_p$basis_cell else NULL)
    convergence_cell <- c(if (!is.null(extracted_k)) extracted_k$convergence_cell else NULL,
                          if (!is.null(extracted_p)) extracted_p$convergence_cell else NULL)
    LLK_cell <- c(if (!is.null(extracted_k)) extracted_k$LLK_cell else NULL,
                  if (!is.null(extracted_p)) extracted_p$LLK_cell else NULL)

    AIC_res <- rbind(AICk_res, AICp_res)
  }

  # Create results table
  AIC_res <- if (all_kp) {
    cbind(p = grid$p, k = grid$k, AIC = AICkp, convergence = convergence_cell,
          LLK = LLK_cell)
  } else {
    AIC_res
  }

  # Find best model
  valid_indices <- which(is.finite(AICkp) & !is.na(AICkp))
  if (length(valid_indices) == 0) stop("No valid AIC values found")

  best_idx <- valid_indices[which.min(AICkp[valid_indices])]

  best_p <- if (all_kp) grid$p[best_idx] else best_p
  best_k <- if (all_kp) grid$k[best_idx] else best_k

  # Refit with best parameters using binData
  final_result <- run_optimization(
    p_val = best_p,
    k_val = best_k,
    ini = ini,
    SiCell = binData$SiCell,
    newData = binData$newData,
    est_pts = binData$est_pts,
    basis_type = basis_type,
    bin_size = bin_size,
    optim_control = list(abstol = 1e-7, maxit = 1000),
    alpha = alpha
  )

  return(list(
    k = best_k,
    p = best_p,
    cvec = final_result$cvec,
    AIC_res = AIC_res,
    basis = final_result$basis,
    convergence = final_result$convergence
  ))
}
