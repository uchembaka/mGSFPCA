input_check <- function(data, p, k, basis_type, optim_trace,
                        maxit, optim_tol, bin_size, nRegGrid,
                        init_coeff, obs_range, mu_nbasis, use_kp_grid, skip_check) {

  if (skip_check) {
    # Skip all checks and define new terms
    evalGrid <- seq(0, 1, length.out = nRegGrid)
    obs_range <- if (is.null(obs_range)) c(0, 1) else obs_range
    aT <- obs_range[1]
    bT <- obs_range[2]

    # Return inputs and new terms as-is
    return(list(
      data = data,
      p = p,
      k = k,
      basis_type = basis_type,
      optim_trace = optim_trace,
      maxit = maxit,
      optim_tol = optim_tol,
      bin_size = bin_size,
      nRegGrid = nRegGrid,
      init_coeff = init_coeff,
      obs_range = obs_range,
      mu_nbasis = mu_nbasis,
      evalGrid = evalGrid,
      aT = aT,
      bT = bT,
      use_kp_grid = use_kp_grid
    ))
  }

  # Check if data is a matrix with 3 columns and n > 10
  if (ncol(data) != 3) {
    stop("data must be a matrix/dataframe with exactly 3 columns")
  }
  if (nrow(data) <= 10) {
    stop("data must have more than 10 rows")
  }

  # Check if p, k, optim_trace, maxit, bin_size, nRegGrid, mu_nbasis are integers and scalars
  integer_params <- list(p = p, k = k, optim_trace = optim_trace, maxit = maxit,
                         bin_size = bin_size, nRegGrid = nRegGrid, mu_nbasis = mu_nbasis)
  for (param_name in names(integer_params)) {
    param <- integer_params[[param_name]]
    if (length(param) > 1) {
      if (param_name %in% c("p", "k")) {
        # Allow vectors for p and k, validated later
        next
      } else {
        stop(sprintf("%s must be a scalar integer", param_name))
      }
    }
    if (!is.numeric(param) || any(param %% 1 != 0)) {
      stop(sprintf("%s must be an integer", param_name))
    }
  }

  # Check p and k relationship
  p <- sort(p)
  k <- sort(k)
  if (length(p) == 1 && length(k) == 1) {
    if (k < p) {
      stop("When p and k are scalars, k must be >= p")
    }
  } else {
    # At least one of p or k is a vector
    param_grid <- expand.grid(p = p, k = k)
    if (!any(param_grid$k >= param_grid$p)) {
      stop("At least one combination in expand.grid(p, k) must have k >= p")
    }
  }

  # Cap optim_trace at 5
  if (optim_trace > 5) {
    warning("optim_trace capped at 5 (maximum allowed value)")
    optim_trace <- 5
  }

  # Check basis_type
  if (!is.character(basis_type) || !basis_type %in% c("fourier", "bspline")) {
    stop("basis_type must be either 'fourier' or 'bspline'")
  }

  # Check optim_tol is numeric and positive
  if (!is.numeric(optim_tol) || optim_tol <= 0) {
    stop("optim_tol must be a positive numeric value")
  }

  # Check obs_range is a vector of length 2
  if (!is.null(obs_range)) {
    if (!is.vector(obs_range) || length(obs_range) != 2) {
      stop("obs_range must be a vector of length 2")
    }
    if (!is.numeric(obs_range) || any(is.na(obs_range))) {
      stop("obs_range must contain two numeric values")
    }
    if (obs_range[2] <= obs_range[1]) {
      stop("obs_range[2] must be greater than obs_range[1]")
    }
  }

  # Check init_coeff is numeric
  if (!is.null(init_coeff) && !is.numeric(init_coeff)) {
    stop("init_coeff must be numeric")
  }

  # Parameter validation for p, k, and init_coeff
  if (length(p) > 1 || length(k) > 1) {
    if (!is.null(init_coeff)) {
      warning("Either p or k is not a single integer, init_coeff will be ignored.")
      init_coeff <- NULL
    }
  } else if (!is.null(init_coeff)) {
    if (k*p + p + 1 != length(init_coeff)) {
      stop("length of init_coeff not equal to k*p + p + 1")
    }
  }

  # Check logicals
  if (!is.logical(skip_check)) {
    warning("skip_check not logical. Set to TRUE")
  }

  # check validity of use_kp_grid
  if (!is.logical(use_kp_grid)) {
    warning("use_kp_grid not logical. Set to TRUE")
    use_kp_grid <- TRUE
  } else {
    if (length(k) == 1 || length(p) == 1) {
      if (!use_kp_grid) {
        warning("use_kp_grid set to TRUE")
        use_kp_grid <- TRUE
      }
    } else if (any(p > k[1]) && use_kp_grid == FALSE) {
      stop("At least one value in p is greater than min k. Set use_kp_grid to TRUE")
    }
  }


  # Create evaluation grid
  evalGrid <- seq(0, 1, length.out = nRegGrid)
  doRegrid <- chngRange(data[,2], c(0,1))

  if (is.null(obs_range)) {
    if (doRegrid) {
      stop("To do estimation using argvals not defined on [0,1], you must provide the obs_range of your argvals.")
    }
    obs_range <- c(0,1)
  }
  aT <- obs_range[1]
  bT <- obs_range[2]

  # Data preparation
  if (basis_type == "fourier") {
    if (length(k) > 1) {
      k <- k[k %% 2 == 1]
    } else {
      k <- ifelse(k %% 2 == 1, k, k+1)
    }

    message("k: ", paste(k, collapse = ", "))
  }

  # Optim controls
  optim_control <- list(
    trace = optim_trace,
    maxit = maxit,
    fnscale = 1,
    reltol = optim_tol
  )


  # Return validated inputs as a list
  return(list(
    data = data,
    p = p,
    k = k,
    basis_type = basis_type,
    optim_trace = optim_trace,
    maxit = maxit,
    optim_tol = optim_tol,
    optim_control = optim_control,
    bin_size = bin_size,
    nRegGrid = nRegGrid,
    init_coeff = init_coeff,
    obs_range = obs_range,
    mu_nbasis = mu_nbasis,
    evalGrid = evalGrid,
    aT = aT,
    bT = bT,
    use_kp_grid = use_kp_grid
  ))
}
