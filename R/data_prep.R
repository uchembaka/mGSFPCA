data_prep <- function(data, k, p, obs_range, aT, bT, bin_size, nRegGrid, basis_type, init_coeff, mu_nbasis) {

  data <- as.matrix(data)
  # Round time values to 4 decimal places
  data[,2] <- round(data[,2], 4)

  # Order data by ID and time
  IDs <- unique(data[,1])
  # n <- length(IDs)
  ord_data <- data[order(data[,1], data[,2]),]
  ord_data[,2] <- (ord_data[,2] - aT)/(bT - aT)

  # Get unique observation points
  new_obs_pts <- sort(unique(ord_data[,2]))

  # Create estimation grid
  # if (bin_size > length(new_obs_pts)) {
  #   bin_size <- length(new_obs_pts)
  #   estGrid <- new_obs_pts
  # } else {
    estGrid <- seq(0, 1, length.out = bin_size)
  # }

  if (bin_size <= 51) {
    gcvGrid <- estGrid
  } else {
    gcvGrid <- seq(0, 1, length.out = 51) # this is fixed
  }

  # Estimate mean function
  mean_result <- get_mean(ord_data, estGrid, basis_type, mu_nbasis)

  gcvData <- get_est_data(ord_data, gcvGrid, mean_result$fit)
  binData <- get_est_data(ord_data, estGrid, mean_result$fit)

  # Initialize parameters
  kp1 <- length(k) == 1 & length(p) == 1
  if (kp1 && !is.null(init_coeff)){
    ini <- init_coeff
  }else{
    ini <- NULL
  }


  # Return required components
  return(list(
    gcvData = gcvData,
    binData = binData,
    bin_size = bin_size,
    ini = ini,
    mu_results = mean_result,
    kp1 = kp1
  ))
}


get_est_data <- function (ord_data, est_grid, mu_fit) {
  IDs <- unique(ord_data[,1])
  # Get unique observation points
  new_obs_pts <- sort(unique(ord_data[,2]))

  # Create gcv grid
  if (length(new_obs_pts) < length(est_grid)) {
    est_grid <- new_obs_pts
  }

  # Prepare centered data
  cY <- ord_data[,3] - mu_fit
  data <- cbind(ord_data, cY)

  # Bin data to est_grid
  est_data <- data
  if (length(new_obs_pts) != length(est_grid) && !all(new_obs_pts %in% est_grid)) {
    est_pts <- sort(unique(binData(new_obs_pts, est_grid)))
    est_data[,2] <- binData(est_data[,2], est_grid)
  } else {
    est_pts <- new_obs_pts
  }

  # Add column for time indices
  est_data <- cbind(est_data, match(est_data[,2], est_pts))

  # Bin data by index
  newData <- do.call(rbind, lapply(IDs, function(id) {
    datai <- est_data[est_data[,1] == id, , drop = FALSE]
    bin_result <- averageByIndex(datai[,c(3,4), drop = FALSE], datai[,5])
    datai <- datai[bin_result$firstOccurrence, , drop = FALSE]
    datai[,c(3,4)] <- bin_result$Y
    datai
  }))

  # Get covariance components
  SiCell <- get_SiCell(newData)

  return(list(SiCell = SiCell,
              newData = newData,
              est_pts = est_pts))
}
