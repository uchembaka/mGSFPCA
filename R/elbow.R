# This function is gotten from Nicolas Casajus Github 
# https://github.com/ahasverus/elbow/blob/master/R/elbow.R


elbow <- function(data){
  
  ## Data transformation ----
  
  data <- data[ , 1:2]
  data <- data[order(data[ , 1]), ]
  
  
  ## Get constant increase/decrease in y ----
  
  constant <- data[c(1, nrow(data)), ]
  colnames(constant) <- c("x", "y")
  
  mod <- stats::lm(y ~ x, data = constant)
  
  data[ , "constant"] <- round(mod$coef[[1]] + mod$coef[[2]] * data[ , 1], 3)
  
  
  ## Detect inflection point ----
  
  pos <- round(nrow(data) / 2)
  
  if (data[pos, "constant"] < data[pos, 2]) { # Concave Down
    
    ymin <- min(data[ , 2])
    data[ , "benefits"] <- ymin + round(data[ , 2] - data[ , "constant"], 3)
    maxi <- data[which.max(data[ , "benefits"]), ]
    
  } else { # Concave Up
    
    ymax <- max(data[ , 2])
    data[ , "benefits"] <- ymax - round(data[ , "constant"] - data[ , 2], 3)
    maxi <- data[which.min(data[ , "benefits"]), ]
  }
  
  
  ## Store results ----
  
  xxx <- list()
  xxx[[1]] <- maxi[1, 1]
  xxx[[2]] <- data
  names(xxx) <- c(paste(colnames(data)[1], "selected", sep = "_"), "data")
  return(xxx)
}
