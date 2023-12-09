split_conformal_quantile_regression <- function(data, alpha,indices_split) {
  
  #devo decidere cosa dare in output, secondo me ha senso solo che riceva anche tanti punti di cui fare prediction e restituisca molti intervalli
  
  #Until 4th quarter of 1992 is training set, test afterwards, last 10 validation
    #data <-X
  NFCI <- rnorm(10000,0,1)
  A191RL1Q225SBEA <- 2* NFCI + rnorm(10000,0,1)
  data <- data.frame(NFCI, A191RL1Q225SBEA)
  indices_split=81
  alpha <- 0.3
  
   
  indices <- sample(1:nrow(data))
  
  # I1 <- indices[1:(indices_split - 1)]
  # I2 <- indices[indices_split:(length(indices)-10)]
  # I3 <- indices[(length(indices)-9) : length(indices)]
  
  # I1 <- indices[1:58]
  # I2 <- indices[59:115]
  # I3 <- indices[116:172]
  
  I1 <- indices[1:333]
  I2 <- indices[334:667]
  I3 <- indices[668:1000]
  
  train <- data[I1,]
  test <- data[I2,]
  validation <- data[I3,]
  
  # Fit quantile regression on the first split (I1)
  
  fit_lo <- rq(train$A191RL1Q225SBEA ~ train$NFCI, tau = 0.3) 
  fit_hi <- rq(train$A191RL1Q225SBEA ~ train$NFCI, tau = alpha+0.3) 
  
  # Extract quantile functions evaluated on test
  X_test <- model.matrix(~ NFCI, data = test) #first column for the intercept
  q_alpha_lo <- X_test %*% coef(fit_lo) 
  q_alpha_hi <- X_test %*% coef(fit_hi) 
  
  #Let's check if we're doing correctly:
  plot(test$NFCI, test$A191RL1Q225SBEA)
  
  # Add the quantile lines to the plot
  lines(test$NFCI, q_alpha_lo, col = "blue", lwd = 2)
  lines(test$NFCI, q_alpha_hi, col = "red", lwd = 2)
  
  
  

  # Compute Ei for each i in the second split (I2)
  E_i <- rep(NA, length(I2))
  for (i in seq_along(I2)) {
    E_i[i] <- max(q_alpha_lo[i] - test$A191RL1Q225SBEA[i], test$A191RL1Q225SBEA[i] - q_alpha_hi[i]) 
  }
  
  # Compute Q(1-??)(E, I2)
  quantile_E <- quantile(E_i, pmin(1, pmax(0, (1 - alpha) * (1 + 1/length(I2))))) #otherwise going out of [0,1]
  
  # Output: Prediction interval C(x) for the whole range of values
  points <- seq(min(data$NFCI), max(data$NFCI), by = 0.001)
  X_eval <- model.matrix(~ NFCI, data = data.frame(NFCI = points))
  
  q_alpha_lo_tot <- X_eval %*% coef(fit_lo) #I extract the predicted 0.01-quantile for new years after 1993 
  q_alpha_hi_tot <- X_eval %*% coef(fit_hi)
  prediction_interval_lower <- q_alpha_lo_tot - quantile_E 
  prediction_interval_upper <- q_alpha_hi_tot + quantile_E
  
  
  
   # Plot the points connected by a line
   plot(data$NFCI, data$A191RL1Q225SBEA, main = "Quantile Regression", xlab = "NFCI", ylab = "A191RL1Q225SBEA")
   points(points, prediction_interval_upper[order(prediction_interval_upper)], type = "l", col = "green", lwd = 2, main = "Connected Line Plot", xlab = "X", ylab = "Y")
   
   points(points, prediction_interval_lower[order(prediction_interval_lower)], type = "l", col = "blue", lwd = 2, main = "Connected Line Plot", xlab = "X", ylab = "Y")
   
   
  
  
  interpolated_low_quantile_f <- approx(points, prediction_interval_lower, xout = validation$NFCI)$y
  #na_indices <- which(is.na(interpolated_low_quantile_f))
  # Calculate the percentage of points in the dataset that fall below the interpolated curve
  percentage_below <- mean(validation$A191RL1Q225SBEA < interpolated_low_quantile_f) * 100
  
  cat("Percentage of points below the  below curve:", percentage_below, "%\n")
  
  
  interpolated_high_quantile_f <- approx(points, prediction_interval_upper, xout = validation$NFCI)$y
  na_indices <- which(is.na(interpolated_high_quantile_f))
  
  # Calculate the percentage of points in the dataset that fall below the interpolated curve
  percentage_below_upperquantile <- mean(validation$A191RL1Q225SBEA[-na_indices] < interpolated_high_quantile_f[-na_indices]) * 100
  
  cat("Percentage of points below the curve:", percentage_below_upperquantile, "%\n")
  
  
  
  return(percentage_below_upperquantile)
}

