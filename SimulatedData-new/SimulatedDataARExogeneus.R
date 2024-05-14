#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#   RESULT SUMMARY

           


library(caret)
library(sn)
library(quantreg)
library(readxl)
library(forecast)
library(SuppDists)
library(readr)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)
library(progressr)
library(future)
library(future.apply)
library(progressr)
library(mvtnorm)
library(quantregForest)



# da qui comincia la simulazione

run_simulation <- function(){
  
  # Set parameters
  n <- 1000       # Number of data points for training
  n2 <- 100      # Number of data points for testing
  p <- 0.5 * n
  phi_ar2 <- c(0.5, -0.2)   # AR(2) coefficients
  num_exog_vars <-  p     # Parameter to control the number of exogenous variables
  beta <- runif(num_exog_vars, 0.1, 1) # Randomly generated coefficients for exogenous variables
  means <- rnorm(p)
  sigma <- diag(runif(p))
  
  # Initialize variables
  Y_ar2 <- numeric(n)  # Series for AR(2) data
  exog_vars <- matrix(nrow = n + 1, ncol = num_exog_vars)  # Matrix for exogenous variables
  
  exog_vars <- rmvnorm( n = n + 1, mean = means, sigma = sigma)
# 
#    for (i in 1:num_exog_vars) {
#      segment_length <- num_exog_vars / 8  # Using 8 diverse segment types
#      
#      if (i <= segment_length) {
#        # Linear trends with a cap
#        exog_vars[, i] <- pmin(seq(from = 1, to = n + 1, length.out = n + 1) * runif(1, 0.5, 1.5), 1000)
#        
#      } else if (i <= segment_length * 2) {
#        # Controlled exponential decay avoiding Inf at the end
#        decay_rate <- runif(1, 0.005, 0.02)
#        exog_vars[, i] <- pmin(exp(-decay_rate * (1:(n + 1))), 1000)
#        
#      } else if (i <= segment_length * 3) {
#        # Logarithmic trends, capped
#        shift <- abs(rnorm(1, 10 * i, 5)) + 0.01
#        scale <- rnorm(1, i, 0.5)
#        exog_vars[, i] <- pmin(log(shift + scale * (1:(n + 1))), 1000)
#        
#      } else if (i <= segment_length * 4) {
#        # Sinusoidal with increasing amplitude, capped
#        frequency <- runif(1, 15, 30)
#        phase <- runif(1, 0, 2 * pi)
#        amplitude_growth <- runif(1, 0.0001, 0.0005)
#        exog_vars[, i] <- pmin((sin(2 * pi * (1:(n + 1)) / frequency + phase) + 1) * (1 + amplitude_growth * (1:(n + 1))), 1000)
#        
#      } else {
#        # Polynomial, AR, or other types, use similar strategies to avoid extremes
#        # This covers the remaining cases, ensuring all exogenous variables are generated
#        exog_vars[, i] <- pmin((1:(n + 1))^runif(1, -2, -0.5), 1000)  # Polynomial decreasing
#      }
#    }

  
  # Initialize Y values
  Y_ar2[1:2] <- rnorm(2)  # Initial values for the AR(2) process
  
  # Simulate AR(2) with exogenous variables for training data
  for (t in 3:n) {
    Y_ar2[t] <- phi_ar2[1] * Y_ar2[t - 1] + phi_ar2[2] * Y_ar2[t - 2] +
      sum(beta * exog_vars[t, ]) + rt(n = 1, df = 2)
  }
  
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])
  
  # Prepare dataset for model excluding the first two NA values
  data_ar2 <- data.frame(I = 1, Y = Y_ar2[-c(1:2)], Y_lag1 = Y_lag1_ar2[-c(1:2)], Y_lag2 = Y_lag2_ar2[-c(1:2)])
  for (j in 1:num_exog_vars) {
    data_ar2[paste0("X", j)] <- exog_vars[-c(1:2, n+1), j]
  }
  
  # Extract exogenous variables at time n+1 for testing
  exog_vars_n_plus_1 <- exog_vars[n+1, ]
  
  # Simulate test data using the same exogenous variables at n+1
  Y_test <- phi_ar2[1] * Y_ar2[n] + phi_ar2[2] * Y_ar2[n -1] +
    sum(beta * exog_vars_n_plus_1) + rt(n = n2, df = 2)
  
  
  make_formula <- function(num_lags, num_exog) {
    lag_terms <- paste(c("Y_lag1", if (num_lags >= 2) "Y_lag2"), collapse = " + ")
    exog_terms <- paste(paste0("X", 1:num_exog), collapse = " + ")
    formula <- as.formula(paste("Y ~", lag_terms, "+", exog_terms))
    return(formula)
  }
  
#--------------- QUANTILE REGRESSION 
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)

  QuantAR1_OOS <- matrix(0, length(QQ))
  QuantAR2_OOS <- matrix(NA, length(QQ))
  #qrf_model2 <- quantregForest(x = as.matrix(data_ar2[,-c(1,2)]), y = as.matrix(data_ar2[,2]))
  

    for (jq in 1:length(QQ)) {  
       
      # QR1 <- rq(make_formula(1, num_exog_vars), data = data_ar2, tau=QQ[jq])     
      # QuantAR1_OOS[jq] <- c(1, Y_ar2[n], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR1)
      
      
      
      # qrf_model1 <- quantregForest(x = as.matrix(data_ar2[,-c(1,2,4)]), y = as.matrix(data_ar2[,2]))
      # QuantAR1_OOS[jq] <- predict(qrf_model1, newdata = matrix(rep(c(Y_ar2[n],exog_vars_n_plus_1),each = n2), nrow = n2), what = QQ[jq])[1]
      # 
      
      QR2 <- rq(make_formula(2, num_exog_vars), data = data_ar2, tau=QQ[jq])
      QuantAR2_OOS[jq] <- c(1, Y_ar2[n], Y_ar2[n-1], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR2)
      
      #QuantAR2_OOS[jq] <- predict(qrf_model2, newdata = matrix(rep(c(Y_ar2[n],Y_ar2[n-1],exog_vars_n_plus_1),each = n2), nrow = n2), what = QQ[jq])[1]
       
      
    }
  
  

  
  
#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  CQuantAR1_OOS <- matrix(0, length(QQ))
  CQuantAR2_OOS <- matrix(NA, length(QQ))
  Q_low <- matrix(NA, nrow(data_ar2), length(QQ))
  Q_high <- matrix(NA, nrow(data_ar2), length(QQ))
  full_length <- nrow(data_ar2)
  test_length = floor(full_length*66/100)
  train_indices <- sample(1:test_length)
  # data_ar2_1 <- data_ar2[1:test_length,] #train
  # data_ar2_2 <- data_ar2[(test_length + 1) : full_length,] #calibration
  data_ar2_1 <- data_ar2[train_indices,] #train
  data_ar2_2 <- data_ar2[-train_indices,] #calibration
  
  #qrf_model2 <- quantregForest(x = as.matrix(data_ar2_1[,-c(1,2)]), y = as.matrix(data_ar2_1[,2]))
  

    for (jq in 1:length(QQ)) {  
      
      #AR(1)
      # Q_low[1:nrow(data_ar2_2), jq] <- -Inf 
      # # QR <- rq(make_formula(1, num_exog_vars), data = data_ar2_1, tau=QQ[jq])    
      # # Q_high[1 : nrow(data_ar2_2), jq] <- as.matrix(data_ar2_2[,-c(2,4)]) %*% coef(QR) 
      # 
      # qrf_model1 <- quantregForest(x = as.matrix(data_ar2_1[,-c(1,2,4)]), y = as.matrix(data_ar2_1[,2]))
      # Q_high[1 : nrow(data_ar2_2), jq] <- predict(qrf_model1, newdata = as.matrix(data_ar2_2[,-c(1,2,4)]), what = QQ[jq])
      # 
      # # Initialize a vector for errors
      # E_i <- rep(NA, nrow(data_ar2_2))
      # 
      # # Calculate errors for each point in the test set I2
      # for (i in 1:length(E_i)) {
      #   E_i[i] <- max(Q_low[i, jq] - data_ar2_2[i,2], data_ar2_2[i,2] - Q_high[i, jq])
      # }
      # 
      # 
      # # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
      # quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(data_ar2_2)))
      # 
      # #CQuantAR1_OOS[jq] <- c(1, Y_ar2[n], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR) + quantile_E
      # CQuantAR1_OOS[jq] <- predict(qrf_model1, newdata = matrix(rep(c(Y_ar2[n],exog_vars_n_plus_1),each = n2), nrow = n2), what = QQ[jq])[1] + quantile_E
      
       
      #AR(2)
       Q_low[1:nrow(data_ar2_2), jq] <- -Inf 
       QR <- rq(make_formula(2, num_exog_vars), data = data_ar2_1, tau= QQ[jq])
       Q_high[1 : nrow(data_ar2_2), jq] <- as.matrix(data_ar2_2[,-c(2)]) %*% coef(QR) 
       # 
       #Q_high[1 : nrow(data_ar2_2), jq] <- predict(qrf_model2, newdata = as.matrix(data_ar2_2[,-c(1,2)]), what = QQ[jq])

       # Initialize a vector for errors
       E_i <- rep(NA, nrow(data_ar2_2))
       
       # Calculate errors for each point in the test set I2
       for (i in 1:length(E_i)) {
         E_i[i] <- max(Q_low[i, jq] - data_ar2_2[i,2], data_ar2_2[i,2] - Q_high[i, jq])
       }
       
       # Compute Q(QQ[jq])(E, I2) N.B 1 - ?? = QQ[jq]
       quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(data_ar2_2)))
       
       CQuantAR2_OOS[jq] <- c(1, Y_ar2[n], Y_ar2[n-1], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR) + quantile_E
       #CQuantAR2_OOS[jq] <- predict(qrf_model2, newdata = matrix(rep(c(Y_ar2[n],Y_ar2[n-1],exog_vars_n_plus_1),each = n2), nrow = n2), what = QQ[jq])[1] + quantile_E
       
     
    }
  
  
  
  #-------------------------------------------------------------------------------------------
  
  #CALIBRATION OF QR AND CQR
  
  # This function returns the cumulative probability for a given value X, assuming a step function as CDF
  coverageQRAR1 <- rep(NA,length(QQ))
  coverageQRAR2 <- rep(NA,length(QQ))
  coverageCQRAR1 <- rep(NA,length(QQ))
  coverageCQRAR2 <- rep(NA,length(QQ))
  
  
  
  for (jq in 1:length(QQ)) {
  PitQROOSAR1 <- rep(NA, length(Y_test))
  PitQROOSAR2 <- rep(NA, length(Y_test))
  
  PitCQROOSAR1 <- rep(NA, length(Y_test))
  PitCQROOSAR2 <- rep(NA, length(Y_test))
  
  for (i in 1:length(Y_test)){
    if (Y_test[i] <= QuantAR1_OOS[jq]) {
      PitQROOSAR1[i]  <- 1 
    }
    else 
      PitQROOSAR1[i]  <- 0
    
    if (Y_test[i] <= QuantAR2_OOS[jq]) {
      PitQROOSAR2[i]  <- 1 
    }
    else 
      PitQROOSAR2[i]  <- 0
  }
  
  coverageQRAR1[jq] <- sum(PitQROOSAR1) / length(PitQROOSAR1)
  coverageQRAR2[jq] <- sum(PitQROOSAR2) / length(PitQROOSAR2)
  

  for (i in 1:length(Y_test)){
    if (Y_test[i] <= CQuantAR1_OOS[jq]) {
      PitCQROOSAR1[i]  <- 1 
    }
    else 
      PitCQROOSAR1[i]  <- 0
    
    if (Y_test[i] <= CQuantAR2_OOS[jq]) {
      PitCQROOSAR2[i]  <- 1 
    }
    else 
      PitCQROOSAR2[i]  <- 0
  }
  
  coverageCQRAR1[jq] <- sum(PitCQROOSAR1) / length(PitCQROOSAR1)
  coverageCQRAR2[jq] <- sum(PitCQROOSAR2) / length(PitCQROOSAR2)
  
  }
  
    
  ret <- list(coverageQRAR1, coverageQRAR2, coverageCQRAR1, coverageCQRAR2)
  
  return(ret)
}




n2 = 100
n3 = 100 
n_simul <- n3
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(0.01, QQ)

# Setup parallel cluster
cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
clusterExport(cl, varlist=c("run_simulation")) # Export the simulation function to each cluster node
clusterEvalQ(cl, { 
  library(readxl)
  library(quantreg)
  library(quantregForest)
  library(mvtnorm)
}
) # Load required libraries in each cluster node, repeat as necessary for other libraries

# Run simulations in parallel
results <- parLapply(cl, seeds, function(seed) {
  set.seed(seed)
  run_simulation()
})


# Stop the cluster
stopCluster(cl)

# Extract results
coveragetotQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotQRAR2 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR2 <- matrix(NA,n3,length(QQ))



index = 1
for(res in results){
  coveragetotQRAR1[index,] <- res[[1]]
  coveragetotQRAR2[index,] <- res[[2]]
  coveragetotCQRAR1[index,] <- res[[3]]
  coveragetotCQRAR2[index,] <- res[[4]]
  index <- index + 1
}



# RESULTS EVALUATION

# Calculate average coverage across all iterations and for each quantile
average_coverageQRAR1 <- rep(NA,length(QQ))
average_coverageQRAR2 <- rep(NA,length(QQ))
average_coverageCQRAR1 <- rep(NA,length(QQ))
average_coverageCQRAR2 <- rep(NA,length(QQ))

for (jq in 1:length(QQ)) { 
  average_coverageQRAR1[jq] <- mean(coveragetotQRAR1[,jq])
  average_coverageQRAR2[jq] <- mean(coveragetotQRAR2[,jq])
  
  average_coverageCQRAR1[jq] <- mean(coveragetotCQRAR1[,jq])
  average_coverageCQRAR2[jq] <- mean(coveragetotCQRAR2[,jq])
}


wilson_score_interval <- function(x, n, conf.level = 0.95) {
  # x: number of successes
  # n: number of trials
  # conf.level: confidence level (e.g., 0.95 for 95% confidence interval)
  
  # Calculate point estimate for proportion
  p_hat <- x / n
  
  # Find the Z value for the specified confidence level
  z <- qnorm(1 - (1 - conf.level) / 2)
  
  # Wilson score interval formula
  factor <- z^2 / (2 * n)
  denominator <- 1 + z^2 / n
  center_adjustment <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
  
  lower_bound <- (p_hat + factor - center_adjustment) / denominator
  upper_bound <- (p_hat + factor + center_adjustment) / denominator
  
  # Adjust bounds to ensure they remain within [0, 1]
  lower_bound <- max(0, lower_bound)
  upper_bound <- min(1, upper_bound)
  
  return(c(lower = lower_bound, upper = upper_bound))
}

is_within_ci <- function(success_rate, ci_lower, ci_upper) {
  return(success_rate >= ci_lower & success_rate <= ci_upper)
}



#----------------------- QR AR(1)
{
# resultsPitSTQRAR1 <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
# 
# for (jq in 1:length(QQ)) {
#   
#   success_rate <- average_coverageQRAR1[jq]
#   
#   n <- n2*n3
#   
#   successes <- success_rate*n
#   
#   # Calculate Wilson score interval
#   ci <- wilson_score_interval(successes, n)
#   
#   # Add to results data frame
#   resultsPitSTQRAR1 <- rbind(resultsPitSTQRAR1, data.frame(Quantile = QQ[jq], SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
# }
# 
# 
# # Apply the function to each row and add the result as a new column
# resultsPitSTQRAR1$IsWithinCI <- mapply(is_within_ci, resultsPitSTQRAR1$Quantile, resultsPitSTQRAR1$CI_Lower, resultsPitSTQRAR1$CI_Upper)
# 
# # Calculate the percentage of quantiles within the CI
# percentage_within_ci <- mean(resultsPitSTQRAR1$IsWithinCI) * 100
# 
# # Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
# count_above <- sum(resultsPitSTQRAR1$SuccessRate > resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI)
# 
# # Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
# count_below <- sum(resultsPitSTQRAR1$SuccessRate < resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI)
# 
# # Print the result
# print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
# print(paste("Number of times SuccessRate is above the Quantile:", count_above))
# print(paste("Number of times SuccessRate is below the Quantile:", count_below))
# print(paste("MAE of QR:", mean(abs(resultsPitSTQRAR1$Quantile-resultsPitSTQRAR1$SuccessRate))))
}


#----------------------- QR AR(2)
{
resultsPitSTQRAR2 <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

for (jq in 1:length(QQ)) {
  
  success_rate <- average_coverageQRAR2[jq]
  
  n <- n2*n3
  
  successes <- success_rate*n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsPitSTQRAR2 <- rbind(resultsPitSTQRAR2, data.frame(Quantile = QQ[jq], SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsPitSTQRAR2$IsWithinCI <- mapply(is_within_ci, resultsPitSTQRAR2$Quantile, resultsPitSTQRAR2$CI_Lower, resultsPitSTQRAR2$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsPitSTQRAR2$IsWithinCI) * 100

# Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsPitSTQRAR2$SuccessRate > resultsPitSTQRAR2$Quantile & !resultsPitSTQRAR2$IsWithinCI)

# Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsPitSTQRAR2$SuccessRate < resultsPitSTQRAR2$Quantile & !resultsPitSTQRAR2$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))
print(paste("MAE of QR:", mean(abs(resultsPitSTQRAR2$Quantile-resultsPitSTQRAR2$SuccessRate))))
}


#----------------------- CQR AR(1)
{
# resultsPitSTCQRAR1 <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
# 
# for (jq in 1:length(QQ)) {
#   
#   success_rate <- average_coverageCQRAR1[jq]
#   
#   n <- n2*n3
#   
#   successes <- success_rate*n
#   
#   # Calculate Wilson score interval
#   ci <- wilson_score_interval(successes, n)
#   
#   # Add to results data frame
#   resultsPitSTCQRAR1 <- rbind(resultsPitSTCQRAR1, data.frame(Quantile = QQ[jq], SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
# }
# 
# 
# # Apply the function to each row and add the result as a new column
# resultsPitSTCQRAR1$IsWithinCI <- mapply(is_within_ci, resultsPitSTCQRAR1$Quantile, resultsPitSTCQRAR1$CI_Lower, resultsPitSTCQRAR1$CI_Upper)
# 
# # Calculate the percentage of quantiles within the CI
# percentage_within_ci <- mean(resultsPitSTCQRAR1$IsWithinCI) * 100
# 
# # Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
# count_above <- sum(resultsPitSTCQRAR1$SuccessRate > resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI)
# 
# # Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
# count_below <- sum(resultsPitSTCQRAR1$SuccessRate < resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI)
# 
# # Print the result
# print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
# print(paste("Number of times SuccessRate is above the Quantile:", count_above))
# print(paste("Number of times SuccessRate is below the Quantile:", count_below))
# print(paste("MAE of CQR:", mean(abs(resultsPitSTCQRAR1$Quantile-resultsPitSTCQRAR1$SuccessRate))))
 }


#----------------------- CQR AR(2)
{
resultsPitSTCQRAR2 <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

for (jq in 1:length(QQ)) {
  
  success_rate <- average_coverageCQRAR2[jq]
  
  n <- n2*n3
  
  successes <- success_rate*n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsPitSTCQRAR2 <- rbind(resultsPitSTCQRAR2, data.frame(Quantile = QQ[jq], SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsPitSTCQRAR2$IsWithinCI <- mapply(is_within_ci, resultsPitSTCQRAR2$Quantile, resultsPitSTCQRAR2$CI_Lower, resultsPitSTCQRAR2$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsPitSTCQRAR2$IsWithinCI) * 100

# Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsPitSTCQRAR2$SuccessRate > resultsPitSTCQRAR2$Quantile & !resultsPitSTCQRAR2$IsWithinCI)

# Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsPitSTCQRAR2$SuccessRate < resultsPitSTCQRAR2$Quantile & !resultsPitSTCQRAR2$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))
print(paste("MAE of CQR:", mean(abs(resultsPitSTCQRAR2$Quantile-resultsPitSTCQRAR2$SuccessRate))))
}



















# filename <- paste("Sim500_1", ".RData",sep="")
# cat(paste("Saving results to file", filename, "\n"))
# 
# # Save all the variables to the .RData file
# save(
#   PitST_OOStot,
#   CPitST_OOStot,
#   PitSTGDPonly_OOStot,
#   CPitSTGDPonly_OOStot,
#   file=filename
# )
# 


# #-------------------------------- per caricare dati salvati -----------------------
# 
# filename <- paste("Sim500_1", ".RData", sep="")
# cat(paste("Loading results from file", filename, "\n"))
# 
# load(filename)
# 
# #-----------------------------------------------------------------------------------





