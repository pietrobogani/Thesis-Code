source("temp1.R")
source("full.R")

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


run_simulation <- function(){
  
  n <- 204
  jtFirstOOS <- 180 #First index for out-of-sample computations
  source("temp1.R")
  source("full.R")
  funs=rq.funs()
  
  predict.fun=funs$predict.fun
  train.fun=funs$train.fun
  
  
  
  #----Generate AR(2)
  phi_ar2 <- c(0.5, -0.2)  # AR coefficients for AR(2)
  Y_ar2 <- numeric(2*n)
  Y_ar2[1] <- 1
  Y_ar2[2] <- -2
  for (i in 3:(2*n)) {
    Y_ar2[i] <- phi_ar2[1] * Y_ar2[i-1] + phi_ar2[2] * Y_ar2[i-2] + rt(n = 1, df = 2) #heavy tails
  }
  Y_ar2 <- Y_ar2[n: (2*n - 1)] #I implemented a burn-in
  
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])
  Y_lag3_ar2 <- c(NA, NA, NA, Y_ar2[1:(length(Y_ar2)-3)])
  
  
  # Prepare the dataset for AR(3) (excluding the first three NA values)
  data_ar2 <- data.frame(I = 1, Y = Y_ar2[-c(1:3)], Y_lag1 = Y_lag1_ar2[-c(1:3)], Y_lag2 = Y_lag2_ar2[-c(1:3)], Y_lag3 = Y_lag3_ar2[-c(1:3)])
  
  
  
  #QUANTILE REGRESSION (CORRETTA, UGUALE AL PAPER ORIGINALE)
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ, 0.99)
  
  QuantAR1_OOS <- matrix(NA, nrow(data_ar2), length(QQ))
  QuantAR2_OOS <- matrix(NA, nrow(data_ar2), length(QQ))
  QuantAR3_OOS <- matrix(NA, nrow(data_ar2), length(QQ))
  
  
  # QUANTILE REGRESSION OUT OF SAMPLE
  for (jt in jtFirstOOS:nrow(data_ar2)) {
    for (jq in 1:length(QQ)) {  
      
      QR1 <- rq(Y[1:(jt-1)] ~ Y_lag1[1:(jt - 1)], data = data_ar2, tau=QQ[jq])
      QuantAR1_OOS[jt, jq] <- as.matrix(data_ar2[jt, - c(2,4,5) ]) %*% coef(QR1)
      
      # QR2 <- rq(Y[1:(jt-1)] ~ Y_lag1[1:(jt - 1)] + Y_lag2[1:(jt - 1)], data = data_ar2, tau=QQ[jq])
      # QuantAR2_OOS[jt, jq] <- as.matrix(data_ar2[jt, - c(2,5) ]) %*% coef(QR2)
      # 
      # QR3 <- rq(Y[1:(jt-1)] ~ Y_lag1[1:(jt - 1)] + Y_lag2[1:(jt - 1)] + Y_lag3[1:(jt - 1)], data = data_ar2, tau=QQ[jq])
      # QuantAR3_OOS[jt, jq] <- as.matrix(data_ar2[jt, - 2 ]) %*% coef(QR3)
      
    }
  }
  
  #-------------------------------------------------------------------------------------------
  
  
  # Distributional conformal Prediction OUT OF SAMPLE
  DQuantAR1_OOS <- matrix(NA, nrow(data_ar2), length(QQ))
  
  
  for (jt in jtFirstOOS:nrow(data_ar2)) {
    
    Y <- data_ar2$Y[1:(jt-1)]
    X <- data_ar2$Y_lag1[1:(jt-1)]
    Xnew <- data_ar2$Y_lag1[jt]
    ytrial <- seq(-max(abs(Y)), max(abs(Y)), 0.01)
    
    for (jq in 1:length(QQ)) {  
      test = dist.conformal.pred(X,Y,Xnew,train.fun = train.fun, predict.fun = predict.fun, verbose = T, alpha = 1 - QQ[jq])
      DQuantAR1_OOS[jt,jq] <- test$up[1,1]
    }
  }
  
  
  #-------------------------------------------------------------------------------------------
  
  #CALIBRATION OF QUANTILE REGRESSION
  
  # This function returns the cumulative probability for a given value X, assuming a step function as CDF
  
  cumulative_prob <- function(X, quantiles, values) {
    if (length(values[values <= X]) == 0) {
      return(0) # If X is less than the smallest value, the cumulative probability is 0
    }
    
    max_value_not_exceeding_X <- max(values[values <= X], na.rm = TRUE)
    
    quantile_index <- which(values == max_value_not_exceeding_X)
    if (length(quantile_index) == 0) {
      return(1) # If X is greater than all values, return the maximum cumulative probability
    }
    
    return(quantiles[min(quantile_index)])
  }
  
  
  PitSTAR1_OOS <- rep(NA, nrow(data_ar2))

  DPitSTAR1_OOS <- rep(NA, nrow(data_ar2))

  
  
  
  for (jt in jtFirstOOS: nrow(data_ar2) ){
    
    YhRealized <- data_ar2[jt,2]
    
    PitSTAR1_OOS[jt] <- cumulative_prob(YhRealized, QQ, QuantAR1_OOS[jt, ])

    DPitSTAR1_OOS[jt] <- cumulative_prob(YhRealized, QQ, DQuantAR1_OOS[jt, ])

    
  }
  ret <- list(PitSTAR1_OOS, DPitSTAR1_OOS)
  
  return(ret)
}




n_simul <- 10
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

# Setup parallel cluster
cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
clusterExport(cl, varlist=c("run_simulation")) # Export the simulation function to each cluster node
clusterEvalQ(cl, { 
  library(readxl)
  library(quantreg)}
) # Load required libraries in each cluster node, repeat as necessary for other libraries

# Run simulations in parallel
results <- parLapply(cl, seeds, function(seed) {
  #set.seed(seed)
  run_simulation()
})


# Stop the cluster
stopCluster(cl)

# Extract results
PitSTAR1_OOStot <- numeric(0)
PitSTAR2_OOStot <- numeric(0)
PitSTAR3_OOStot <- numeric(0)
DPitSTAR1_OOStot <- numeric(0)
DPitSTAR2_OOStot <- numeric(0)
DPitSTAR3_OOStot <- numeric(0)




for(res in results){
  PitSTAR1_OOStot <- c(PitSTAR1_OOStot, na.omit(res[[1]]))
  PitSTAR2_OOStot <- c(PitSTAR2_OOStot,  na.omit(res[[2]]))
  PitSTAR3_OOStot <- c(PitSTAR3_OOStot, na.omit(res[[3]]))
  DPitSTAR1_OOStot <- c(DPitSTAR1_OOStot, na.omit(res[[4]]))
  DPitSTAR2_OOStot <- c(DPitSTAR2_OOStot, na.omit(res[[5]]))
  DPitSTAR3_OOStot <- c(DPitSTAR3_OOStot, na.omit(res[[6]]))
  
}






#------------------- WILSON SCORES


wilson_score_interval <- function(x, n, conf.level = 0.95) { #fake wilson score, è un altro
  # x: number of successes
  # n: number of trials
  # conf.level: confidence level (e.g., 0.95 for 95% confidence interval)
  
  # Calculate point estimate for proportion
  p_hat <- x / n
  
  # Calculate standard error using the normal approximation
  se <- sqrt(p_hat * (1 - p_hat) / n)
  
  # Find the Z value for the specified confidence level
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)
  
  # Calculate confidence interval
  lower_bound <- p_hat - z * se
  upper_bound <- p_hat + z * se
  
  return(c(lower = max(0, lower_bound), upper = min(1, upper_bound)))
}


is_within_ci <- function(success_rate, ci_lower, ci_upper) {
  return(success_rate >= ci_lower & success_rate <= ci_upper)
}





resultsPitST <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

PitSTAR1_OOStot <- na.omit(PitSTAR1_OOStot)
for (quantile in c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)) {
  successes <- sum(PitSTAR1_OOStot < quantile)
  n <- length(PitSTAR1_OOStot)
  
  # Calculate success rate
  success_rate <- successes / n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsPitST <- rbind(resultsPitST, data.frame(Quantile = quantile, SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}



# Apply the function to each row and add the result as a new column
resultsPitST$IsWithinCI <- mapply(is_within_ci, resultsPitST$Quantile, resultsPitST$CI_Lower, resultsPitST$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsPitST$IsWithinCI) * 100

# Print the result
print(paste("Percentage of quantiles within the confidence interval:", percentage_within_ci))









resultsDPitST <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

DPitSTAR1_OOStot <- na.omit(DPitSTAR1_OOStot)
n <- length(DPitSTAR1_OOStot)

for (quantile in c(0.01, seq(0.05, 0.95, by = 0.025), 0.99)) {
  successes <- sum(DPitSTAR1_OOStot < quantile)
  
  # Calculate success rate
  success_rate <- successes / n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsDPitST <- rbind(resultsDPitST, data.frame(Quantile = quantile, SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsDPitST$IsWithinCI <- mapply(is_within_ci, resultsDPitST$Quantile, resultsDPitST$CI_Lower, resultsDPitST$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsDPitST$IsWithinCI) * 100

# Print the result
print(paste("Percentage of quantiles within the confidence interval:", percentage_within_ci))


