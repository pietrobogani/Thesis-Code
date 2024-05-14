#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#   RESULT SUMMARY

#   n1 = 98:
# QR AR(1): 9 - 2 - 11 , MAE : 0.009415
# CQR AR(1): 12 - 5 - 3 , MAE : 0.007505

# QR AR(2): 7 - 7 - 6 , MAE : 0.009155
# CQR AR(2): 8 - 12 - 0 , MAE : 0.014175


#   n1 = 198:
# QR AR(1): 16 - 1 - 3 , MAE: 0.00486
# CQR AR(1): 17 - 3 - 0 , MAE : 0.004535

# QR AR(2): 14 - 5 - 1 , MAE: 0.00659
# CQR AR(2): 4 - 16 - 0 , MAE : 0.013395


#   n1 = 998:
# QR AR(1): 19 - 1 - 0 , MAE : 0.00329
# CQR AR(1): 15 - 5 - 0 , MAE : 0.006015    MAE AUMENTA RISPETTO A N1 = 198!?!?!?!?!?!?!

# QR AR(2): 19 - 1 - 0 , MAE : 0.002145
# CQR AR(2): 19 - 0 - 1 , MAE :  0.002845
           


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


# da qui comincia la simulazione

run_simulation <- function(){
  
  n <- 201 #3 elements will be lost in the DGP
  n2 <- 100

#-------------------Generate n1 points for the AR(2) model, this is our DGP
  
  phi_ar2 <- c(0.5, -0.2)  # AR coefficients for AR(2)
  Y_ar2 <- numeric(2*n)
  Y_ar2[1] <- 1
  Y_ar2[2] <- -2
  for (i in 3:(2*n)) {
    Y_ar2[i] <- phi_ar2[1] * Y_ar2[i-1] + phi_ar2[2] * Y_ar2[i-2] + rnorm(1,0,1)            #rt(n = 1, df = 2) #heavy tails
  }
  Y_ar2 <- Y_ar2[(n + 1): (2*n)] #I implemented a burn-in

  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])
  Y_lag3_ar2 <- c(NA, NA, NA, Y_ar2[1:(length(Y_ar2)-3)])
  
  Y_test <- phi_ar2[1] * Y_ar2[n] + phi_ar2[2] * Y_ar2[n-1] + rnorm(n2,0,1) #rt(n = n2, df = 2) #I create n2 = 100 points for testing
  
  # Prepare the dataset for AR(3) (excluding the first three NA values)
  data_ar2 <- data.frame(I = 1, Y = Y_ar2[-c(1:3)], Y_lag1 = Y_lag1_ar2[-c(1:3)], Y_lag2 = Y_lag2_ar2[-c(1:3)], Y_lag3 = Y_lag3_ar2[-c(1:3)])
  

  
#--------------- QUANTILE REGRESSION 
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)

  QuantAR1_OOS <- matrix(NA, length(QQ))
  QuantAR2_OOS <- matrix(NA, length(QQ))
  QuantAR3_OOS <- matrix(NA, length(QQ))
  
    for (jq in 1:length(QQ)) {  
      
      QR1 <- rq(Y ~ Y_lag1, data = data_ar2, tau=QQ[jq])       

      QuantAR1_OOS[jq] <- (c(1, Y_ar2[n])) %*% coef(QR1)
      
      QR2 <- rq(Y ~ Y_lag1 + Y_lag2, data = data_ar2, tau=QQ[jq])
      QuantAR2_OOS[jq] <- (c(1, Y_ar2[n], Y_ar2[n-1])) %*% coef(QR2)
       
      # QR3 <- rq(Y ~ Y_lag1 + Y_lag2 + Y_lag3, data = data_ar2, tau=QQ[jq])
      # QuantAR3_OOS[jq] <- (c(1, Y_ar2[n], Y_ar2[n-1], Y_ar2[n-2])) %*% coef(QR3)
      
    }
  
  

  
  
#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  CQuantAR1_OOS <- matrix(NA, length(QQ))
  CQuantAR2_OOS <- matrix(NA, length(QQ))
  CQuantAR3_OOS <- matrix(NA, length(QQ))
  Q_low <- matrix(NA, nrow(data_ar2), length(QQ))
  Q_high <- matrix(NA, nrow(data_ar2), length(QQ))
  full_length <- nrow(data_ar2)
  test_length = floor(full_length*66/100)
  data_ar2_1 <- data_ar2[1:test_length,] #train
  data_ar2_2 <- data_ar2[(test_length + 1) : full_length,] #calibration


    for (jq in 1:length(QQ)) {  
      
      #AR(1)
      Q_low[1:nrow(data_ar2_2), jq] <- -Inf 
      QR <- rq(Y ~ Y_lag1, data = data_ar2_1, tau= QQ[jq])
      Q_high[1 : nrow(data_ar2_2), jq] <- as.matrix(data_ar2_2[,-c(2,4,5)]) %*% coef(QR) 

      # Initialize a vector for errors
      E_i <- rep(NA, nrow(data_ar2_2))
      
      # Calculate errors for each point in the test set I2
      for (i in 1:length(E_i)) {
        E_i[i] <- max(Q_low[i, jq] - data_ar2_2[i,2], data_ar2_2[i,2] - Q_high[i, jq])
      }

      
      # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
      quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(data_ar2_2)))
      
      CQuantAR1_OOS[jq] <- (c(1, Y_ar2[n])) %*% coef(QR) + quantile_E
      
      
      #AR(2)
       Q_low[1:nrow(data_ar2_2), jq] <- -Inf 
       QR <- rq(Y ~ Y_lag1 + Y_lag2, data = data_ar2_1, tau= QQ[jq])
       Q_high[1 : nrow(data_ar2_2), jq] <- as.matrix(data_ar2_2[,-c(2,5)]) %*% coef(QR) 

       # Initialize a vector for errors
       E_i <- rep(NA, nrow(data_ar2_2))
       
       # Calculate errors for each point in the test set I2
       for (i in 1:length(E_i)) {
         E_i[i] <- max(Q_low[i, jq] - data_ar2_2[i,2], data_ar2_2[i,2] - Q_high[i, jq])
       }
       
       # Compute Q(QQ[jq])(E, I2) N.B 1 - ?? = QQ[jq]
       quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(data_ar2_2)))
       
       CQuantAR2_OOS[jq] <- (c(1, Y_ar2[n], Y_ar2[n-1])) %*% coef(QR) + quantile_E
       
      
      #AR(3)
      # Q_low[1:nrow(data_ar2_2), jq] <- -Inf 
      # QR <- rq(Y ~ Y_lag1 + Y_lag2 + Y_lag3, data = data_ar2_1, tau= QQ[jq])
      # Q_high[1 : nrow(data_ar2_2), jq] <- as.matrix(data_ar2_2[,-2]) %*% coef(QR) 
      # Q_high[1 : nrow(data_ar2_2), jq] <- sort(Q_high[1 : nrow(data_ar2_2), jq])
      # # Initialize a vector for errors
      # E_i <- rep(NA, nrow(data_ar2_2))
      # 
      # # Calculate errors for each point in the test set I2
      # for (i in 1:length(E_i)) {
      #   E_i[i] <- max(Q_low[i, jq] - data_ar2_2[i,2], data_ar2_2[i,2] - Q_high[i, jq])
      # }
      # 
      # # Compute Q(QQ[jq])(E, I2) N.B 1 - ?? = QQ[jq]
      # quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/nrow(data_ar2_2)))))
      # 
      # CQuantAR3_OOS[jt,jq] <- (c(1, Y_ar2[n], Y_ar2[n-1], Y_ar2[n-1])) %*% coef(QR) + quantile_E
     
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
  library(quantreg)}
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
resultsPitSTQRAR1 <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

for (jq in 1:length(QQ)) {
  
  success_rate <- average_coverageQRAR1[jq]
  
  n <- n2*n3
  
  successes <- success_rate*n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsPitSTQRAR1 <- rbind(resultsPitSTQRAR1, data.frame(Quantile = QQ[jq], SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsPitSTQRAR1$IsWithinCI <- mapply(is_within_ci, resultsPitSTQRAR1$Quantile, resultsPitSTQRAR1$CI_Lower, resultsPitSTQRAR1$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsPitSTQRAR1$IsWithinCI) * 100

# Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsPitSTQRAR1$SuccessRate > resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI)

# Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsPitSTQRAR1$SuccessRate < resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))
print(paste("MAE of QR:", mean(abs(resultsPitSTQRAR1$Quantile-resultsPitSTQRAR1$SuccessRate))))
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
resultsPitSTCQRAR1 <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

for (jq in 1:length(QQ)) {
  
  success_rate <- average_coverageCQRAR1[jq]
  
  n <- n2*n3
  
  successes <- success_rate*n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsPitSTCQRAR1 <- rbind(resultsPitSTCQRAR1, data.frame(Quantile = QQ[jq], SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsPitSTCQRAR1$IsWithinCI <- mapply(is_within_ci, resultsPitSTCQRAR1$Quantile, resultsPitSTCQRAR1$CI_Lower, resultsPitSTCQRAR1$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsPitSTCQRAR1$IsWithinCI) * 100

# Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsPitSTCQRAR1$SuccessRate > resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI)

# Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsPitSTCQRAR1$SuccessRate < resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))
print(paste("MAE of CQR:", mean(abs(resultsPitSTCQRAR1$Quantile-resultsPitSTCQRAR1$SuccessRate))))
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





