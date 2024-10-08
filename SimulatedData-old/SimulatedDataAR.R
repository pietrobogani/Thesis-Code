#PROBLEMI. 
#1) Se tolgo il noise, ho collinearità perfetta e quindi QR fallisce
#2) Se aggiungo noise peṛ non riesco ad avere esattamento 0 errore



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
  
  n <- 1000
  jtFirstOOS <- 800 #First index for out-of-sample computations
  
  
#-------------------Generate AR(2) model, this is our DGP
  
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
  

  
#--------------- QUANTILE REGRESSION 
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ, 0.99)

  QuantAR1_OOS <- matrix(NA, nrow(data_ar2), length(QQ))
  QuantAR2_OOS <- matrix(NA, nrow(data_ar2), length(QQ))
  QuantAR3_OOS <- matrix(NA, nrow(data_ar2), length(QQ))
  
  for (jt in jtFirstOOS:nrow(data_ar2)) {
    for (jq in 1:length(QQ)) {  
      
      QR1 <- rq(Y[1:(jt-1)] ~ Y_lag1[1:(jt - 1)], data = data_ar2, tau=QQ[jq])
      QuantAR1_OOS[jt, jq] <- as.matrix(data_ar2[jt, - c(2,4,5) ]) %*% coef(QR1)
      
      QR2 <- rq(Y[1:(jt-1)] ~ Y_lag1[1:(jt - 1)] + Y_lag2[1:(jt - 1)], data = data_ar2, tau=QQ[jq])
      QuantAR2_OOS[jt, jq] <- as.matrix(data_ar2[jt, - c(2,5) ]) %*% coef(QR2)
       
      # QR3 <- rq(Y[1:(jt-1)] ~ Y_lag1[1:(jt - 1)] + Y_lag2[1:(jt - 1)] + Y_lag3[1:(jt - 1)], data = data_ar2, tau=QQ[jq])
      # QuantAR3_OOS[jt, jq] <- as.matrix(data_ar2[jt, - 2 ]) %*% coef(QR3)
      
    }
  }
  

  
  
#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  CQuantAR1_OOS <- matrix(NA, nrow(data_ar2), length(QQ))
  CQuantAR2_OOS <- matrix(NA, nrow(data_ar2), length(QQ))
  CQuantAR3_OOS <- matrix(NA, nrow(data_ar2), length(QQ))
  Q_low <- matrix(NA, nrow(data_ar2), length(QQ))
  Q_high <- matrix(NA, nrow(data_ar2), length(QQ))
  
  
  for (jt in jtFirstOOS:nrow(data_ar2)) {
    
    
    full_length <- nrow(data_ar2[1:(jt - 1),])
    test_length = floor(full_length*50/100)
    data_ar2_1 <- data_ar2[1:test_length,]
    data_ar2_2 <- data_ar2[(test_length + 1) : (jt - 1),]


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
      quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/nrow(data_ar2_2)))))
      
      CQuantAR1_OOS[jt,jq] <- as.matrix(data_ar2[jt, -c(2,4,5)]) %*% coef(QR) + quantile_E
      
      
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
       quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/nrow(data_ar2_2)))))
       
       CQuantAR2_OOS[jt,jq] <- as.matrix(data_ar2[jt, -c(2,5)]) %*% coef(QR) + quantile_E
       
      
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
      # CQuantAR3_OOS[jt,jq] <- as.matrix(data_ar2[jt, -2]) %*% coef(QR) + quantile_E
     
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
  PitSTAR2_OOS <- rep(NA, nrow(data_ar2))
  PitSTAR3_OOS <- rep(NA, nrow(data_ar2))
  
  CPitSTAR1_OOS <- rep(NA, nrow(data_ar2))
  CPitSTAR2_OOS <- rep(NA, nrow(data_ar2))
  CPitSTAR3_OOS <- rep(NA, nrow(data_ar2))
  

  
  for (jt in jtFirstOOS: nrow(data_ar2) ){
    
    YhRealized <- data_ar2[jt,2]

    PitSTAR1_OOS[jt] <- cumulative_prob(YhRealized, QQ, QuantAR1_OOS[jt, ])
    PitSTAR2_OOS[jt] <- cumulative_prob(YhRealized, QQ, QuantAR2_OOS[jt, ])
    # PitSTAR3_OOS[jt] <- cumulative_prob(YhRealized, QQ, QuantAR3_OOS[jt, ])
    
    CPitSTAR1_OOS[jt] <- cumulative_prob(YhRealized, QQ, CQuantAR1_OOS[jt, ])
    CPitSTAR2_OOS[jt] <- cumulative_prob(YhRealized, QQ, CQuantAR2_OOS[jt, ])
    # CPitSTAR3_OOS[jt] <- cumulative_prob(YhRealized, QQ, CQuantAR3_OOS[jt, ])
    
  }
  ret <- list(PitSTAR1_OOS, PitSTAR2_OOS, PitSTAR3_OOS, CPitSTAR1_OOS, CPitSTAR2_OOS, CPitSTAR3_OOS)
  
  return(ret)
}






n_simul <- 200
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
  set.seed(seed)
  run_simulation()
})


# Stop the cluster
stopCluster(cl)

# Extract results
PitSTAR1_OOStot <- numeric(0)
PitSTAR2_OOStot <- numeric(0)
PitSTAR3_OOStot <- numeric(0)
CPitSTAR1_OOStot <- numeric(0)
CPitSTAR2_OOStot <- numeric(0)
CPitSTAR3_OOStot <- numeric(0)




for(res in results){
  PitSTAR1_OOStot <- c(PitSTAR1_OOStot, na.omit(res[[1]]))
  PitSTAR2_OOStot <- c(PitSTAR2_OOStot,  na.omit(res[[2]]))
  PitSTAR3_OOStot <- c(PitSTAR3_OOStot, na.omit(res[[3]]))
  CPitSTAR1_OOStot <- c(CPitSTAR1_OOStot, na.omit(res[[4]]))
  CPitSTAR2_OOStot <- c(CPitSTAR2_OOStot, na.omit(res[[5]]))
  CPitSTAR3_OOStot <- c(CPitSTAR3_OOStot, na.omit(res[[6]]))

}


# 
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



#------------------- WILSON SCORES

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

# Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsPitST$SuccessRate > resultsPitST$Quantile & !resultsPitST$IsWithinCI)

# Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsPitST$SuccessRate < resultsPitST$Quantile & !resultsPitST$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for QR AR(1):", percentage_within_ci))
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))
#--------------------------------




resultsPitST <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

PitSTAR2_OOStot <- na.omit(PitSTAR2_OOStot)
for (quantile in c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)) {
  successes <- sum(PitSTAR2_OOStot < quantile)
  n <- length(PitSTAR2_OOStot)
  
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

# Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsPitST$SuccessRate > resultsPitST$Quantile & !resultsPitST$IsWithinCI)

# Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsPitST$SuccessRate < resultsPitST$Quantile & !resultsPitST$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for QR AR(2):", percentage_within_ci))
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))

#------------------------------------------------



resultsCPitST <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

CPitSTAR1_OOStot <- na.omit(CPitSTAR1_OOStot)
n <- length(CPitSTAR1_OOStot)

for (quantile in c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)) {
  successes <- sum(CPitSTAR1_OOStot < quantile)
  
  # Calculate success rate
  success_rate <- successes / n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsCPitST <- rbind(resultsCPitST, data.frame(Quantile = quantile, SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsCPitST$IsWithinCI <- mapply(is_within_ci, resultsCPitST$Quantile, resultsCPitST$CI_Lower, resultsCPitST$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsCPitST$IsWithinCI) * 100

# Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsCPitST$SuccessRate > resultsCPitST$Quantile & !resultsCPitST$IsWithinCI)

# Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsCPitST$SuccessRate < resultsCPitST$Quantile & !resultsCPitST$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for CQR AR(1):", percentage_within_ci))
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))

#-----------------------------------------------------


resultsCPitST <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

CPitSTAR2_OOStot <- na.omit(CPitSTAR2_OOStot)
n <- length(CPitSTAR2_OOStot)

for (quantile in c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)) {
  successes <- sum(CPitSTAR2_OOStot < quantile)
  
  # Calculate success rate
  success_rate <- successes / n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsCPitST <- rbind(resultsCPitST, data.frame(Quantile = quantile, SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsCPitST$IsWithinCI <- mapply(is_within_ci, resultsCPitST$Quantile, resultsCPitST$CI_Lower, resultsCPitST$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsCPitST$IsWithinCI) * 100

# Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsCPitST$SuccessRate > resultsCPitST$Quantile & !resultsCPitST$IsWithinCI)

# Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsCPitST$SuccessRate < resultsCPitST$Quantile & !resultsCPitST$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for CQR AR(2):", percentage_within_ci))
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))




























PITtest_env <- new.env()
source("PITtest.r",local = PITtest_env)
rstestboot_env <- new.env()
source("rstestboot.r",local = rstestboot_env)

# PROBABILITY INTEGRAL TRANSFORM (PIT) FOR QUANTILE REGRESSION, AR1
rvec <- seq(0, 1, by = 0.001)
zST_ecdf1 <- PITtest_env$PITtest(PitSTAR1_OOStot, rvec)


# PROBABILITY INTEGRAL TRANSFORM (PIT) FOR CONFORMAL QUANTILE REGRESSION, AR1

CzST_ecdf1 <- PITtest_env$PITtest(CPitSTAR1_OOStot, rvec)


# Plot PIT for quantile regression vs. conformal quantile regression 
plot(rvec, CzST_ecdf1, type = 'l', col = 'blue', xlab = 'tau', ylab = 'Empirical CDF', main = 'Full Conditional Model')
lines(rvec, zST_ecdf1, type = 'l', col = 'red')
lines(rvec, rvec , col = 'black',lty=2)

legend('bottomright', legend = c('Conformal Quantile Regression AR1', 'Quantile regression AR1'), cex = 1,fill = c('blue', 'red'))


# PROBABILITY INTEGRAL TRANSFORM (PIT) FOR QUANTILE REGRESSION, AR2
zST_ecdf2 <- PITtest_env$PITtest(PitSTAR2_OOStot, rvec)


# PROBABILITY INTEGRAL TRANSFORM (PIT) FOR CONFORMAL QUANTILE REGRESSION, AR2
CzST_ecdf2 <- PITtest_env$PITtest(CPitSTAR2_OOStot, rvec)


# Plot PIT for quantile regression vs. conformal quantile regression 
plot(rvec, CzST_ecdf2, type = 'l', col = 'blue', xlab = 'tau', ylab = 'Empirical CDF', main = 'Full Conditional Model')
lines(rvec, zST_ecdf2, type = 'l', col = 'red')
lines(rvec, rvec , col = 'black',lty=2)

legend('bottomright', legend = c('Conformal Quantile Regression AR2', 'Quantile regression AR2'), cex = 1,fill = c('blue', 'red'))




# PROBABILITY INTEGRAL TRANSFORM (PIT) FOR QUANTILE REGRESSION, AR3
zST_ecdf3 <- PITtest_env$PITtest(PitSTAR3_OOStot, rvec)


# PROBABILITY INTEGRAL TRANSFORM (PIT) FOR CONFORMAL QUANTILE REGRESSION, AR3
CzST_ecdf3 <- PITtest_env$PITtest(CPitSTAR3_OOStot, rvec)


# Plot PIT for quantile regression vs. conformal quantile regression 
plot(rvec, CzST_ecdf3, type = 'l', col = 'blue', xlab = 'tau', ylab = 'Empirical CDF', main = 'Full Conditional Model')
lines(rvec, zST_ecdf3, type = 'l', col = 'red')
lines(rvec, rvec , col = 'black',lty=2)

legend('bottomright', legend = c('Conformal Quantile Regression AR3', 'Quantile regression AR3'), cex = 1,fill = c('blue', 'red'))

















QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(0.01, QQ, 0.99)
quantAR1_est_cqr <- unique(CzST_ecdf1)[2: (length(unique(CzST_ecdf1))-1)]
quantAR1_est_qr <- unique(zST_ecdf1)[2: (length(unique(zST_ecdf1))-1)]

quantAR2_est_cqr <- unique(CzST_ecdf2)[2: (length(unique(CzST_ecdf2))-1)]
quantAR2_est_qr <- unique(zST_ecdf2)[2: (length(unique(zST_ecdf2))-1)]

quantAR3_est_cqr <- unique(CzST_ecdf3)[2: (length(unique(CzST_ecdf3))-1)]
quantAR3_est_qr <- unique(zST_ecdf3)[2: (length(unique(zST_ecdf3))-1)]



# Compute RMSE
rmseAR1_cqr  <- sqrt(mean((quantAR1_est_cqr  - QQ)^2))
rmseAR1_qr  <- sqrt(mean((quantAR1_est_qr  - QQ)^2))

rmseAR2_cqr  <- sqrt(mean((quantAR2_est_cqr  - QQ)^2))
rmseAR2_qr  <- sqrt(mean((quantAR2_est_qr  - QQ)^2))

rmseAR3_cqr  <- sqrt(mean((quantAR3_est_cqr  - QQ)^2))
rmseAR3_qr  <- sqrt(mean((quantAR3_est_qr  - QQ)^2))

# Compare RMSE and MAE
list(RMSEAR1_cqr  = rmseAR1_cqr, RMSEAR1_qr  = rmseAR1_qr, RMSEAR2_cqr  = rmseAR2_cqr, RMSEAR2_qr  = rmseAR2_qr,RMSEAR3_cqr  = rmseAR3_cqr, RMSEAR3_qr  = rmseAR3_qr )




#   #  T <- 204   jtFirstOOS <- 80 simul 200  
#[1] "Percentage of quantiles within the confidence interval CQR:  23.8095238095238 " above 14, below 2 RMSE_cqr 0.01827088
#[1] "Percentage of quantiles within the confidence interval QR:  33.3333333333333 " above 8 below 6 RMSE_qr 0.01365091

#  T <- 204   jtFirstOOS <- 80 simul 500  
#[1] "Percentage of quantiles within the confidence interval CQR: 14.2857142857143 above 15, below 3 RMSE_cqr 0.01689259
#[1] "Percentage of quantiles within the confidence interval QR: 19.047619047619" above 9, below 8 RMSE_qr 0.01460333

#  T <- 1004   jtFirstOOS <- 800 simul 200  
# "Percentage of quantiles within the confidence interval CQR: 23.8095238095238 above 9, below 7 RMSE_cqr 0.01101005
# "Percentage of quantiles within the confidence interval QR: 28.5714285714286  above 8, below 7 RMSE_qr 0.01379268

#  T <- 1004   jtFirstOOS <- 800 simul 500  
#[1] "Percentage of quantiles within the confidence interval CQR:  23.8095238095238 above 10, below 6 RMSE_cqr 0.008377895
#[1] "Percentage of quantiles within the confidence interval QR: 23.8095238095238 above 11, below 5  RMSE_qr 0.01045755

# T <- 2004 ,  jtFirstOOS <- 1800 #simul 200  
#[1] "Percentage of quantiles within the confidence interval CQR: 
#[1] "Percentage of quantiles within the confidence interval QR: 

# T <- 2004 ,  jtFirstOOS <- 1800 #simul 500  
#[1] "Percentage of quantiles within the confidence interval CQR: 
#[1] "Percentage of quantiles within the confidence interval QR: 

#T <- 3004 #time length  jtFirstOOS <- 2800 simul 200   
#[1] "Percentage of quantiles within the confidence interval CQR: 
#[1] "Percentage of quantiles within the confidence interval QR: 

#T <- 3004 #time length  jtFirstOOS <- 2800 simul 500  
#[1] "Percentage of quantiles within the confidence interval CQR: 4.76190476190476 above 12, below 8 RMSE_cqr 0.01146384
#[1] "Percentage of quantiles within the confidence interval QR: 0 above 12, below 9 RMSE_qr 0.01607059


