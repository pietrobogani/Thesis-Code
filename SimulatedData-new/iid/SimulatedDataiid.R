# - che errore uso? per ora sto usando code pesanti
# - visto che è iid e senza componente temporale, ha senso prevedere 100 punti con covariate diverse? Se sì, questo è un problema per DCP
#   che ci metterà tantissimo


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
library(ggplot2)



# da qui comincia la simulazione

run_simulation <- function(){
  
  # Set parameters
  n <- 101       # Number of data points for training
  n2 <- 100      # Number of data points for testing
  p <- round(0.1* n)
  num_exog_vars <-  p     # Parameter to control the number of exogenous variables
  beta <- runif(num_exog_vars, 0.1, 1) # Randomly generated coefficients for exogenous variables
  means <- rnorm(p)
  sigma <- diag(runif(p))
  
  # Initialize variables
  Y <- numeric(n)  
  exog_vars <- rmvnorm( n = n + n2, mean = means, sigma = sigma)

  # Simulate i.i.d data with exogenous variables 
  for (t in 1:(n + n2)) {
    Y[t] <- sum(beta * exog_vars[t, ]) + rt(n = 1, df = 2)
  }
  

  # Prepare dataset for model 
  data_ar2 <- data.frame(I = 1, Y = Y[-c((n+1):(n+n2))])
  for (j in 1:num_exog_vars) {
    data_ar2[paste0("X", j)] <- exog_vars[-c((n+1):(n+n2)), j]
  }
  
  data_ar2_test <- data.frame(I = 1, Y = Y[c((n+1):(n+n2))])
  for (j in 1:num_exog_vars) {
    data_ar2_test[paste0("X", j)] <- exog_vars[c((n+1):(n+n2)), j]
  }
  
  Y_test <- Y[c((n+1):(n+n2))]
  
  
  acf_values <- acf(Y, plot = FALSE)
  
  # Step 2: Plot ACF
  plot(acf_values, main="Autocorrelation Function")
  
  
  
  
  
  
  # Function to generate the formula
  make_formula <- function(num_exog_vars) {
    # Create a character vector for the exogenous variable names
    exog_vars <- paste0("X", 1:num_exog_vars)
    # Construct the right-hand side of the formula
    rhs <- paste(exog_vars, collapse = " + ")
    # Combine to form the full formula as a string
    formula_str <- paste("Y ~", rhs)
    # Convert the string to a formula object
    formula_obj <- as.formula(formula_str)
    
    return(formula_obj)
  }
  
  source("temp1.R")
  source("full.R")
  funs=rq.funs()
  
  predict.fun=funs$predict.fun
  train.fun=funs$train.fun
  
#--------------- QUANTILE REGRESSION 
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)

  QuantAR_OOS <- matrix(NA, n2, length(QQ))
  #qrf_model <- quantregForest(x = as.matrix(data_ar2[,-c(1,2)]), y = as.matrix(data_ar2[,2]))
  

    for (jq in 1:length(QQ)) {  
       
       QR <- rq(make_formula(num_exog_vars), data = data_ar2, tau = QQ[jq])     
       QuantAR_OOS[,jq] <- as.matrix(data_ar2_test[,-2]) %*% coef(QR)
      
       #QuantAR_OOS[,jq] <- predict(qrf_model, newdata = as.matrix(data_ar2_test[,-c(1,2)]), what = QQ[jq])
       
      
    }

  
  
#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  CQuant_OOS <- matrix(NA, n2, length(QQ))
  Q_low <- matrix(NA, nrow(data_ar2), length(QQ))
  Q_high <- matrix(NA, nrow(data_ar2), length(QQ))
  full_length <- nrow(data_ar2)
  test_length = floor(full_length*50/100)
  data_ar2_1 <- data_ar2[1:test_length,] #train
  data_ar2_2 <- data_ar2[(test_length + 1) : full_length,] #calibration

  #qrf_model <- quantregForest(x = as.matrix(data_ar2_1[,-c(1,2)]), y = as.matrix(data_ar2_1[,2]))
  

    for (jq in 1:length(QQ)) {  
       
       Q_low[1:nrow(data_ar2_2), jq] <- -Inf 
       QR <- rq(make_formula(num_exog_vars), data = data_ar2_1, tau = QQ[jq])
       Q_high[1 : nrow(data_ar2_2), jq] <- as.matrix(data_ar2_2[,-2]) %*% coef(QR) 
        
       #Q_high[1 : nrow(data_ar2_2), jq] <- predict(qrf_model, newdata = as.matrix(data_ar2_2[,-c(1,2)]), what = QQ[jq])

       # Initialize a vector for errors
       E_i <- rep(NA, nrow(data_ar2_2))
       
       # Calculate errors for each point in the test set I2
       for (i in 1:length(E_i)) {
         E_i[i] <- max(Q_low[i, jq] - data_ar2_2[i,2], data_ar2_2[i,2] - Q_high[i, jq])
       }
       
       # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
       quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(data_ar2_2)))
       
       CQuant_OOS[,jq] <- as.matrix(data_ar2_test[,-2]) %*% coef(QR) + quantile_E
       #CQuant_OOS[,jq] <- predict(qrf_model, newdata = as.matrix(data_ar2_test[,-c(1,2)]), what = QQ[jq]) + quantile_E
     
    }
  

  
  #-------------------------------------------------------------------------------------------
  
  #CALIBRATION OF QR AND CQR
  
  # This function returns the cumulative probability for a given value X, assuming a step function as CDF
  coverageQR <- rep(NA,length(QQ))
  coverageCQR <- rep(NA,length(QQ))
  
  PitQROOS <- rep(NA, length(Y_test))
  PitCQROOS <- rep(NA, length(Y_test))
  
  
  for (jq in 1:length(QQ)) {
    
  for (i in 1:length(Y_test)){
    if (Y_test[i] <= QuantAR_OOS[i,jq]) {
      PitQROOS[i]  <- 1 
    }
    else 
      PitQROOS[i]  <- 0
  }
  
  coverageQR[jq] <- sum(PitQROOS) / length(PitQROOS)


  for (i in 1:length(Y_test)){
    if (Y_test[i] <= CQuant_OOS[i,jq]) {
      PitCQROOS[i]  <- 1 
    }
    else 
      PitCQROOS[i]  <- 0
}
  
  coverageCQR[jq] <- sum(PitCQROOS) / length(PitCQROOS)

  }
  
    
  ret <- list(coverageQR, coverageCQR)
  
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
coveragetotQR <- matrix(NA,n3,length(QQ))
coveragetotCQR <- matrix(NA,n3,length(QQ))



index = 1
for(res in results){
  coveragetotQR[index,] <- res[[1]]
  coveragetotCQR[index,] <- res[[2]]
  index <- index + 1
}



# RESULTS EVALUATION

# Calculate average coverage across all iterations and for each quantile
average_coverageQR <- rep(NA,length(QQ))
average_coverageCQR <- rep(NA,length(QQ))

for (jq in 1:length(QQ)) { 
  average_coverageQR[jq] <- mean(coveragetotQR[,jq])

  average_coverageCQR[jq] <- mean(coveragetotCQR[,jq])
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



#----------------------- QR 
{
 resultsPitSTQR <- data.frame(Quantile = numeric(), EmpiricalCoverage = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
 
 for (jq in 1:length(QQ)) {
   
   success_rate <- average_coverageQR[jq]
   
   n <- n2*n3
   
   successes <- success_rate*n
   
   # Calculate Wilson score interval
   ci <- wilson_score_interval(successes, n)
   
   # Add to results data frame
   resultsPitSTQR <- rbind(resultsPitSTQR, data.frame(Quantile = QQ[jq], EmpiricalCoverage = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
 }
 
 
 # Apply the function to each row and add the result as a new column
 resultsPitSTQR$IsWithinCI <- mapply(is_within_ci, resultsPitSTQR$Quantile, resultsPitSTQR$CI_Lower, resultsPitSTQR$CI_Upper)
 
 # Calculate the percentage of quantiles within the CI
 percentage_within_ci <- mean(resultsPitSTQR$IsWithinCI) * 100
 
 # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
 count_above <- sum(resultsPitSTQR$EmpiricalCoverage > resultsPitSTQR$Quantile & !resultsPitSTQR$IsWithinCI)
 
 # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
 count_below <- sum(resultsPitSTQR$EmpiricalCoverage < resultsPitSTQR$Quantile & !resultsPitSTQR$IsWithinCI)
 
 # Print the result
 print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
 print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
 print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
 print(paste("MAE of QR:", mean(abs(resultsPitSTQR$Quantile-resultsPitSTQR$EmpiricalCoverage))))
}

#----------------------- CQR 
{
 resultsPitSTCQR <- data.frame(Quantile = numeric(), EmpiricalCoverage = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
 
 for (jq in 1:length(QQ)) {
   
   success_rate <- average_coverageCQR[jq]
   
   n <- n2*n3
   
   successes <- success_rate*n
   
   # Calculate Wilson score interval
   ci <- wilson_score_interval(successes, n)
   
   # Add to results data frame
   resultsPitSTCQR <- rbind(resultsPitSTCQR, data.frame(Quantile = QQ[jq], EmpiricalCoverage = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
 }
 
 
 # Apply the function to each row and add the result as a new column
 resultsPitSTCQR$IsWithinCI <- mapply(is_within_ci, resultsPitSTCQR$Quantile, resultsPitSTCQR$CI_Lower, resultsPitSTCQR$CI_Upper)
 
 # Calculate the percentage of quantiles within the CI
 percentage_within_ci <- mean(resultsPitSTCQR$IsWithinCI) * 100
 
 # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
 count_above <- sum(resultsPitSTCQR$EmpiricalCoverage > resultsPitSTCQR$Quantile & !resultsPitSTCQR$IsWithinCI)
 
 # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
 count_below <- sum(resultsPitSTCQR$EmpiricalCoverage < resultsPitSTCQR$Quantile & !resultsPitSTCQR$IsWithinCI)
 
 # Print the result
 print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
 print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
 print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
 print(paste("MAE of CQR:", mean(abs(resultsPitSTCQR$Quantile-resultsPitSTCQR$EmpiricalCoverage))))
 }



{
  x1 <- resultsPitSTCQR$Quantile
  y1 <- resultsPitSTCQR$EmpiricalCoverage
  x2 <- resultsPitSTQR$Quantile
  y2 <- resultsPitSTQR$EmpiricalCoverage
  
  # Create the plot
  plot(x1, y1, type = "n", xlab = "Quantile", ylab = "Empirical Coverage", main = "Staircase Plots with Diagonal Line")
  
  # Add the staircase lines for the first dataset
  for (i in 2:length(x1)) {
    segments(x1[i-1], y1[i-1], x1[i-1], y1[i], col = "blue") # Vertical segment
    segments(x1[i-1], y1[i], x1[i], y1[i], col = "blue")     # Horizontal segment
  }
  
  # Add the points for the first dataset
  points(x1, y1)
  
  # Add the staircase lines for the second dataset
  for (i in 2:length(x2)) {
    segments(x2[i-1], y2[i-1], x2[i-1], y2[i], col = "red") # Vertical segment
    segments(x2[i-1], y2[i], x2[i], y2[i], col = "red")     # Horizontal segment
  }
  
  # Add the points for the second dataset
  points(x2, y2, col = "red")
  
  # Add the diagonal line
  abline(a = 0, b = 1, col = "black", lty = 2)
  
  # Optional: Add a legend
  legend("bottomright", legend = c("CQR AR1", "QR AR1"), col = c("blue", "red"), pch = c(16, 16), bty = "n", cex = 1.5, pt.cex = 1.5)
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





