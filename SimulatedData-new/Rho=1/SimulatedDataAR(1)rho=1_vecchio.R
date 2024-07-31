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
library(quantregForest)




run_simulation <- function(){
  
  n <- 1001 #3 elements will be lost in the DGP
  n2 <- 100
  
  source("temp1.R")
  source("full.R")
  funs=rq.funs()
  
  predict.fun=funs$predict.fun
  train.fun=funs$train.fun

#-------------------Generate n1 points for the AR(1) model, this is our DGP
  
  phi_ar2 <- 1.05  # AR coefficients for AR(1)
  Y_ar2 <- numeric(n)
  Y_ar2[1] <- 1
  for (i in 2:n) {
    Y_ar2[i] <- phi_ar2 * Y_ar2[i-1] + rt(n = 1, df = 2) #heavy tails
  }
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])

  Y_test <- phi_ar2 * Y_ar2[n] + rt(n = n2, df = 2) #I create n2 = 100 points for testing
  
  # Prepare the dataset for AR(3) (excluding the first three NA values)
  data_ar2 <- data.frame(I = 1, Y = Y_ar2[-c(1:3)], Y_lag1 = Y_lag1_ar2[-c(1:3)])
  

  
#--------------- QUANTILE REGRESSION 
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)
  
  qrf_model <- quantregForest(x = data_ar2[, 3], y = data_ar2[,2])
  

  QuantAR1_OOS <- matrix(NA, length(QQ))

    for (jq in 1:length(QQ)) {  
      
        #QR1 <- rq(Y ~ Y_lag1, data = data_ar2, tau=QQ[jq])
        #QuantAR1_OOS[jq] <- (c(1, Y_ar2[n])) %*% coef(QR1)
      
      QuantAR1_OOS[jq] <- predict(qrf_model, newdata = as.matrix(rep(Y_ar2[n],n2)), what = QQ[jq])[1]
      
    }
  

  
  
#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  CQuantAR1_OOS <- matrix(NA, length(QQ))
  full_length <- nrow(data_ar2)
  test_length = floor(full_length*50/100)
  data_ar2_1 <- data_ar2[1:test_length,] #train
  data_ar2_2 <- data_ar2[(test_length + 1) : full_length,] #calibration
  Q_low <- matrix(NA, nrow(data_ar2_2), length(QQ))
  Q_high <- matrix(NA, nrow(data_ar2_2), length(QQ))
  
  coverageCQRAR1 <- rep(NA,length(QQ))
  qrf_model <- quantregForest(x = data_ar2_1[, 3], y = data_ar2_1[,2])
  


    for (jq in 1:length(QQ)) {  
      
      #AR(1)
       Q_low[1:nrow(data_ar2_2), jq] <- -Inf 

       #QR1 <- rq(Y ~ Y_lag1, data = data_ar2_1, tau=QQ[jq])
      #Q_high[,jq] <- as.matrix(data_ar2_2[, c(1,3)]) %*% coef(QR1)
       Q_high[,jq] <- predict(qrf_model, newdata = as.matrix(data_ar2_2[, 3]), what = QQ[jq])
       
       # Initialize a vector for errors
       E_i <- rep(NA, nrow(data_ar2_2))
       
       # Calculate errors for each point in the test set I2
       for (i in 1:length(E_i)) {
         E_i[i] <- max(Q_low[i, jq] - data_ar2_2[i,2], data_ar2_2[i,2] - Q_high[i, jq])
       }
       
       # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
       quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(data_ar2_2)))
       
       #CQuantAR1_OOS[jq] <-  (c(1, Y_ar2[n])) %*% coef(QR1) + quantile_E
       CQuantAR1_OOS[jq] <- predict(qrf_model, newdata = as.matrix(rep(Y_ar2[n],n2)), what = QQ[jq])[1] + quantile_E

    }
  
  
  #--------------- DISTRIBUTIONAL CONFORMAL PREDICTION
   
   DQuantAR1_OOS <- matrix(0, length(QQ))
   
   
    # y <- data_ar2[,2]
    # x <- data_ar2[,3]
    # x0 <- Y_ar2[n]
    # 
    #  for (jq in 1:length(QQ)) {
    #    test = dist.conformal.pred(x,y,x0,train.fun = train.fun, predict.fun = predict.fun, verbose = T, alpha = 1 - QQ[jq])
    #    DQuantAR1_OOS[jq] <- test$up[1,1]
    #    cat("Upper level of CI:", test$up[1, 1], "- Lower level of CI:", test$lo[1, 1], "\n")
    #    cat("Minimum of vector Y:", min(y), "\n")
    #  }
    # 
    # 
  #-------------------------------------------------------------------------------------------
  
  #CALIBRATION OF QR AND CQR
  
  # This function returns the cumulative probability for a given value X, assuming a step function as CDF
  coverageQRAR1 <- rep(NA,length(QQ))
  coverageDCPAR1 <- rep(NA,length(QQ))
  PitDCPOOSAR1 <- rep(NA, length(Y_test))
  
  

  
  for (jq in 1:length(QQ)) {
  PitQROOSAR1 <- rep(NA, length(Y_test))

  PitCQROOSAR1 <- rep(NA, length(Y_test))
  
  PitDCPOOSAR1 <- rep(NA, length(Y_test))
  

  for (i in 1:length(Y_test)){
    if (Y_test[i] <= QuantAR1_OOS[jq]) {
      PitQROOSAR1[i]  <- 1 
    }
    else 
      PitQROOSAR1[i]  <- 0
  }
  
  coverageQRAR1[jq] <- sum(PitQROOSAR1) / length(PitQROOSAR1)


   for (i in 1:length(Y_test)){
     if (Y_test[i] <= CQuantAR1_OOS[jq]) {
       PitCQROOSAR1[i]  <- 1 
     }
     else 
       PitCQROOSAR1[i]  <- 0
   }
  
  coverageCQRAR1[jq] <- sum(PitCQROOSAR1) / length(PitCQROOSAR1)
  
  
   for (i in 1:length(Y_test)){
     if (Y_test[i] <= DQuantAR1_OOS[jq]) {
       PitDCPOOSAR1[i]  <- 1 
     }
     else 
       PitDCPOOSAR1[i]  <- 0
   }
   
  coverageDCPAR1[jq] <- sum(PitDCPOOSAR1) / length(PitDCPOOSAR1)

  }
  
    
  ret <- list(coverageQRAR1, coverageCQRAR1,coverageDCPAR1)
  
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
  library(quantregForest)}
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
coveragetotCQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotDCPAR1 <- matrix(NA,n3,length(QQ))




index = 1
for(res in results){
  coveragetotQRAR1[index,] <- res[[1]]
  coveragetotCQRAR1[index,] <- res[[2]]
  coveragetotDCPAR1[index,] <- res[[3]]
  
  index <- index + 1
}



# RESULTS EVALUATION

# Calculate average coverage across all iterations and for each quantile
average_coverageQRAR1 <- rep(NA,length(QQ))

average_coverageDCPAR1 <- rep(NA,length(QQ))

average_coverageCQRAR1 <- rep(NA,length(QQ))

for (jq in 1:length(QQ)) { 
  average_coverageQRAR1[jq] <- mean(coveragetotQRAR1[,jq])

  average_coverageCQRAR1[jq] <- mean(coveragetotCQRAR1[,jq])
  
  average_coverageDCPAR1[jq] <- mean(coveragetotDCPAR1[,jq])
  
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
resultsPitSTQRAR1 <- data.frame(Quantile = numeric(), EmpiricalCoverage = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

for (jq in 1:length(QQ)) {
  
  success_rate <- average_coverageQRAR1[jq]
  
  n <- n2*n3
  
  successes <- success_rate*n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsPitSTQRAR1 <- rbind(resultsPitSTQRAR1, data.frame(Quantile = QQ[jq], EmpiricalCoverage = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsPitSTQRAR1$IsWithinCI <- mapply(is_within_ci, resultsPitSTQRAR1$Quantile, resultsPitSTQRAR1$CI_Lower, resultsPitSTQRAR1$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsPitSTQRAR1$IsWithinCI) * 100

# Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsPitSTQRAR1$EmpiricalCoverage > resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI)

# Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsPitSTQRAR1$EmpiricalCoverage < resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
print(paste("MAE of QR:", mean(abs(resultsPitSTQRAR1$Quantile-resultsPitSTQRAR1$EmpiricalCoverage))))
}


#----------------------- CQR AR(1)
{
resultsPitSTCQRAR1 <- data.frame(Quantile = numeric(), EmpiricalCoverage = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

for (jq in 1:length(QQ)) {
  
  success_rate <- average_coverageCQRAR1[jq]
  
  n <- n2*n3
  
  successes <- success_rate*n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsPitSTCQRAR1 <- rbind(resultsPitSTCQRAR1, data.frame(Quantile = QQ[jq], EmpiricalCoverage = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsPitSTCQRAR1$IsWithinCI <- mapply(is_within_ci, resultsPitSTCQRAR1$Quantile, resultsPitSTCQRAR1$CI_Lower, resultsPitSTCQRAR1$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsPitSTCQRAR1$IsWithinCI) * 100

# Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsPitSTCQRAR1$EmpiricalCoverage > resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI)

# Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsPitSTCQRAR1$EmpiricalCoverage < resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
print(paste("MAE of CQR:", mean(abs(resultsPitSTCQRAR1$Quantile-resultsPitSTCQRAR1$EmpiricalCoverage))))
}



#----------------------- DCP AR(1)

{
  resultsPitSTDCPAR1 <- data.frame(Quantile = numeric(), EmpiricalCoverage = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
  
  for (jq in 1:length(QQ)) {
    
    success_rate <- average_coverageDCPAR1[jq]
    
    n <- n2*n3
    
    successes <- success_rate*n
    
    # Calculate Wilson score interval
    ci <- wilson_score_interval(successes, n)
    
    # Add to results data frame
    resultsPitSTDCPAR1 <- rbind(resultsPitSTDCPAR1, data.frame(Quantile = QQ[jq], EmpiricalCoverage = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
  }
  
  
  # Apply the function to each row and add the result as a new column
  resultsPitSTDCPAR1$IsWithinCI <- mapply(is_within_ci, resultsPitSTDCPAR1$Quantile, resultsPitSTDCPAR1$CI_Lower, resultsPitSTDCPAR1$CI_Upper)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTDCPAR1$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTDCPAR1$EmpiricalCoverage > resultsPitSTDCPAR1$Quantile & !resultsPitSTDCPAR1$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTDCPAR1$EmpiricalCoverage < resultsPitSTDCPAR1$Quantile & !resultsPitSTDCPAR1$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for DCP:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of DCP:", mean(abs(resultsPitSTDCPAR1$Quantile-resultsPitSTDCPAR1$EmpiricalCoverage))))
}





#------------- PLOT CALIBRATION CURVES OF QR and CQR and DCP with AR(1) MODEL

{
  x1 <- resultsPitSTCQRAR1$Quantile
  y1 <- resultsPitSTCQRAR1$EmpiricalCoverage
  x2 <- resultsPitSTQRAR1$Quantile
  y2 <- resultsPitSTQRAR1$EmpiricalCoverage
  x3 <- resultsPitSTDCPAR1$Quantile
  y3 <- resultsPitSTDCPAR1$EmpiricalCoverage
  
  # Create the plot
  plot(x1, y1, type = "n", xlab = "Quantile", ylab = "Empirical Coverage", main = "Staircase Plots with Diagonal Line")
  
  # Add the staircase lines for the first dataset
  for (i in 2:length(x1)) {
    segments(x1[i-1], y1[i-1], x1[i-1], y1[i], col = "blue") # Vertical segment
    segments(x1[i-1], y1[i], x1[i], y1[i], col = "blue")     # Horizontal segment
  }
  
  # Add the points for the first dataset
  points(x1, y1, col = "blue")
  
  # Add the staircase lines for the second dataset
  for (i in 2:length(x2)) {
    segments(x2[i-1], y2[i-1], x2[i-1], y2[i], col = "red") # Vertical segment
    segments(x2[i-1], y2[i], x2[i], y2[i], col = "red")     # Horizontal segment
  }
  
  # Add the points for the second dataset
  points(x2, y2, col = "red")
  
  # Add the staircase lines for the third dataset
  for (i in 2:length(x3)) {
    segments(x3[i-1], y3[i-1], x3[i-1], y3[i], col = "green") # Vertical segment
    segments(x3[i-1], y3[i], x3[i], y3[i], col = "green")     # Horizontal segment
  }
  
  # Add the points for the third dataset
  points(x3, y3, col = "green")
  
  # Add the diagonal line
  abline(a = 0, b = 1, col = "black", lty = 2)
  
  # Optional: Add a legend
  legend("bottomright", legend = c("CQR AR1", "QR AR1", "DCP AR1"), col = c("blue", "red", "green"), pch = c(16, 17, 18), bty = "n", cex = 1.5, pt.cex = 1.5)
}















 filename <- paste("rho95dcp98", ".RData",sep="")
 cat(paste("Saving results to file", filename, "\n"))
# 
 # Save all the variables to the .RData file
 save(
   resultsPitSTDCPAR1,
   resultsPitSTCQRAR1,
   resultsPitSTQRAR1,   
   file=filename
 )
 


# #-------------------------------- per caricare dati salvati -----------------------
# 
# filename <- paste("Sim500_1", ".RData", sep="")
# cat(paste("Loading results from file", filename, "\n"))
# 
# load(filename)
# 
# #-----------------------------------------------------------------------------------





