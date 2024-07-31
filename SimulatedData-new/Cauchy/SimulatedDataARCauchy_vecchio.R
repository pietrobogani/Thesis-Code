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
library(ggplot2)
source("functions.R")

run_simulation <- function(){
  
  n <- 101 #3 elements will be lost in the DGP
  n2 <- 100
  
  #-------------------Generate n1 points for the AR(2) model, this is our DGP
  
  phi_ar2 <- c(0.5, -0.2)  # AR coefficients for AR(2)
  Y_ar2 <- numeric(2*n)
  Y_ar2[1] <- 1
  Y_ar2[2] <- -2
  for (i in 3:(2*n)) {
    Y_ar2[i] <- phi_ar2[1] * Y_ar2[i-1] + phi_ar2[2] * Y_ar2[i-2] + rt(n = 1, df = 1) #CAUCHY ERRORS
  }
  Y_ar2 <- Y_ar2[(n + 1): (2*n)] #I implemented a burn-in
  
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])
  Y_lag3_ar2 <- c(NA, NA, NA, Y_ar2[1:(length(Y_ar2)-3)])
  
  Y_test <- phi_ar2[1] * Y_ar2[n] + phi_ar2[2] * Y_ar2[n-1] + rt(n = n2, df = 1) #I create n2 = 100 points for testing
  
  # Prepare the dataset (excluding the first three NA values)
  data_ar2 <- data.frame(I = 1, Y = Y_ar2[-c(1:3)], Y_lag1 = Y_lag1_ar2[-c(1:3)], Y_lag2 = Y_lag2_ar2[-c(1:3)], Y_lag3 = Y_lag3_ar2[-c(1:3)])
  
  QQ <- seq(0.005,0.965, by = 0.005) #vector of quantiles I'll do quantile regression on
  #QQ <- c(0.01, QQ)
  
  #--------------- QUANTILE REGRESSION 
  
  QuantAR1_OOS <- matrix(NA, length(QQ))
  QuantAR2_OOS <- matrix(NA, length(QQ))
  QuantAR3_OOS <- matrix(NA, length(QQ))
  
  {
    for (jq in 1:length(QQ)) {  
      QR1 <- rq(Y ~ Y_lag1, data = data_ar2, tau=QQ[jq])
      QuantAR1_OOS[jq] <- (c(1, Y_ar2[n])) %*% coef(QR1)
      
      QR2 <- rq(Y ~ Y_lag1 + Y_lag2, data = data_ar2, tau=QQ[jq])
      QuantAR2_OOS[jq] <- (c(1, Y_ar2[n], Y_ar2[n-1])) %*% coef(QR2)
      
      QR3 <- rq(Y ~ Y_lag1 + Y_lag2 + Y_lag3, data = data_ar2, tau=QQ[jq])
      QuantAR3_OOS[jq] <- (c(1, Y_ar2[n], Y_ar2[n-1], Y_ar2[n-2])) %*% coef(QR3)
    }
  }
  
  #--------------- QUANTILE RANDOM FOREST 
  
  {
    # qrf_model1 <- quantregForest(x = as.matrix(data_ar2[, 3]), y = data_ar2[,2]) #random forest AR1
    # qrf_model2 <- quantregForest(x = data_ar2[, c(3,4)], y = data_ar2[,2]) #random forest AR2
    # qrf_model3 <- quantregForest(x = data_ar2[, c(3,4,5)], y = data_ar2[,2]) #random forest AR3
    # 
    # for (jq in 1:length(QQ)) {  
    #   QuantAR1_OOS[jq] <- predict(qrf_model1, newdata = as.matrix(Y_ar2[n]), what = QQ[jq])
    #   QuantAR2_OOS[jq] <- predict(qrf_model2, newdata = c(Y_ar2[n],Y_ar2[n-1]), what = QQ[jq])
    #   QuantAR3_OOS[jq] <- predict(qrf_model3, newdata = c(Y_ar2[n],Y_ar2[n-1], Y_ar2[n-2]), what = QQ[jq])
    # }
  }
  
  #--------------- CONFORMALIZED QUANTILE REGRESSION
  
  full_length <- nrow(data_ar2)
  test_length = floor(full_length*50/100)
  data_ar2_1 <- data_ar2[1:test_length,] #train
  data_ar2_2 <- data_ar2[(test_length + 1) : full_length,] #calibration
  
  {
    x0 <- data_ar2_1[,-c(2,4,5)]
    y0 <- data_ar2_1[,2]
    x1 <- data_ar2_2[,-c(2,4,5)]
    y1 <- data_ar2_2[,2]
    x_test <- c(1,Y_ar2[n])
    
    CQuantAR1_OOS <- qCQR(Y ~ Y_lag1, x0, y0, x1, y1, x_test, QQ) 
    
    x0 <- data_ar2_1[,-c(2,5)]
    y0 <- data_ar2_1[,2]
    x1 <- data_ar2_2[,-c(2,5)]
    y1 <- data_ar2_2[,2]
    x_test <- c(1,Y_ar2[n],Y_ar2[n-1])
    
    CQuantAR2_OOS <- qCQR(Y ~ Y_lag1 + Y_lag2, x0, y0, x1, y1, x_test, QQ) 
    
    x0 <- data_ar2_1[,-2]
    y0 <- data_ar2_1[,2]
    x1 <- data_ar2_2[,-2]
    y1 <- data_ar2_2[,2]
    x_test <- c(1,Y_ar2[n],Y_ar2[n-1],Y_ar2[n-2])
    
    CQuantAR3_OOS <- qCQR(Y ~ Y_lag1 + Y_lag2 + Y_lag3, x0, y0, x1, y1, x_test, QQ) 
  }
  
  #--------------- CONFORMALIZED QUANTILE RANDOM FOREST
  
  {
    # x0 <- as.matrix(data_ar2_1[,3])
    # y0 <- data_ar2_1[,2]
    # x1 <- as.matrix(data_ar2_2[,3])
    # y1 <- data_ar2_2[,2]
    # x_test <- as.matrix(Y_ar2[n])
    # 
    # CQuantAR1_OOS_copy <- qCQRF(x0, y0, x1, y1, x_test, QQ) 
    # 
    # x0 <- data_ar2_1[,-c(1,2,5)]
    # y0 <- data_ar2_1[,2]
    # x1 <- data_ar2_2[,-c(1,2,5)]
    # y1 <- data_ar2_2[,2]
    # x_test <- c(Y_ar2[n],Y_ar2[n-1])
    # 
    # CQuantAR2_OOS <- qCQRF(x0, y0, x1, y1, x_test, QQ) 
    # 
    # x0 <- data_ar2_1[,-c(1,2)]
    # y0 <- data_ar2_1[,2]
    # x1 <- data_ar2_2[,-c(1,2)]
    # y1 <- data_ar2_2[,2]
    # x_test <- c(Y_ar2[n],Y_ar2[n-1],Y_ar2[n-2])
    # 
    # CQuantAR3_OOS <- qCQRF(x0, y0, x1, y1, x_test, QQ) 
  }
  
  #--------------- DISTRIBUTIONAL CONFORMAL PREDICTION
  
  {
    DQuantAR2_OOS <- matrix(0, length(QQ))
    
    # source("temp1.R")
    # source("full.R")
    # funs <- rq.funs()
    # 
    # predict.fun <- funs$predict.fun
    # train.fun <- funs$train.fun
    
    # y <- data_ar2[,2]
    # x <- data_ar2[,c(3,4)]
    # x0 <- c(Y_ar2[n], Y_ar2[n-1])
    
    # for (jq in 1:length(QQ)) {
    #   test <- dist.conformal.pred(x, y, x0, train.fun = train.fun, predict.fun = predict.fun, verbose = TRUE, alpha = 1 - QQ[jq])
    #   DQuantAR2_OOS[jq] <- test$up[1, 1]
    #   cat("Upper level of CI:", test$up[1, 1], "- Lower level of CI:", test$lo[1, 1], "\n")
    #   cat("Minimum of vector Y:", min(y), "\n")
    # }
  }
  
  #--------------- CALIBRATION

  coverageQRAR1 <- compute_coverage(Y_test, QuantAR1_OOS, QQ)
  coverageQRAR2 <- compute_coverage(Y_test, QuantAR2_OOS, QQ)
  coverageQRAR3 <- compute_coverage(Y_test, QuantAR3_OOS, QQ)
  
  coverageCQRAR1 <- compute_coverage(Y_test, CQuantAR1_OOS, QQ)
  coverageCQRAR2 <- compute_coverage(Y_test, CQuantAR2_OOS, QQ)
  coverageCQRAR3 <- compute_coverage(Y_test, CQuantAR3_OOS, QQ)
  
  coverageDCPAR2 <- compute_coverage(Y_test, DQuantAR2_OOS, QQ)
  
  
  #------------------- Calculation of CRPS
  
  CRPSQRAR1 <- calculate_crps_from_quantiles(100, QuantAR1_OOS, QQ)
  CRPSQRAR2 <- calculate_crps_from_quantiles(Y_test, QuantAR2_OOS, QQ)
  CRPSQRAR3 <- calculate_crps_from_quantiles(Y_test, QuantAR3_OOS, QQ)
  
  CRPSCQRAR1 <- calculate_crps_from_quantiles(100, CQuantAR1_OOS, QQ)
  CRPSCQRAR2 <- calculate_crps_from_quantiles(Y_test, CQuantAR2_OOS, QQ)
  CRPSCQRAR3 <- calculate_crps_from_quantiles(Y_test, CQuantAR3_OOS, QQ)
  
  CRPSDCPAR2 <- calculate_crps_from_quantiles(Y_test, DQuantAR2_OOS, QQ)
  
  

  ret <- list(coverageQRAR1, coverageQRAR2, coverageQRAR3, coverageCQRAR1, coverageCQRAR2, coverageCQRAR3, coverageDCPAR2, 
              CRPSQRAR1, CRPSQRAR2, CRPSQRAR3, CRPSCQRAR1, CRPSCQRAR2, CRPSCQRAR3, CRPSDCPAR2)
  
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
  library(quantreg)
  library(quantregForest)
  source("functions.R")
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
coveragetotQRAR3 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR2 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR3 <- matrix(NA,n3,length(QQ))
coveragetotDCPAR2 <- matrix(NA,n3,length(QQ))

CRPStotQRAR1 <- rep(NA,n3)
CRPStotQRAR2 <- rep(NA,n3)
CRPStotQRAR3 <- rep(NA,n3)

CRPStotCQRAR1 <- rep(NA,n3)
CRPStotCQRAR2 <- rep(NA,n3)
CRPStotCQRAR3 <- rep(NA,n3)

CRPStotDCPAR2 <- rep(NA,n3)




index = 1
for(res in results){
  coveragetotQRAR1[index,] <- res[[1]]
  coveragetotQRAR2[index,] <- res[[2]]
  coveragetotQRAR3[index,] <- res[[3]]
  
  coveragetotCQRAR1[index,] <- res[[4]]
  coveragetotCQRAR2[index,] <- res[[5]]
  coveragetotCQRAR3[index,] <- res[[6]]
  
  coveragetotDCPAR2[index,] <- res[[7]]
  
  CRPStotQRAR1[index] <- res[[8]]
  CRPStotQRAR2[index] <- res[[9]]
  CRPStotQRAR3[index] <- res[[10]]
  
  CRPStotCQRAR1[index] <- res[[11]]
  CRPStotCQRAR2[index] <- res[[12]]
  CRPStotCQRAR3[index] <- res[[13]]
  
  CRPStotDCPAR2[index] <- res[[14]]
  
  index <- index + 1
}



# RESULTS EVALUATION

# Calculate average coverage across all iterations and for each quantile
average_coverageQRAR1 <- compute_average_coverage(coveragetotQRAR1, QQ) #vector of length = length(QQ)
average_coverageQRAR2 <- compute_average_coverage(coveragetotQRAR2, QQ)
average_coverageQRAR3 <- compute_average_coverage(coveragetotQRAR3, QQ)

average_coverageCQRAR1 <- compute_average_coverage(coveragetotCQRAR1, QQ)
average_coverageCQRAR2 <- compute_average_coverage(coveragetotCQRAR2, QQ)
average_coverageCQRAR3 <- compute_average_coverage(coveragetotCQRAR3, QQ)

average_coverageDCPAR2 <- compute_average_coverage(coveragetotDCPAR2, QQ)



#----------------------- QR AR(1)
{
  resultsPitSTQRAR1 <- compute_results(average_coverageQRAR1, n2*n3, QQ)
  
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
  print(paste("CRPS of QR:", mean(CRPStotQRAR1)))
}


#----------------------- QR AR(2)
{
  resultsPitSTQRAR2 <- compute_results(average_coverageQRAR2, n2*n3, QQ)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTQRAR2$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTQRAR2$EmpiricalCoverage > resultsPitSTQRAR2$Quantile & !resultsPitSTQRAR2$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTQRAR2$EmpiricalCoverage < resultsPitSTQRAR2$Quantile & !resultsPitSTQRAR2$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of QR:", mean(abs(resultsPitSTQRAR2$Quantile-resultsPitSTQRAR2$EmpiricalCoverage))))
  print(paste("CRPS of QR:", mean(CRPStotQRAR2)))
  
}


#----------------------- QR AR(3)
{
  resultsPitSTQRAR3 <- compute_results(average_coverageQRAR3, n2*n3, QQ)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTQRAR3$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTQRAR3$EmpiricalCoverage > resultsPitSTQRAR3$Quantile & !resultsPitSTQRAR3$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTQRAR3$EmpiricalCoverage < resultsPitSTQRAR3$Quantile & !resultsPitSTQRAR3$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of QR:", mean(abs(resultsPitSTQRAR3$Quantile-resultsPitSTQRAR3$EmpiricalCoverage))))
  print(paste("CRPS of QR:", mean(CRPStotQRAR3)))
  
}


#----------------------- CQR AR(1)
{
  resultsPitSTCQRAR1 <- compute_results(average_coverageCQRAR1, n2*n3, QQ)
  
  
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
  print(paste("CRPS of CQR:", mean(CRPStotCQRAR1)))
  
}


#----------------------- CQR AR(2)
{
  resultsPitSTCQRAR2 <- compute_results(average_coverageCQRAR2, n2*n3, QQ)
  
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTCQRAR2$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTCQRAR2$EmpiricalCoverage > resultsPitSTCQRAR2$Quantile & !resultsPitSTCQRAR2$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTCQRAR2$EmpiricalCoverage < resultsPitSTCQRAR2$Quantile & !resultsPitSTCQRAR2$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of CQR:", mean(abs(resultsPitSTCQRAR2$Quantile-resultsPitSTCQRAR2$EmpiricalCoverage))))
  print(paste("CRPS of CQR:", mean(CRPStotCQRAR2)))
  
}


#----------------------- CQR AR(3)
{
  resultsPitSTCQRAR3 <- compute_results(average_coverageCQRAR3, n2*n3, QQ)
  
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTCQRAR3$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTCQRAR3$EmpiricalCoverage > resultsPitSTCQRAR3$Quantile & !resultsPitSTCQRAR3$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTCQRAR3$EmpiricalCoverage < resultsPitSTCQRAR3$Quantile & !resultsPitSTCQRAR3$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of CQR:", mean(abs(resultsPitSTCQRAR3$Quantile-resultsPitSTCQRAR3$EmpiricalCoverage))))
  print(paste("CRPS of CQR:", mean(CRPStotCQRAR3)))
  
}



#----------------------- DCP AR(2)

{
  resultsPitSTDCPAR2 <- compute_results(average_coverageDCPAR2, n2*n3, QQ)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTDCPAR2$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTDCPAR2$EmpiricalCoverage > resultsPitSTDCPAR2$Quantile & !resultsPitSTDCPAR2$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTDCPAR2$EmpiricalCoverage < resultsPitSTDCPAR2$Quantile & !resultsPitSTDCPAR2$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for DCP:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of DCP:", mean(abs(resultsPitSTDCPAR2$Quantile-resultsPitSTDCPAR2$EmpiricalCoverage))))
  print(paste("CRPS of DCP:", mean(CRPStotDCPAR2)))
  
}









#------------- PLOT CALIBRATION CURVES OF QR and CQR with AR(1) MODEL

{
  x1 <- resultsPitSTCQRAR1$Quantile
  y1 <- resultsPitSTCQRAR1$EmpiricalCoverage
  x2 <- resultsPitSTQRAR1$Quantile
  y2 <- resultsPitSTQRAR1$EmpiricalCoverage
  
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
  legend("bottomright", legend = c("CQR AR1", "QR AR1"), col = c("blue", "red"), pch = c(16, 17), bty = "n", cex = 1.5, pt.cex = 1.5)
}


#------------- PLOT CALIBRATION CURVES OF QR and CQR and DCP with AR(2) MODEL

{
  x1 <- resultsPitSTCQRAR2$Quantile
  y1 <- resultsPitSTCQRAR2$EmpiricalCoverage
  x2 <- resultsPitSTQRAR2$Quantile
  y2 <- resultsPitSTQRAR2$EmpiricalCoverage
  x3 <- resultsPitSTDCPAR2$Quantile
  y3 <- resultsPitSTDCPAR2$EmpiricalCoverage
  
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
  legend("bottomright", legend = c("CQR AR2", "QR AR2", "DCP AR2"), col = c("blue", "red", "green"), pch = c(16, 17, 18), bty = "n", cex = 1.5, pt.cex = 1.5)
}



#------------- PLOT CALIBRATION CURVES OF QR and CQR with AR(3) MODEL
{
x1 <- resultsPitSTCQRAR3$Quantile
y1 <- resultsPitSTCQRAR3$EmpiricalCoverage
x2 <- resultsPitSTQRAR3$Quantile
y2 <- resultsPitSTQRAR3$EmpiricalCoverage

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
legend("bottomright", legend = c("CQR AR3", "QR AR3"), col = c("blue", "red"), pch = c(16, 17), bty = "n", cex = 1.5, pt.cex = 1.5)
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





