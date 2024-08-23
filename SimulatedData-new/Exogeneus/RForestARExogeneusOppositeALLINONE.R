#sia con errore normale sia con alfa a controllare errore, ar2 non performa meglio di ar1 e ar3

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
library(openxlsx)
source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")

file_path <- "ALLINONE_Exogenous_Results_Opposite_RForest.xlsx"
#file_path <- "ALLINONE_Exogenous_Results_Opposite_RForesttemp.xlsx"


# Check if the file exists
if (!file.exists(file_path)) {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add a worksheet named "Data"
  addWorksheet(wb, "Data")
  
  # Save the workbook (this creates the file)
  saveWorkbook(wb, file_path, overwrite = TRUE)
} else {
  # Load the existing workbook
  wb <- loadWorkbook(file_path)
}




run_simulation <- function(n,ratio_p_n,alpha){
  
  # Set parameters
  #n <- 31      # Number of data points for training. 3 will always be lost!
  n2 <- 100    # Number of data points for testing
  p <- ratio_p_n * n
  phi_ar2 <- c(0.5, -0.2)   # AR(2) coefficients
  num_exog_vars <-  p     # Parameter to control the number of exogenous variables
  beta <- runif(num_exog_vars, 0, 1) # Randomly generated coefficients for exogenous variables
  means <- rnorm(p)
  sigma <- diag(runif(p,0,10))
  
  # Initialize variables
  Y_ar2 <- numeric(n + n2)  # Series for AR(2) data
  exog_vars <- matrix(nrow = n + n2, ncol = num_exog_vars)  # Matrix for exogenous variables
  
  exog_vars <- rmvnorm( n = n + n2, mean = means, sigma = sigma)

  # Initialize Y values
  Y_ar2[1:2] <- rnorm(2)  
  
  # Simulate AR(2) with exogenous variables for training data
  for (t in 3:(n + n2)) {
    Y_ar2[t] <- phi_ar2[1] * Y_ar2[t - 1] + phi_ar2[2] * Y_ar2[t - 2] +
      sum(beta * exog_vars[t, ]) + alpha * rt(1, df = 2)
  }
  
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])
  Y_lag3_ar2 <- c(NA, NA, NA, Y_ar2[1:(length(Y_ar2)-3)])
  
  # Prepare dataset for model excluding the first three NA values
  indices <- c((n+1):(n+n2))
  data <- data.frame(I = 1, Y = Y_ar2[-c(1:3,indices)], Y_lag1 = Y_lag1_ar2[-c(1:3,indices)], Y_lag2 = Y_lag2_ar2[-c(1:3,indices)], Y_lag3 = Y_lag3_ar2[-c(1:3,indices)])
  for (j in 1:num_exog_vars) {
    data[paste0("X", j)] <- exog_vars[-c(1:3,indices), j]
  }
  
  
  data_test <- data.frame(I = 1, Y = Y_ar2[indices], Y_lag1 = Y_lag1_ar2[indices], Y_lag2 = Y_lag2_ar2[indices], Y_lag3 = Y_lag3_ar2[indices])
  for (j in 1:num_exog_vars) {
    data_test[paste0("X", j)] <- exog_vars[indices, j]
  }
  
  Y_test <- Y_ar2[c((n+1):(n+n2))]
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)
  
  
  #--------------- QUANTILE REGRESSION 
  
  QuantAR0_OOS <- matrix(0, n2, length(QQ))
  QuantAR1_OOS <- matrix(0, n2, length(QQ))
  QuantAR2_OOS <- matrix(0, n2, length(QQ))
  QuantAR3_OOS <- matrix(0, n2, length(QQ))
# 
 
  #--------------- QUANTILE RANDOM FOREST 
  
   {
    qrf_model1 <- quantregForest(x = (data[,-c(1,2,4,5)]), y = (data[,2]))
    qrf_model2 <- quantregForest(x = (data[,-c(1,2,5)]), y = (data[,2]))
    qrf_model3 <- quantregForest(x = (data[,-c(1,2)]), y = (data[,2]))
   
   for (jq in 1:length(QQ)) {  
     
     QuantAR1_OOS[,jq] <- predict(qrf_model1, newdata = as.matrix(data_test[,-c(1,2,4,5)]), what = (1-QQ[jq]))
     QuantAR2_OOS[,jq] <- predict(qrf_model2, newdata = as.matrix(data_test[,-c(1,2,5)]), what = (1-QQ[jq]))
     QuantAR3_OOS[,jq] <- predict(qrf_model3, newdata = as.matrix(data_test[,-c(1,2)]), what = (1-QQ[jq]))
     
   }
   }
  
  #--------------- CONFORMAL QUANTILE REGRESSION
  
  full_length <- nrow(data)
  test_length = floor(full_length*50/100)
  data_1 <- data[1:test_length,] #train
  data_2 <- data[(test_length + 1) : full_length,] #calibration
  
  CQuantAR0_OOS <- matrix(0, n2, length(QQ))
  CQuantAR1_OOS <- matrix(0, n2, length(QQ))
  CQuantAR2_OOS <- matrix(0, n2, length(QQ))
  CQuantAR3_OOS <- matrix(0, n2, length(QQ))

  
  #--------------- CONFORMAL QUANTILE RANDOM FOREST
  
   {
   #  
     x0 <- data_1[,-c(1,2,4,5)]
     y0 <- data_1[,2]
     x1 <- data_2[,-c(1,2,4,5)]
     y1 <- data_2[,2]
     x_test <- data_test[,-c(1,2,4,5)]
     
     CQuantAR1_OOS <- qCQRF_opposite(x0,y0,x1,y1,x_test,QQ) 
   #  
       x0 <- data_1[,-c(1,2,5)]
       y0 <- data_1[,2]
       x1 <- data_2[,-c(1,2,5)]
       y1 <- data_2[,2]
       x_test <- data_test[,-c(1,2,5)]
       
       CQuantAR2_OOS <- qCQRF_opposite(x0,y0,x1,y1,x_test,QQ) 
   #  
     x0 <- data_1[,-c(1,2)]
     y0 <- data_1[,2]
     x1 <- data_2[,-c(1,2)]
     y1 <- data_2[,2]
     x_test <- data_test[,-c(1,2)]
     
     CQuantAR3_OOS <- qCQRF_opposite(x0,y0,x1,y1,x_test,QQ) 
    }


  #-------------------------------------------------------------------------------------------
  
  # CALIBRATION
  coverageQRAR0 <- compute_coverage(Y_test, QuantAR0_OOS, QQ)
  coverageQRAR1 <- compute_coverage(Y_test, QuantAR1_OOS, QQ)
  coverageQRAR2 <- compute_coverage(Y_test, QuantAR2_OOS, QQ)
  coverageQRAR3 <- compute_coverage(Y_test, QuantAR3_OOS, QQ)
  
  coverageCQRAR0 <- compute_coverage(Y_test, CQuantAR0_OOS, QQ)
  coverageCQRAR1 <- compute_coverage(Y_test, CQuantAR1_OOS, QQ)
  coverageCQRAR2 <- compute_coverage(Y_test, CQuantAR2_OOS, QQ)
  coverageCQRAR3 <- compute_coverage(Y_test, CQuantAR3_OOS, QQ)
  


  
  
  
  ret <- list(coverageQRAR1, coverageQRAR2, coverageQRAR3, coverageCQRAR1, coverageCQRAR2, coverageCQRAR3,
              coverageQRAR0, coverageCQRAR0)
  
  return(ret)
}



vector_n <- c(101,201,1001)
vector_p_n <- c(0.1,0.2,0.3,0.4)
vector_alpha <- c(1) 


MAEQRAR2 <- c()
MAECQRAR2 <- c()

n2 = 100
n3 = 100
n_simul <- n3
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(0.01, QQ)

count <- 2 #so I leave first row empty









for (n in vector_n){
  for(p_n in vector_p_n){
    for(alpha in vector_alpha){
      
# Setup parallel cluster
cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
clusterExport(cl, varlist=c("run_simulation","n","p_n","alpha")) # Export the simulation function to each cluster node
clusterEvalQ(cl, { 
  library(readxl)
  library(quantreg)
  library(quantregForest)
  library(mvtnorm)
  source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")
}
) # Load required libraries in each cluster node, repeat as necessary for other libraries

# Run simulations in parallel
results <- parLapply(cl, seeds, function(seed) {
  set.seed(seed)
  run_simulation(n,p_n,alpha)
})


# Stop the cluster
stopCluster(cl)

# Extract results
coveragetotQRAR0 <- matrix(NA,n3,length(QQ))
coveragetotQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotQRAR2 <- matrix(NA,n3,length(QQ))
coveragetotQRAR3 <- matrix(NA,n3,length(QQ))

coveragetotCQRAR0 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR2 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR3 <- matrix(NA,n3,length(QQ))






index = 1
for(res in results){
  coveragetotQRAR1[index,] <- res[[1]]
  coveragetotQRAR2[index,] <- res[[2]]
  coveragetotQRAR3[index,] <- res[[3]]
  
  coveragetotCQRAR1[index,] <- res[[4]]
  coveragetotCQRAR2[index,] <- res[[5]]
  coveragetotCQRAR3[index,] <- res[[6]]
  

  coveragetotQRAR0[index,] <- res[[7]]
  coveragetotCQRAR0[index,] <- res[[8]]
  
  
  index <- index + 1
}



# RESULTS EVALUATION


# Calculate average coverage across all iterations and for each quantile
average_coverageQRAR1 <- compute_average_coverage(coveragetotQRAR1, QQ) 
average_coverageQRAR2 <- compute_average_coverage(coveragetotQRAR2, QQ)
average_coverageQRAR3 <- compute_average_coverage(coveragetotQRAR3, QQ)

average_coverageCQRAR1 <- compute_average_coverage(coveragetotCQRAR1, QQ) #CQR should have on the right at least (1-alpha) = QQ probability 
average_coverageCQRAR2 <- compute_average_coverage(coveragetotCQRAR2, QQ)
average_coverageCQRAR3 <- compute_average_coverage(coveragetotCQRAR3, QQ)

average_coverageQRAR0 <- compute_average_coverage(coveragetotQRAR0, QQ) #vector of length = length(QQ)
average_coverageCQRAR0 <- compute_average_coverage(coveragetotCQRAR0, QQ) #vector of length = length(QQ)


resultsPitSTQRAR0 <- compute_results(average_coverageQRAR0, n2*n3, (1-QQ))
resultsPitSTQRAR1 <- compute_results(average_coverageQRAR1, n2*n3, (1-QQ))
resultsPitSTQRAR2 <- compute_results(average_coverageQRAR2, n2*n3, (1-QQ))
resultsPitSTQRAR3 <- compute_results(average_coverageQRAR3, n2*n3, (1-QQ))
resultsPitSTCQRAR0 <- compute_results(average_coverageCQRAR0, n2*n3, (1-QQ))
resultsPitSTCQRAR1 <- compute_results(average_coverageCQRAR1, n2*n3, (1-QQ)) #CQR should have its empirical coverage <= quantile target
resultsPitSTCQRAR2 <- compute_results(average_coverageCQRAR2, n2*n3, (1-QQ))
resultsPitSTCQRAR3 <- compute_results(average_coverageCQRAR3, n2*n3, (1-QQ))
  



#------------ Write Results in the excel file


wb <- loadWorkbook(file_path)

# Access the worksheet 
sheet <- "Data"

#In the first column, write the parameters of the simulation:
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", p/n = ", p_n,  ", alpha = ", alpha, ", QR AR(1)" ), startCol = 1, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", p/n = ", p_n,  ", alpha = ", alpha, ", CQR AR(1)" ), startCol = 1, startRow = count+1, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", p/n = ", p_n,  ", alpha = ", alpha, ", QR AR(2)" ), startCol = 1, startRow = count+2, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", p/n = ", p_n,  ", alpha = ", alpha, ", CQR AR(2)" ), startCol = 1, startRow = count+3, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", p/n = ", p_n,  ", alpha = ", alpha, ", QR AR(3)" ), startCol = 1, startRow = count+4, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", p/n = ", p_n,  ", alpha = ", alpha, ", CQR AR(3)" ), startCol = 1, startRow = count+5, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", p/n = ", p_n,  ", alpha = ", alpha, ", QR AR(0)" ), startCol = 1, startRow = count+6, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", p/n = ", p_n,  ", alpha = ", alpha, ", CQR AR(0)" ), startCol = 1, startRow = count+7, colNames = FALSE)

#In the 2°,3°,4° columns, put how many inside, below and above CI
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR1$IsWithinCI) , startCol = 2, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR1$EmpiricalCoverage > resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI) , startCol = 3, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR1$EmpiricalCoverage < resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI) , startCol = 4, startRow = count, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR1$IsWithinCI) , startCol = 2, startRow = count + 1, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR1$EmpiricalCoverage > resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI) , startCol = 3, startRow = count + 1, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR1$EmpiricalCoverage < resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI) , startCol = 4, startRow = count + 1, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR2$IsWithinCI) , startCol = 2, startRow = count + 2, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR2$EmpiricalCoverage > resultsPitSTQRAR2$Quantile & !resultsPitSTQRAR2$IsWithinCI) , startCol = 3, startRow = count + 2, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR2$EmpiricalCoverage < resultsPitSTQRAR2$Quantile & !resultsPitSTQRAR2$IsWithinCI) , startCol = 4, startRow = count + 2, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR2$IsWithinCI) , startCol = 2, startRow = count + 3, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR2$EmpiricalCoverage > resultsPitSTCQRAR2$Quantile & !resultsPitSTCQRAR2$IsWithinCI) , startCol = 3, startRow = count + 3, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR2$EmpiricalCoverage < resultsPitSTCQRAR2$Quantile & !resultsPitSTCQRAR2$IsWithinCI) , startCol = 4, startRow = count + 3, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR3$IsWithinCI) , startCol = 2, startRow = count + 4, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR3$EmpiricalCoverage > resultsPitSTQRAR3$Quantile & !resultsPitSTQRAR3$IsWithinCI) , startCol = 3, startRow = count + 4, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR3$EmpiricalCoverage < resultsPitSTQRAR3$Quantile & !resultsPitSTQRAR3$IsWithinCI) , startCol = 4, startRow = count + 4, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR3$IsWithinCI) , startCol = 2, startRow = count + 5, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR3$EmpiricalCoverage > resultsPitSTCQRAR3$Quantile & !resultsPitSTCQRAR3$IsWithinCI) , startCol = 3, startRow = count + 5, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR3$EmpiricalCoverage < resultsPitSTCQRAR3$Quantile & !resultsPitSTCQRAR3$IsWithinCI) , startCol = 4, startRow = count + 5, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR0$IsWithinCI) , startCol = 2, startRow = count + 6, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR0$EmpiricalCoverage > resultsPitSTCQRAR0$Quantile & !resultsPitSTCQRAR0$IsWithinCI) , startCol = 3, startRow = count + 6, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR0$EmpiricalCoverage < resultsPitSTCQRAR0$Quantile & !resultsPitSTCQRAR0$IsWithinCI) , startCol = 4, startRow = count + 6, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR0$IsWithinCI) , startCol = 2, startRow = count + 7, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR0$EmpiricalCoverage > resultsPitSTCQRAR0$Quantile & !resultsPitSTCQRAR0$IsWithinCI) , startCol = 3, startRow = count + 7, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR0$EmpiricalCoverage < resultsPitSTCQRAR0$Quantile & !resultsPitSTCQRAR0$IsWithinCI) , startCol = 4, startRow = count + 7, colNames = FALSE)

#In the 5° column, the MAE will be placed
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR1$Quantile-resultsPitSTQRAR1$EmpiricalCoverage)), startCol = 5, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR1$Quantile-resultsPitSTCQRAR1$EmpiricalCoverage)), startCol = 5, startRow = count+1, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR2$Quantile-resultsPitSTQRAR2$EmpiricalCoverage)), startCol = 5, startRow = count+2, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR2$Quantile-resultsPitSTCQRAR2$EmpiricalCoverage)), startCol = 5, startRow = count+3, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR3$Quantile-resultsPitSTQRAR3$EmpiricalCoverage)), startCol = 5, startRow = count+4, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR3$Quantile-resultsPitSTCQRAR3$EmpiricalCoverage)), startCol = 5, startRow = count+5, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR0$Quantile-resultsPitSTQRAR0$EmpiricalCoverage)), startCol = 5, startRow = count+6, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR0$Quantile-resultsPitSTCQRAR0$EmpiricalCoverage)), startCol = 5, startRow = count+7, colNames = FALSE)


#In the 6°,7° column, the MAE+ and MAE- will be placed

filtered_results <- resultsPitSTQRAR1[resultsPitSTQRAR1$EmpiricalCoverage > resultsPitSTQRAR1$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count, colNames = FALSE)
filtered_results <- resultsPitSTCQRAR1[resultsPitSTCQRAR1$EmpiricalCoverage > resultsPitSTCQRAR1$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count+1, colNames = FALSE)
filtered_results <- resultsPitSTQRAR2[resultsPitSTQRAR2$EmpiricalCoverage > resultsPitSTQRAR2$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count+2, colNames = FALSE)
filtered_results <- resultsPitSTCQRAR2[resultsPitSTCQRAR2$EmpiricalCoverage > resultsPitSTCQRAR2$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count+3, colNames = FALSE)
filtered_results <- resultsPitSTQRAR3[resultsPitSTQRAR3$EmpiricalCoverage > resultsPitSTQRAR3$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count+4, colNames = FALSE)
filtered_results <- resultsPitSTCQRAR3[resultsPitSTCQRAR3$EmpiricalCoverage > resultsPitSTCQRAR3$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count+5, colNames = FALSE)
filtered_results <- resultsPitSTQRAR0[resultsPitSTQRAR0$EmpiricalCoverage > resultsPitSTQRAR0$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count+6, colNames = FALSE)
filtered_results <- resultsPitSTCQRAR0[resultsPitSTCQRAR0$EmpiricalCoverage > resultsPitSTCQRAR0$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count+7, colNames = FALSE)

filtered_results <- resultsPitSTQRAR1[resultsPitSTQRAR1$EmpiricalCoverage < resultsPitSTQRAR1$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count, colNames = FALSE)
filtered_results <- resultsPitSTCQRAR1[resultsPitSTCQRAR1$EmpiricalCoverage < resultsPitSTCQRAR1$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count+1, colNames = FALSE)
filtered_results <- resultsPitSTQRAR2[resultsPitSTQRAR2$EmpiricalCoverage < resultsPitSTQRAR2$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count+2, colNames = FALSE)
filtered_results <- resultsPitSTCQRAR2[resultsPitSTCQRAR2$EmpiricalCoverage < resultsPitSTCQRAR2$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count+3, colNames = FALSE)
filtered_results <- resultsPitSTQRAR3[resultsPitSTQRAR3$EmpiricalCoverage < resultsPitSTQRAR3$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count+4, colNames = FALSE)
filtered_results <- resultsPitSTCQRAR3[resultsPitSTCQRAR3$EmpiricalCoverage < resultsPitSTCQRAR3$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count+5, colNames = FALSE)
filtered_results <- resultsPitSTQRAR0[resultsPitSTQRAR0$EmpiricalCoverage < resultsPitSTQRAR0$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count+6, colNames = FALSE)
filtered_results <- resultsPitSTCQRAR0[resultsPitSTCQRAR0$EmpiricalCoverage < resultsPitSTCQRAR0$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count+7, colNames = FALSE)


saveWorkbook(wb, file_path, overwrite = TRUE)


count <- count + 9 #lascio una riga vuota


#------------- PLOT CALIBRATION CURVES OF QR and CQR for AR(2) MODEL
# Create data frames for both datasets
{
df1 <- data.frame(
  Quantile = rev(resultsPitSTCQRAR2$Quantile),
  EmpiricalCoverage = rev(resultsPitSTCQRAR2$EmpiricalCoverage),
  Group = "CQR QRF"
)

df2 <- data.frame(
  Quantile = rev(resultsPitSTQRAR2$Quantile),
  EmpiricalCoverage = rev(resultsPitSTQRAR2$EmpiricalCoverage),
  Group = "QRF"
)

# Combine the data frames
df <- bind_rows(df1, df2)

# Define the colors
cqr_colors <- "#0000FF"
qr_colors <- "#FF0000"

# Define the legend position based on the value of n - 3
legend_position <- if ((n - 3 == 98) || (n - 3 == 998 && p_n == 0.1)) {
  c(0.85, 0.15)
} else {
  "none"
}

p <- ggplot(df, aes(x = Quantile, y = EmpiricalCoverage, color = Group)) +
  geom_step(direction = "vh", linewidth = 1) + # Use `linewidth` instead of `size` for lines
  geom_point(size = 3) + # Add points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) + # Use `linewidth` for the diagonal line
  scale_color_manual(values = c("CQR QRF" = cqr_colors, "QRF" = qr_colors)) + # Manual color scale
  labs(title = paste("n1 =", n - 3, "p/n =", p_n), x = "Quantile Levels", y = "Empirical Coverage") + # Add labels and title
  theme_minimal() +
  theme(
    legend.position = legend_position, # Set legend position based on condition
    legend.title = element_blank(),
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)
  )

# Print the plot
print(p)
# Save the plot as a PDF file
ggsave(filename = paste0("Exogeneus_calibration_n", n - 3,"_", p_n, "_AR(2).pdf"), plot = p, width = 7, height = 5)
}

MAEQRAR2 <- c(MAEQRAR2,mean(abs(resultsPitSTQRAR2$Quantile-resultsPitSTQRAR2$EmpiricalCoverage)))
MAECQRAR2 <- c(MAECQRAR2,mean(abs(resultsPitSTCQRAR2$Quantile-resultsPitSTCQRAR2$EmpiricalCoverage)))

  }
  }
}






# Calculate the difference between the two vectors
difference <- MAEQRAR2 - MAECQRAR2

# Create a data frame
df <- data.frame(
  x = vector_n,
  MAEQRAR2 = MAEQRAR2,
  MAECQRAR2 = MAECQRAR2,
  Difference = difference
)

ggplot(df) +
  geom_point(aes(x = x, y = MAEQRAR2), color = "blue") +
  geom_line(aes(x = x, y = MAEQRAR2), color = "blue", linetype = "dashed") +
  geom_point(aes(x = x, y = MAECQRAR2), color = "red") +
  geom_line(aes(x = x, y = MAECQRAR2), color = "red", linetype = "solid") +
  geom_point(aes(x = x, y = Difference), color = "green") +
  geom_line(aes(x = x, y = Difference), color = "green", linetype = "dotted") +
  scale_x_reverse() +
  labs(title = "Plot of MAEQR, MAECQR, and Their Difference",
       x = "Numbers of observations",
       y = "MAE") +
  theme_minimal()
