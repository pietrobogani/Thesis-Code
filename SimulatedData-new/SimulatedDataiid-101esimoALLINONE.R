# qui ora testo generando n2 volte il punto n+1. Nell' altro file, generavo tutti i punti tra n e n2 e testavo su loro, non solo 
#su n+1

#Il risultato è che nonostante i dati siano exchangeable, non ho copertura

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

source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")

file_path <- "ALLINONE_iid_101esimo.xlsx"
# Create a new workbook
wb <- createWorkbook()

# Add a worksheet named "Data"
addWorksheet(wb, "Data")

# Save the workbook (this creates the file)
saveWorkbook(wb, file_path, overwrite = TRUE)

wb <- loadWorkbook(file_path)


# da qui comincia la simulazione

run_simulation <- function(n,ratio_p_n){
  
  # Set parameters
  n2 <- 100      # Number of data points for testing
  p <- ratio_p_n * n
  num_exog_vars <-  p     # Parameter to control the number of exogenous variables
  beta <- runif(num_exog_vars, 0.1, 1) # Randomly generated coefficients for exogenous variables
  means <- rnorm(p)
  sigma <- diag(runif(p))
  
  # Initialize variables
  Y <- numeric(n)  
  exog_vars <- rmvnorm( n = n + 1, mean = means, sigma = sigma)

  # Simulate i.i.d data with exogenous variables 
  for (t in 1:n) {
    Y[t] <- sum(beta * exog_vars[t, ]) + rt(n = 1, df = 2)
  }
  

  # Prepare dataset for model 
  data_ar2 <- data.frame(I = 1, Y = Y)
  for (j in 1:num_exog_vars) {
    data_ar2[paste0("X", j)] <- exog_vars [1:n, j]
  }
  

  
  Y_test <- sum(beta * exog_vars[n+1, ]) + rt(n = n2, df = 2)
  

#--------------- QUANTILE REGRESSION 
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)

  Quant_OOS <- matrix(NA, length(QQ))
  #qrf_model <- quantregForest(x = as.matrix(data_ar2[,-c(1,2)]), y = as.matrix(data_ar2[,2]))
  

    for (jq in 1:length(QQ)) {  
       
       QR <- rq(make_formula(0,num_exog_vars), data = data_ar2, tau = QQ[jq])     
       Quant_OOS[jq] <- c(1,exog_vars[n+1, ]) %*% coef(QR)
      
       #QuantAR_OOS[,jq] <- predict(qrf_model, newdata = as.matrix(data_ar2_test[,-c(1,2)]), what = QQ[jq])
       
      
    }

  
  
#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  CQuant_OOS <- matrix(NA, length(QQ))
  Q_low <- matrix(NA, nrow(data_ar2), length(QQ))
  Q_high <- matrix(NA, nrow(data_ar2), length(QQ))
  full_length <- nrow(data_ar2)
  test_length = floor(full_length*50/100)
  data_ar2_1 <- data_ar2[1:test_length,] #train
  data_ar2_2 <- data_ar2[(test_length + 1) : full_length,] #calibration

  x0 <- data_ar2_1[-2]
  y0 <- data_ar2_1[,2]
  x1 <- data_ar2_2[,-2]
  y1 <- data_ar2_2[,2]
  x_test <- c(1,exog_vars[n+1, ])
  
  CQuant_OOS <- qCQR(make_formula(0, num_exog_vars),x0,y0,x1,y1,x_test,QQ) 



  
  #-------------------------------------------------------------------------------------------
  
  #CALIBRATION OF QR AND CQR
  

  coverageQR <- compute_coverage(Y_test, Quant_OOS, QQ)
  coverageCQR <- compute_coverage(Y_test, CQuant_OOS, QQ)
  
    
  ret <- list(coverageQR, coverageCQR)
  
  return(ret)
}

vector_n <- c(101,201,1001)
vector_p_n <- c(0.1,0.2,0.3,0.4)

n2 = 100
n3 = 100 
n_simul <- n3
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(0.01, QQ)




count <- 2 #so I leave first row empty

for (n in vector_n){
  for(p_n in vector_p_n){
    
    
# Setup parallel cluster
cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
clusterExport(cl, varlist=c("run_simulation","n","p_n")) # Export the simulation function to each cluster node
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
  run_simulation(n,p_n)
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
average_coverageQR <- compute_average_coverage(coveragetotQR, QQ) #vector of length = length(QQ)
average_coverageCQR <- compute_average_coverage(coveragetotCQR, QQ) #vector of length = length(QQ)



#----------------------- QR 
{
  resultsPitSTQR <- compute_results(average_coverageQR, n2*n3, QQ)
  
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
  resultsPitSTCQR <- compute_results(average_coverageCQR, n2*n3, QQ)
  
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



#------------ Write Results in the excel file


wb <- loadWorkbook(file_path)

# Access the worksheet 
sheet <- "Data"

#In the first column, write the parameters of the simulation:
writeData(wb, sheet = sheet, x = paste("n1 = ", n,", p/n = ", p_n, ", QR" ), startCol = 1, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n,", p/n = ", p_n, ", CQR" ), startCol = 1, startRow = count+1, colNames = FALSE)


#In the 2°,3°,4° columns, put how many inside, below and above CI
writeData(wb, sheet = sheet, x = sum(resultsPitSTQR$IsWithinCI) , startCol = 2, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQR$EmpiricalCoverage > resultsPitSTQR$Quantile & !resultsPitSTQR$IsWithinCI) , startCol = 3, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQR$EmpiricalCoverage < resultsPitSTQR$Quantile & !resultsPitSTQR$IsWithinCI) , startCol = 4, startRow = count, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQR$IsWithinCI) , startCol = 2, startRow = count + 1, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQR$EmpiricalCoverage > resultsPitSTCQR$Quantile & !resultsPitSTCQR$IsWithinCI) , startCol = 3, startRow = count + 1, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQR$EmpiricalCoverage < resultsPitSTCQR$Quantile & !resultsPitSTCQR$IsWithinCI) , startCol = 4, startRow = count + 1, colNames = FALSE)

#In the 5° column, the MAE will be placed
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQR$Quantile-resultsPitSTQR$EmpiricalCoverage)), startCol = 5, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQR$Quantile-resultsPitSTCQR$EmpiricalCoverage)), startCol = 5, startRow = count+1, colNames = FALSE)

#In the 6°,7° column, the MAE+ and MAE- will be placed

filtered_results <- resultsPitSTQR[resultsPitSTQR$EmpiricalCoverage > resultsPitSTQR$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count, colNames = FALSE)
filtered_results <- resultsPitSTCQR[resultsPitSTCQR$EmpiricalCoverage > resultsPitSTCQR$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count+1, colNames = FALSE)

filtered_results <- resultsPitSTQR[resultsPitSTQR$EmpiricalCoverage < resultsPitSTQR$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count, colNames = FALSE)
filtered_results <- resultsPitSTCQR[resultsPitSTCQR$EmpiricalCoverage < resultsPitSTCQR$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count+1, colNames = FALSE)


saveWorkbook(wb, file_path, overwrite = TRUE)


count <- count + 3 #lascio una riga vuota


  }
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





