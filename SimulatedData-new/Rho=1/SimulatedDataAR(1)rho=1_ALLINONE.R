
#Code updated to have 100 steps ahead in the test set and to produce also plots
#random forest missing for now


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
library(openxlsx)


source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")


file_path <- "ALLINONE_rho=1_Results.xlsx"
# Create a new workbook
wb <- createWorkbook()

# Add a worksheet named "Data"
addWorksheet(wb, "Data")

# Save the workbook (this creates the file)
saveWorkbook(wb, file_path, overwrite = TRUE)

wb <- loadWorkbook(file_path)




run_simulation <- function(n, phi_ar2){
  
  #n <- 1001 #3 elements will be lost in the DGP
  n2 <- 100
  
#------------------- Generate n + n2 points for the AR(1) model, this is our DGP
  
  Y <- numeric(n+n2)
  Y[1] <- rnorm(1)
  for (i in 2:(n+n2)) {
    Y[i] <- phi_ar2 * Y[i-1] + rt(n = 1, df = 2) #heavy tails
  }
  # Create lagged variables
  Y_lag1 <- c(NA, Y[-length(Y)])

  # Prepare dataset for model excluding the first three NA values
  indices <- c((n+1):(n+n2))
  data <- data.frame(I = 1, Y = Y[-c(1:3,indices)], Y_lag1 = Y_lag1[-c(1:3,indices)])

  
  data_test <- data.frame(I = 1, Y = Y[indices], Y_lag1 = Y_lag1[indices])

  Y_test <- Y[c((n+1):(n+n2))]
  
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)

  
#--------------- QUANTILE REGRESSION 
  
  QuantAR1_OOS <- matrix(0, n2, length(QQ))
  
    for (jq in 1:length(QQ)) {  
      
        QR1 <- rq(Y ~ Y_lag1, data = data, tau=QQ[jq])
        QuantAR1_OOS[,jq] <- as.matrix(data_test[,-2]) %*% coef(QR1)
      
    }
  
#--------------- QUANTILE RANDOM FOREST 
  
    # qrf_model <- quantregForest(x = data_ar2[, 3], y = data_ar2[,2])
    # QuantAR1_OOS[jq] <- predict(qrf_model, newdata = as.matrix(rep(Y_ar2[n],n2)), what = QQ[jq])[1]
  
#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  CQuantAR1_OOS <- matrix(NA, n2, length(QQ))
  full_length <- nrow(data)
  test_length = floor(full_length*50/100)
  data_1 <- data[1:test_length,] #train
  data_2 <- data[(test_length + 1) : full_length,] #calibration

  x0 <- data_1[,-2]
  y0 <- data_1[,2]
  x1 <- data_2[,-2]
  y1 <- data_2[,2]
  x_test <- data_test[,-2]
  formula <- make_formula(1, 0)
  
  CQuantAR1_OOS <- qCQR(formula,x0,y0,x1,y1,x_test,QQ)
  
  #--------------- CONFORMAL QUANTILE RANDOM FOREST
  
  #qrf_model <- quantregForest(x = data_1[, 3], y = data_1[,2])
  

  
  #---------------- CALIBRATION OF QR AND CQR
  
  coverageQRAR1 <- compute_coverage(Y_test, QuantAR1_OOS, QQ)
  coverageCQRAR1 <- compute_coverage(Y_test, CQuantAR1_OOS, QQ)
  
  
    
  ret <- list(coverageQRAR1, coverageCQRAR1)
  
  return(ret)
}




vector_n <- c(101,201,1001)
vector_phi <- c(0.95, 1, 1.05)

n2 = 100
n3 = 100 
n_simul <- n3
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(0.01, QQ)

count <- 2 #so I leave first row empty

for (n in vector_n){
  for(phi in vector_phi){



# Setup parallel cluster
cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
clusterExport(cl, varlist=c("run_simulation","n","phi")) # Export the simulation function to each cluster node
clusterEvalQ(cl, { 
  library(readxl)
  library(quantreg)
  library(quantregForest)
  source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")
}
) # Load required libraries in each cluster node, repeat as necessary for other libraries

# Run simulations in parallel
results <- parLapply(cl, seeds, function(seed) {
  set.seed(seed)
  run_simulation(n, phi)
})


# Stop the cluster
stopCluster(cl)

# Extract results
coveragetotQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR1 <- matrix(NA,n3,length(QQ))




index = 1
for(res in results){
  coveragetotQRAR1[index,] <- res[[1]]
  coveragetotCQRAR1[index,] <- res[[2]]

  index <- index + 1
}

average_coverageQRAR1 <- compute_average_coverage(coveragetotQRAR1, QQ) #vector of length = length(QQ)
average_coverageCQRAR1 <- compute_average_coverage(coveragetotCQRAR1, QQ) #vector of length = length(QQ)

resultsPitSTQRAR1 <- compute_results(average_coverageQRAR1, n2*n3, QQ)
resultsPitSTCQRAR1 <- compute_results(average_coverageCQRAR1, n2*n3, QQ)

#------------ Write Results in the excel file


wb <- loadWorkbook(file_path)

# Access the worksheet 
sheet <- "Data"

#In the first column, write the parameters of the simulation:
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", phi = ", phi, ", QR AR(1)" ), startCol = 1, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", phi = ", phi, ", CQR AR(1)" ), startCol = 1, startRow = count+1, colNames = FALSE)

#In the 2°,3°,4° columns, put how many inside, below and above CI
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR1$IsWithinCI) , startCol = 2, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR1$EmpiricalCoverage > resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI) , startCol = 3, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR1$EmpiricalCoverage < resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI) , startCol = 4, startRow = count, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR1$IsWithinCI) , startCol = 2, startRow = count + 1, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR1$EmpiricalCoverage > resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI) , startCol = 3, startRow = count + 1, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR1$EmpiricalCoverage < resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI) , startCol = 4, startRow = count + 1, colNames = FALSE)

#In the 5° column, the MAE will be placed
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR1$Quantile-resultsPitSTQRAR1$EmpiricalCoverage)), startCol = 5, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR1$Quantile-resultsPitSTCQRAR1$EmpiricalCoverage)), startCol = 5, startRow = count+1, colNames = FALSE)

#In the 6°,7° column, the MAE+ and MAE- will be placed
filtered_results <- resultsPitSTQRAR1[resultsPitSTQRAR1$EmpiricalCoverage > resultsPitSTQRAR1$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count, colNames = FALSE)
filtered_results <- resultsPitSTCQRAR1[resultsPitSTCQRAR1$EmpiricalCoverage > resultsPitSTCQRAR1$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count+1, colNames = FALSE)


filtered_results <- resultsPitSTQRAR1[resultsPitSTQRAR1$EmpiricalCoverage < resultsPitSTQRAR1$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count, colNames = FALSE)
filtered_results <- resultsPitSTCQRAR1[resultsPitSTCQRAR1$EmpiricalCoverage < resultsPitSTCQRAR1$Quantile, ]
writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count+1, colNames = FALSE)

saveWorkbook(wb, file_path, overwrite = TRUE)


count <- count + 3 #lascio una riga vuota



#------------- PLOT CALIBRATION CURVES OF QR and CQR and DCP with AR(1) MODEL

{
  x1 <- resultsPitSTCQRAR1$Quantile
  y1 <- resultsPitSTCQRAR1$EmpiricalCoverage
  x2 <- resultsPitSTQRAR1$Quantile
  y2 <- resultsPitSTQRAR1$EmpiricalCoverage
  
  
  # Create the plot
  plot(x1, y1, type = "n", xlab = "Quantile", ylab = "Empirical Coverage", main = paste("n1 = ", n-3,", phi = ", phi))
  
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
  
  
  
  # Add the diagonal line
  abline(a = 0, b = 1, col = "black", lty = 2)
  
  # Optional: Add a legend
  legend("bottomright", legend = c("CQR AR1", "QR AR1"), col = c("blue", "red"), pch = c(16, 17, 18), bty = "n", cex = 1.5, pt.cex = 1.5)
}




  }
}
