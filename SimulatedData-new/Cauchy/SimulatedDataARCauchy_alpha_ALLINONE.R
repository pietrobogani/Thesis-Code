#also random forest and alpha


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
library(ggplot2)


source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")


file_path <- "ALLINONE_Cauchy_alpha_Results.xlsx" #for quantile regression
#file_path <- "ALLINONE_Cauchy_ResultsRForest.xlsx"

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




run_simulation <- function(n,alpha){
  
  #n <- 1001 #3 elements will be lost in the DGP
  n2 <- 100
  
#------------------- Generate n + n2 points for the AR(1) model, this is our DGP
  
  phi_ar2 <- c(0.5, -0.2)  # AR coefficients for AR(2)
  
  Y_ar2 <- numeric(n+n2)
  Y_ar2[1:2] <- rnorm(2)
  for (i in 3:(n+n2)) {
    Y_ar2[i] <- phi_ar2[1] * Y_ar2[i-1] + phi_ar2[2] * Y_ar2[i-2] + alpha*rt(n = 1, df = 1) #CAUCHY ERRORS!
  }
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])
  Y_lag3_ar2 <- c(NA, NA, NA, Y_ar2[1:(length(Y_ar2)-3)])

  # Prepare dataset for model excluding the first three NA values
  indices <- c((n+1):(n+n2))
  data <- data.frame(I = 1, Y = Y_ar2[-c(1:3,indices)], Y_lag1 = Y_lag1_ar2[-c(1:3,indices)], Y_lag2 = Y_lag2_ar2[-c(1:3,indices)], Y_lag3 = Y_lag3_ar2[-c(1:3,indices)])

  
  data_test <- data.frame(I = 1, Y = Y_ar2[indices], Y_lag1 = Y_lag1_ar2[indices], Y_lag2 = Y_lag2_ar2[indices], Y_lag3 = Y_lag3_ar2[indices])

  Y_test <- Y_ar2[c((n+1):(n+n2))]
  
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)

  
#--------------- QUANTILE REGRESSION 
  
  QuantAR1_OOS <- matrix(0, n2, length(QQ))
  QuantAR2_OOS <- matrix(0, n2, length(QQ))
  QuantAR3_OOS <- matrix(0, n2, length(QQ))
  
      for (jq in 1:length(QQ)) {  
        
          QR1 <- rq(Y ~ Y_lag1, data = data, tau=QQ[jq])
          QuantAR1_OOS[,jq] <- as.matrix(data_test[,-c(2,4,5)]) %*% coef(QR1)
          
          QR2 <- rq(Y ~ Y_lag1 + Y_lag2 , data = data, tau=QQ[jq])
          QuantAR2_OOS[,jq] <- as.matrix(data_test[,-c(2,5)]) %*% coef(QR2)
          
          QR3 <- rq(Y ~ Y_lag1 + Y_lag2 + Y_lag3, data = data, tau=QQ[jq])
          QuantAR3_OOS[,jq] <- as.matrix(data_test[,-2]) %*% coef(QR3)
        
      }
      
#--------------- QUANTILE RANDOM FOREST 
  
 #  qrf_model1 <- quantregForest(x = as.matrix(data[, 3]), y = data[,2])
 #  qrf_model2 <- quantregForest(x = as.matrix(data[, c(3,4)]), y = data[,2])
 #  qrf_model3 <- quantregForest(x = as.matrix(data[, c(3,4,5)]), y = data[,2])
 #  
 #    for (jq in 1:length(QQ)) {  
 #   
 #      QuantAR1_OOS[,jq] <- predict(qrf_model1, newdata = as.matrix(data_test[,3]), what = QQ[jq])
 #      QuantAR2_OOS[,jq] <- predict(qrf_model2, newdata = as.matrix(data_test[,c(3,4)]), what = QQ[jq])      
 #      QuantAR3_OOS[,jq] <- predict(qrf_model3, newdata = as.matrix(data_test[,c(3,4,5)]), what = QQ[jq])
 #      
 # }
#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  full_length <- nrow(data)
  test_length = floor(full_length*50/100)
  data_1 <- data[1:test_length,] #train
  data_2 <- data[(test_length + 1) : full_length,] #calibration

  
   #AR1
   x0 <- data_1[,-c(2,4,5)]
   y0 <- data_1[,2]
   x1 <- data_2[,-c(2,4,5)]
   y1 <- data_2[,2]
   x_test <- data_test[,-c(2,4,5)]
   formula <- make_formula(1, 0)
   
   CQuantAR1_OOS <- qCQR(formula,x0,y0,x1,y1,x_test,QQ)
   
   #AR2
   x0 <- data_1[,-c(2,5)]
   y0 <- data_1[,2]
   x1 <- data_2[,-c(2,5)]
   y1 <- data_2[,2]
   x_test <- data_test[,-c(2,5)]
   formula <- make_formula(2, 0)
   
   CQuantAR2_OOS <- qCQR(formula,x0,y0,x1,y1,x_test,QQ)
   
   #AR2
   x0 <- data_1[,-2]
   y0 <- data_1[,2]
   x1 <- data_2[,-2]
   y1 <- data_2[,2]
   x_test <- data_test[,-2]
   formula <- make_formula(3, 0)
   
   CQuantAR3_OOS <- qCQR(formula,x0,y0,x1,y1,x_test,QQ)
  
  #--------------- CONFORMAL QUANTILE RANDOM FOREST
  
  # #AR1
  # x0 <- data_1[,-c(2,4,5)]
  # y0 <- data_1[,2]
  # x1 <- data_2[,-c(2,4,5)]
  # y1 <- data_2[,2]
  # x_test <- data_test[,-c(2,4,5)]
  # 
  # CQuantAR1_OOS <- qCQRF(x0,y0,x1,y1,x_test,QQ)
  # 
  # #AR2
  # x0 <- data_1[,-c(2,5)]
  # y0 <- data_1[,2]
  # x1 <- data_2[,-c(2,5)]
  # y1 <- data_2[,2]
  # x_test <- data_test[,-c(2,5)]
  # 
  # CQuantAR2_OOS <- qCQRF(x0,y0,x1,y1,x_test,QQ)
  # 
  # #AR2
  # x0 <- data_1[,-2]
  # y0 <- data_1[,2]
  # x1 <- data_2[,-2]
  # y1 <- data_2[,2]
  # x_test <- data_test[,-2]
  # 
  # CQuantAR3_OOS <- qCQRF(x0,y0,x1,y1,x_test,QQ)
  # 
  # 
  #---------------- CALIBRATION OF QR AND CQR
  
  coverageQRAR1 <- compute_coverage(Y_test, QuantAR1_OOS, QQ)
  coverageCQRAR1 <- compute_coverage(Y_test, CQuantAR1_OOS, QQ)
  
  coverageQRAR2 <- compute_coverage(Y_test, QuantAR2_OOS, QQ)
  coverageCQRAR2 <- compute_coverage(Y_test, CQuantAR2_OOS, QQ)
  
  coverageQRAR3 <- compute_coverage(Y_test, QuantAR3_OOS, QQ)
  coverageCQRAR3 <- compute_coverage(Y_test, CQuantAR3_OOS, QQ)
    
  ret <- list(coverageQRAR1, coverageCQRAR1, coverageQRAR2, coverageCQRAR2, coverageQRAR3, coverageCQRAR3)
  
  return(ret)
}



vector_alpha <- c(1,0.1,0.01) 
vector_n <- c(101,201,1001)

n2 = 100
n3 = 100 
n_simul <- n3
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(0.01, QQ)

count <- 2 #so I leave first row empty

for (n in vector_n){
  for (alpha in vector_alpha){
    


 # Setup parallel cluster
 cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
 clusterExport(cl, varlist=c("run_simulation","n","alpha")) # Export the simulation function to each cluster node
 clusterEvalQ(cl, { 
   library(readxl)
   library(quantreg)
   library(quantregForest)
   source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")
 }
 ) # Load required libraries in each cluster node, repeat as necessary for other libraries

# Run simulations in parallel
 results <- parLapply(cl, seeds, function(seed) {
#results <- lapply(seeds, function(seed) {
  set.seed(seed)
  run_simulation(n,alpha)
})


# Stop the cluster
stopCluster(cl)

# Extract results
coveragetotQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotQRAR2 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR2 <- matrix(NA,n3,length(QQ))
coveragetotQRAR3 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR3 <- matrix(NA,n3,length(QQ))



index = 1
for(res in results){
  coveragetotQRAR1[index,] <- res[[1]]
  coveragetotCQRAR1[index,] <- res[[2]]
  coveragetotQRAR2[index,] <- res[[3]]
  coveragetotCQRAR2[index,] <- res[[4]]
  coveragetotQRAR3[index,] <- res[[5]]
  coveragetotCQRAR3[index,] <- res[[6]]
  
  index <- index + 1
}

average_coverageQRAR1 <- compute_average_coverage(coveragetotQRAR1, QQ) #vector of length = length(QQ)
average_coverageQRAR2 <- compute_average_coverage(coveragetotQRAR2, QQ)
average_coverageQRAR3 <- compute_average_coverage(coveragetotQRAR3, QQ)

average_coverageCQRAR1 <- compute_average_coverage(coveragetotCQRAR1, QQ)
average_coverageCQRAR2 <- compute_average_coverage(coveragetotCQRAR2, QQ)
average_coverageCQRAR3 <- compute_average_coverage(coveragetotCQRAR3, QQ)

resultsPitSTQRAR1 <- compute_results(average_coverageQRAR1, n2*n3, QQ)
resultsPitSTQRAR2 <- compute_results(average_coverageQRAR2, n2*n3, QQ)
resultsPitSTQRAR3 <- compute_results(average_coverageQRAR3, n2*n3, QQ)
resultsPitSTCQRAR1 <- compute_results(average_coverageCQRAR1, n2*n3, QQ)
resultsPitSTCQRAR2 <- compute_results(average_coverageCQRAR2, n2*n3, QQ)
resultsPitSTCQRAR3 <- compute_results(average_coverageCQRAR3, n2*n3, QQ)

#------------ Write Results in the excel file


wb <- loadWorkbook(file_path)

# Access the worksheet 
sheet <- "Data"

#In the first column, write the parameters of the simulation:
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3, ", QR AR(1)", ",alpha =", alpha), startCol = 1, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3, ", CQR AR(1)", ",alpha =", alpha ), startCol = 1, startRow = count+1, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3, ", QR AR(2)", ",alpha =", alpha ), startCol = 1, startRow = count+2, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3, ", CQR AR(2)", ",alpha =", alpha ), startCol = 1, startRow = count+3, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3, ", QR AR(3)", ",alpha =", alpha ), startCol = 1, startRow = count+4, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3, ", CQR AR(3)", ",alpha =", alpha ), startCol = 1, startRow = count+5, colNames = FALSE)

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

#In the 5° column, the MAE will be placed
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR1$Quantile-resultsPitSTQRAR1$EmpiricalCoverage)), startCol = 5, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR1$Quantile-resultsPitSTCQRAR1$EmpiricalCoverage)), startCol = 5, startRow = count+1, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR2$Quantile-resultsPitSTQRAR2$EmpiricalCoverage)), startCol = 5, startRow = count+2, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR2$Quantile-resultsPitSTCQRAR2$EmpiricalCoverage)), startCol = 5, startRow = count+3, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR3$Quantile-resultsPitSTQRAR3$EmpiricalCoverage)), startCol = 5, startRow = count+4, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR3$Quantile-resultsPitSTCQRAR3$EmpiricalCoverage)), startCol = 5, startRow = count+5, colNames = FALSE)

# #In the 6°,7° column, the MAE+ and MAE- will be placed
# filtered_results <- resultsPitSTQRAR1[resultsPitSTQRAR1$EmpiricalCoverage > resultsPitSTQRAR1$Quantile, ]
# writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count, colNames = FALSE)
# filtered_results <- resultsPitSTCQRAR1[resultsPitSTCQRAR1$EmpiricalCoverage > resultsPitSTCQRAR1$Quantile, ]
# writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 6, startRow = count+1, colNames = FALSE)
# 

# filtered_results <- resultsPitSTQRAR1[resultsPitSTQRAR1$EmpiricalCoverage < resultsPitSTQRAR1$Quantile, ]
# writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count, colNames = FALSE)
# filtered_results <- resultsPitSTCQRAR1[resultsPitSTCQRAR1$EmpiricalCoverage < resultsPitSTCQRAR1$Quantile, ]
# writeData(wb, sheet = sheet, x = mean(abs(filtered_results$Quantile-filtered_results$EmpiricalCoverage)), startCol = 7, startRow = count+1, colNames = FALSE)

saveWorkbook(wb, file_path, overwrite = TRUE)


count <- count + 7 #lascio una riga vuota



#------------- PLOT CALIBRATION CURVES OF QR and CQR and DCP with AR(1) MODEL

# Create data frames for both datasets
df1 <- data.frame(
  Quantile = resultsPitSTCQRAR2$Quantile,
  EmpiricalCoverage = resultsPitSTCQRAR2$EmpiricalCoverage,
  Group = "CQR AR2"
)

df2 <- data.frame(
  Quantile = resultsPitSTQRAR2$Quantile,
  EmpiricalCoverage = resultsPitSTQRAR2$EmpiricalCoverage,
  Group = "QR AR2"
)

# Combine the data frames
df <- bind_rows(df1, df2)

# Define the colors
cqr_colors <- "#0000FF"
qr_colors <- "#FF0000"

# Create staircase plot
p <- ggplot(df, aes(x = Quantile, y = EmpiricalCoverage, color = Group)) +
  geom_step(direction = "vh", size = 1) + # Add staircase lines
  geom_point(size = 3) + # Add points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) + # Add diagonal line
  scale_color_manual(values = c("CQR AR2" = cqr_colors, "QR AR2" = qr_colors)) + # Manual color scale
  labs(title = paste("n1 =", n - 3, ", AR(2) Model",", alpha = ", alpha), x = "Quantile Levels", y = "Empirical Coverage") + # Add labels and title
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.15), # Position the legend at the bottom right
    legend.title = element_blank(),
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)
  )

# Print the plot
print(p)
# Save the plot as a PDF file
ggsave(filename = paste0("Cauchy_calibration_n", n - 3, "_AR(2)_QRegression.pdf"), plot = p, width = 7, height = 5)
}

}
