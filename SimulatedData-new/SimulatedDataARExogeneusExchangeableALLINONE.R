
#Qualunque sia il passo di sampling, non ottengo copertura, a meno di smettere di prevedere il 101-esimo punto e prevedere
#da 100 a 200


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

 
file_path <- "ALLINONE_Exogenous_Results_Exchangeable.xlsx"
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add a worksheet named "Data"
  addWorksheet(wb, "Data")
  
  # Save the workbook (this creates the file)
  saveWorkbook(wb, file_path, overwrite = TRUE)
 
 wb <- loadWorkbook(file_path)


run_simulation <- function(n,ratio_p_n,B){
  
  # Set parameters
  #n <- 31      # Number of data points for training. 3 will always be lost!
  n2 <- 100 * B    # Number of data points for testing
  p <- ratio_p_n * n / B
  phi_ar2 <- c(0.5, -0.2)   # AR(2) coefficients
  num_exog_vars <-  p     # Parameter to control the number of exogenous variables
  beta <- runif(num_exog_vars, 0.1, 1) # Randomly generated coefficients for exogenous variables
  means <- rnorm(p)
  sigma <- diag(runif(p))
  
  # Initialize variables
  Y_ar2 <- numeric(n+n2)  # Series for AR(2) data
  exog_vars <- matrix(nrow = n + n2, ncol = num_exog_vars)  # Matrix for exogenous variables
  
  exog_vars <- rmvnorm( n = n + n2, mean = means, sigma = sigma)

  # Initialize Y values
  Y_ar2[1:2] <- rnorm(2)  
  
  # Simulate AR(2) with exogenous variables for training data
  for (t in 3:(n+n2)) {
    Y_ar2[t] <- phi_ar2[1] * Y_ar2[t - 1] + phi_ar2[2] * Y_ar2[t - 2] +
      sum(beta * exog_vars[t, ]) + rt(1, df = 2)
  }
  
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])
  Y_lag3_ar2 <- c(NA, NA, NA, Y_ar2[1:(length(Y_ar2)-3)])
  
  

  # Prepare dataset for model
  indices <- c((1:3),(n+1):(n+n2))
  
   data_ar2 <- data.frame(I = 1, Y = Y_ar2[-indices], Y_lag1 = Y_lag1_ar2[-indices], Y_lag2 = Y_lag2_ar2[-indices], Y_lag3 = Y_lag3_ar2[-indices])
   for (j in 1:num_exog_vars) {
     data_ar2[paste0("X", j)] <- exog_vars[-indices, j]
   }
   
   data_ar2_test <- data.frame(I = 1, Y = Y_ar2[indices], Y_lag1 = Y_lag1_ar2[indices], Y_lag2 = Y_lag2_ar2[indices], Y_lag3 = Y_lag3_ar2[indices])
   for (j in 1:num_exog_vars) {
     data_ar2_test[paste0("X", j)] <- exog_vars[indices, j]
   }
   data_ar2_test <- data_ar2_test[4:nrow(data_ar2_test),]
   
  
  Y_test <- Y_ar2[c((n+1):(n+n2))]
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)
  
   
  
  #--------------- QUANTILE REGRESSION 
  
  QuantAR0_OOS <- matrix(0, n2, length(QQ))
  QuantAR1_OOS <- matrix(0,n2, length(QQ))
  QuantAR2_OOS <- matrix(0, n2, length(QQ))
  QuantAR3_OOS <- matrix(0,n2, length(QQ))

   for (jq in 1:length(QQ)) {  
    
     QR0 <- rq(make_formula(0, num_exog_vars), data = data_ar2, tau=QQ[jq])     
     #QuantAR0_OOS[jq] <- c(1, exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR0)
     QuantAR0_OOS[,jq] <- as.matrix(data_ar2_test[,-c(2:5)]) %*% coef(QR0)
     
     # QR1 <- rq(make_formula(1, num_exog_vars), data = data_ar2, tau=QQ[jq])     
      #[QuantAR1_OOS[jq] <- c(1, Y_ar2[n], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR1)
     #  
      QR2 <- rq(make_formula(2, num_exog_vars), data = data_ar2, tau=QQ[jq])
     #QuantAR2_OOS[jq] <- c(1, Y_ar2[n], Y_ar2[n-1], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR2)
      QuantAR2_OOS[,jq] <- as.matrix(data_ar2_test[,-c(2,5)]) %*% coef(QR2)
      
      
     # QR3 <- rq(make_formula(3, num_exog_vars), data = data_ar2, tau=QQ[jq])
     # QuantAR3_OOS[jq] <- c(1, Y_ar2[n], Y_ar2[n-1], Y_ar2[n-2], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR3)

  
  }
  
  #--------------- QUANTILE RANDOM FOREST 
  
  {
  # qrf_model1 <- quantregForest(x = (data_ar2[,-c(1,2,4,5)]), y = (data_ar2[,2]))
  # qrf_model2 <- quantregForest(x = (data_ar2[,-c(1,2,5)]), y = (data_ar2[,2]))
  # qrf_model3 <- quantregForest(x = (data_ar2[,-c(1,2)]), y = (data_ar2[,2]))
  
  #for (jq in 1:length(QQ)) {  
    
    #QuantAR1_OOS[jq] <- predict(qrf_model1, newdata = c(Y_ar2[n],exog_vars_n_plus_1), what = QQ[jq])
    #QuantAR2_OOS[jq] <- predict(qrf_model2, newdata = c(Y_ar2[n],Y_ar2[n-1],exog_vars_n_plus_1), what = QQ[jq])
    #QuantAR3_OOS[jq] <- predict(qrf_model3, newdata = c(Y_ar2[n],Y_ar2[n-1],Y_ar2[n-2],exog_vars_n_plus_1), what = QQ[jq])
    
  #}
  }
  
  #--------------- CONFORMAL QUANTILE REGRESSION
  
  full_length <- nrow(data_ar2)
  Q_low <- matrix(NA, nrow(data_ar2), length(QQ))
  Q_high <- matrix(NA, nrow(data_ar2), length(QQ))
  test_length = floor(full_length*50/100)
  data_ar2_1 <- data_ar2[1:test_length,] #train
  data_ar2_2 <- data_ar2[(test_length + 1) : full_length,] #calibration
  
  CQuantAR0_OOS <- matrix(0,n2, length(QQ))
  CQuantAR1_OOS <- matrix(0,n2, length(QQ))
  CQuantAR2_OOS <- matrix(0,n2, length(QQ))
  CQuantAR3_OOS <- matrix(0,n2, length(QQ))
  
  #AR0 e AR2
  {  
       x0 <- data_ar2_1[,-c(2,3,4,5)]
       y0 <- data_ar2_1[,2]
       x1 <- data_ar2_2[,-c(2,3,4,5)]
       y1 <- data_ar2_2[,2]
      # x_test <- c(1,exog_vars_n_plus_1)
       x_test <- as.matrix(data_ar2_test[,-c(2,3,4,5)])
       
       CQuantAR0_OOS <- qCQR(make_formula(0, num_exog_vars),x0,y0,x1,y1,x_test,QQ) 
    
     # for (jq in 1:length(QQ)) {  
     #   
     #   Q_low[1:nrow(data_ar2_2), jq] <- -Inf 
     #   QR <- rq(make_formula(0,num_exog_vars), data = data_ar2_1, tau = QQ[jq])
     #   Q_high[1 : nrow(data_ar2_2), jq] <- as.matrix(data_ar2_2[,-c(2:5)]) %*% coef(QR) 
     #   
     #   #Q_high[1 : nrow(data_ar2_2), jq] <- predict(qrf_model, newdata = as.matrix(data_ar2_2[,-c(1,2)]), what = QQ[jq])
     #   
     #   # Initialize a vector for errors
     #   E_i <- rep(NA, nrow(data_ar2_2))
     #   
     #   # Calculate errors for each point in the test set I2
     #   for (i in 1:length(E_i)) {
     #     E_i[i] <- max(Q_low[i, jq] - data_ar2_2[i,2], data_ar2_2[i,2] - Q_high[i, jq])
     #   }
     #   
     #   # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
     #   quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(data_ar2_2)))
     #   
     #   CQuantAR0_OOS[,jq] <- as.matrix(data_ar2_test[,-c(2:5)]) %*% coef(QR) + quantile_E
     #   #CQuant_OOS[,jq] <- predict(qrf_model, newdata = as.matrix(data_ar2_test[,-c(1,2)]), what = QQ[jq]) + quantile_E
     # 
     # }
    
    
    
     
      # x0 <- data_ar2_1[,-c(2,4,5)]
      # y0 <- data_ar2_1[,2]
      # x1 <- data_ar2_2[,-c(2,4,5)]
      # y1 <- data_ar2_2[,2]
      # x_test <- c(1,Y_ar2[n],exog_vars_n_plus_1)
      # 
      # CQuantAR1_OOS <- qCQR(make_formula(1, num_exog_vars),x0,y0,x1,y1,x_test,QQ) 
      # 
     
      x0 <- data_ar2_1[,-c(2,5)]
      y0 <- data_ar2_1[,2]
      x1 <- data_ar2_2[,-c(2,5)]
      y1 <- data_ar2_2[,2]
      x_test <- as.matrix(data_ar2_test[,-c(2,5)])
      
      CQuantAR2_OOS <- qCQR(make_formula(2, num_exog_vars),x0,y0,x1,y1,x_test,QQ) 
    
     
     # for (jq in 1:length(QQ)) {  
     #   
     #   Q_low[1:nrow(data_ar2_2), jq] <- -Inf 
     #   QR <- rq(make_formula(2,num_exog_vars), data = data_ar2_1, tau = QQ[jq])
     #   Q_high[1 : nrow(data_ar2_2), jq] <- as.matrix(data_ar2_2[,-c(2,5)]) %*% coef(QR) 
     #   
     #   #Q_high[1 : nrow(data_ar2_2), jq] <- predict(qrf_model, newdata = as.matrix(data_ar2_2[,-c(1,2)]), what = QQ[jq])
     #   
     #   # Initialize a vector for errors
     #   E_i <- rep(NA, nrow(data_ar2_2))
     #   
     #   # Calculate errors for each point in the test set I2
     #   for (i in 1:length(E_i)) {
     #     E_i[i] <- max(Q_low[i, jq] - data_ar2_2[i,2], data_ar2_2[i,2] - Q_high[i, jq])
     #   }
     #   
     #   # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
     #   quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(data_ar2_2)))
     #   
     #   CQuantAR2_OOS[,jq] <- as.matrix(data_ar2_test[,-c(2,5)]) %*% coef(QR) + quantile_E
     #   #CQuant_OOS[,jq] <- predict(qrf_model, newdata = as.matrix(data_ar2_test[,-c(1,2)]), what = QQ[jq]) + quantile_E
     #   
     # }
    
    
    
    #  x0 <- data_ar2_1[,-2]
    #  y0 <- data_ar2_1[,2]
    #  x1 <- data_ar2_2[,-2]
    #  y1 <- data_ar2_2[,2]
    #  x_test <- c(1,Y_ar2[n],Y_ar2[n-1],Y_ar2[n-2],exog_vars_n_plus_1)
    #  
    #  CQuantAR3_OOS <- qCQR(make_formula(3, num_exog_vars),x0,y0,x1,y1,x_test,QQ) 
  }
  
  #--------------- CONFORMAL QUANTILE RANDOM FOREST
  
  {
  #  
  #  x0 <- data_ar2_1[,-c(1,2,4,5)]
  #  y0 <- data_ar2_1[,2]
  #  x1 <- data_ar2_2[,-c(1,2,4,5)]
  #  y1 <- data_ar2_2[,2]
  #  x_test <- c(Y_ar2[n],exog_vars_n_plus_1)
  #  
  #  CQuantAR1_OOS <- qCQRF(x0,y0,x1,y1,x_test,QQ) 
  #  
  #  x0 <- data_ar2_1[,-c(1,2,5)]
  #  y0 <- data_ar2_1[,2]
  #  x1 <- data_ar2_2[,-c(1,2,5)]
  #  y1 <- data_ar2_2[,2]
  #  x_test <- c(Y_ar2[n],Y_ar2[n-1],exog_vars_n_plus_1)
  #  
  #  CQuantAR2_OOS <- qCQRF(x0,y0,x1,y1,x_test,QQ) 
  #  
  #  x0 <- data_ar2_1[,-c(1,2)]
  #  y0 <- data_ar2_1[,2]
  #  x1 <- data_ar2_2[,-c(1,2)]
  #  y1 <- data_ar2_2[,2]
  #  x_test <- c(Y_ar2[n],Y_ar2[n-1],Y_ar2[n-2],exog_vars_n_plus_1)
  #  
  #  CQuantAR3_OOS <- qCQRF(x0,y0,x1,y1,x_test,QQ) 
   }

  #--------------- DISTRIBUTIONAL CONFORMAL PREDICTION
  
  {
    
    
    # source("temp1.R")
    # source("full.R")
    # funs=rq.funs()
    # 
    # predict.fun=funs$predict.fun
    # train.fun=funs$train.fun
    # 
     DQuantAR2_OOS <- matrix(0, n2,length(QQ))
    # 
    #  y <- data_ar2[,2]
    #  x <- data_ar2[,-c(1,2,5)]
    #  x0 <- c(Y_ar2[n],Y_ar2[n-1],exog_vars[n+1, 1:num_exog_vars])
    # 
    #  for (jq in 1:length(QQ)) {
    #    test = dist.conformal.pred(x,y,x0,train.fun = train.fun, predict.fun = predict.fun, verbose = T, alpha = 1 - QQ[jq])
    #    DQuantAR2_OOS[jq] <- test$up[1,1]
    #    cat("Upper level of CI:", test$up[1, 1], "- Lower level of CI:", test$lo[1, 1], "\n")
    #    cat("Minimum of vector Y:", min(y), "\n")
    #  }
    #  
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
  
  coverageDCPAR2 <- compute_coverage(Y_test, DQuantAR2_OOS, QQ)
  
#AR1
  #  
  #  coverageQRAR0 <- rep(NA,length(QQ))
  #  coverageCQRAR0 <- rep(NA,length(QQ))
  # # 
  #  PitQROOS <- rep(NA, length(Y_test))
  #  PitCQROOS <- rep(NA, length(Y_test))
  #  
  #  
  #  for (jq in 1:length(QQ)) {
  #    
  #    for (i in 1:length(Y_test)){
  #      if (Y_test[i] <= QuantAR0_OOS[i,jq]) {
  #        PitQROOS[i]  <- 1 
  #      }
  #      else 
  #        PitQROOS[i]  <- 0
  #    }
  #    
  #    coverageQRAR0[jq] <- sum(PitQROOS) / length(PitQROOS)
  #    
  #    
  #    for (i in 1:length(Y_test)){
  #      if (Y_test[i] <= CQuantAR0_OOS[i,jq]) {
  #        PitCQROOS[i]  <- 1 
  #      }
  #      else 
  #        PitCQROOS[i]  <- 0
  #    }
  #  
  #    coverageCQRAR0[jq] <- sum(PitCQROOS) / length(PitCQROOS)
  #  
  # 
  #  }
  #  
  #  
  #  #AR2 
  #  
  #  coverageQRAR2 <- rep(NA,length(QQ))
  #  coverageCQRAR2 <- rep(NA,length(QQ))
  #  
  #  PitQROOS <- rep(NA, length(Y_test))
  #  PitCQROOS <- rep(NA, length(Y_test))
  #  
  #  
  #  for (jq in 1:length(QQ)) {
  #    
  #    for (i in 1:length(Y_test)){
  #      if (Y_test[i] <= QuantAR2_OOS[i,jq]) {
  #        PitQROOS[i]  <- 1 
  #      }
  #      else 
  #        PitQROOS[i]  <- 0
  #    }
  #    
  #    coverageQRAR2[jq] <- sum(PitQROOS) / length(PitQROOS)
  #    
  #    for (i in 1:length(Y_test)){
  #      if (Y_test[i] <= CQuantAR2_OOS[i,jq]) {
  #        PitCQROOS[i]  <- 1 
  #      }
  #      else 
  #        PitCQROOS[i]  <- 0
  #    }
  #    
  #    coverageCQRAR2[jq] <- sum(PitCQROOS) / length(PitCQROOS)
  #    
  #  }
  #  
   
  
  
  
  
  ret <- list(coverageQRAR1, coverageQRAR2, coverageQRAR3, coverageCQRAR1, coverageCQRAR2, coverageCQRAR3, coverageDCPAR2,
              coverageQRAR0, coverageCQRAR0)
  
  return(ret)
}

B <- 1
vector_n <- c(101*B,201*B)
vector_p_n <- c(0.05,0.2,0.3,0.4)


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
clusterExport(cl, varlist=c("run_simulation","n","p_n","B")) # Export the simulation function to each cluster node
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
  run_simulation(n,p_n,B)
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

coveragetotDCPAR2 <- matrix(NA,n3,length(QQ))





index = 1
for(res in results){
  coveragetotQRAR1[index,] <- res[[1]]
  coveragetotQRAR2[index,] <- res[[2]]
  coveragetotQRAR3[index,] <- res[[3]]
  
  coveragetotCQRAR1[index,] <- res[[4]]
  coveragetotCQRAR2[index,] <- res[[5]]
  coveragetotCQRAR3[index,] <- res[[6]]
  
  coveragetotDCPAR2[index,] <- res[[7]]
  
  coveragetotQRAR0[index,] <- res[[8]]
  coveragetotCQRAR0[index,] <- res[[9]]
  
  
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

average_coverageQRAR0 <- compute_average_coverage(coveragetotQRAR0, QQ) #vector of length = length(QQ)
average_coverageCQRAR0 <- compute_average_coverage(coveragetotCQRAR0, QQ) #vector of length = length(QQ)


resultsPitSTQRAR0 <- compute_results(average_coverageQRAR0, n2*n3, QQ)
resultsPitSTQRAR1 <- compute_results(average_coverageQRAR1, n2*n3, QQ)
resultsPitSTQRAR2 <- compute_results(average_coverageQRAR2, n2*n3, QQ)
resultsPitSTQRAR3 <- compute_results(average_coverageQRAR3, n2*n3, QQ)
resultsPitSTCQRAR0 <- compute_results(average_coverageCQRAR0, n2*n3, QQ)
resultsPitSTCQRAR1 <- compute_results(average_coverageCQRAR1, n2*n3, QQ)
resultsPitSTCQRAR2 <- compute_results(average_coverageCQRAR2, n2*n3, QQ)
resultsPitSTCQRAR3 <- compute_results(average_coverageCQRAR3, n2*n3, QQ)
  



#------------ Write Results in the excel file


wb <- loadWorkbook(file_path)

# Access the worksheet 
sheet <- "Data"

#In the first column, write the parameters of the simulation:
writeData(wb, sheet = sheet, x = paste("n1 = ", n/B-3,", p/n = ", p_n, ", QR AR(1)" ), startCol = 1, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n/B-3,", p/n = ", p_n, ", CQR AR(1)" ), startCol = 1, startRow = count+1, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n/B-3,", p/n = ", p_n, ", QR AR(2)" ), startCol = 1, startRow = count+2, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n/B-3,", p/n = ", p_n, ", CQR AR(2)" ), startCol = 1, startRow = count+3, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n/B-3,", p/n = ", p_n, ", QR AR(3)" ), startCol = 1, startRow = count+4, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n/B-3,", p/n = ", p_n, ", CQR AR(3)" ), startCol = 1, startRow = count+5, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n/B-3,", p/n = ", p_n, ", QR AR(0)" ), startCol = 1, startRow = count+6, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n/B-3,", p/n = ", p_n, ", CQR AR(0)" ), startCol = 1, startRow = count+7, colNames = FALSE)

#In the 2�,3�,4� columns, put how many inside, below and above CI
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

writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR0$IsWithinCI) , startCol = 2, startRow = count + 6, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR0$EmpiricalCoverage > resultsPitSTQRAR0$Quantile & !resultsPitSTQRAR0$IsWithinCI) , startCol = 3, startRow = count + 6, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR0$EmpiricalCoverage < resultsPitSTQRAR0$Quantile & !resultsPitSTQRAR0$IsWithinCI) , startCol = 4, startRow = count + 6, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR0$IsWithinCI) , startCol = 2, startRow = count + 7, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR0$EmpiricalCoverage > resultsPitSTCQRAR0$Quantile & !resultsPitSTCQRAR0$IsWithinCI) , startCol = 3, startRow = count + 7, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR0$EmpiricalCoverage < resultsPitSTCQRAR0$Quantile & !resultsPitSTCQRAR0$IsWithinCI) , startCol = 4, startRow = count + 7, colNames = FALSE)

#In the 5� column, the MAE will be placed
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR1$Quantile-resultsPitSTQRAR1$EmpiricalCoverage)), startCol = 5, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR1$Quantile-resultsPitSTCQRAR1$EmpiricalCoverage)), startCol = 5, startRow = count+1, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR2$Quantile-resultsPitSTQRAR2$EmpiricalCoverage)), startCol = 5, startRow = count+2, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR2$Quantile-resultsPitSTCQRAR2$EmpiricalCoverage)), startCol = 5, startRow = count+3, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR3$Quantile-resultsPitSTQRAR3$EmpiricalCoverage)), startCol = 5, startRow = count+4, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR3$Quantile-resultsPitSTCQRAR3$EmpiricalCoverage)), startCol = 5, startRow = count+5, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR0$Quantile-resultsPitSTQRAR0$EmpiricalCoverage)), startCol = 5, startRow = count+6, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR0$Quantile-resultsPitSTCQRAR0$EmpiricalCoverage)), startCol = 5, startRow = count+7, colNames = FALSE)


#In the 6�,7� column, the MAE+ and MAE- will be placed

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
  }
}

saveWorkbook(wb, file_path, overwrite = TRUE)

