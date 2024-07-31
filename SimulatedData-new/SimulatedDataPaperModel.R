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


source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")

file_path <- "ALLINONE_Paper_Model_Results.xlsx"
# Create a new workbook
wb <- createWorkbook()

# Add a worksheet named "Data"
addWorksheet(wb, "Data")

# Save the workbook (this creates the file)
saveWorkbook(wb, file_path, overwrite = TRUE)

wb <- loadWorkbook(file_path)
# da qui comincia la simulazione


run_simulation <- function(n, step){


  n2 <- 1000
  tot_length <- n + n2
  if (step == 1)
  {
    sigmah <- 0.04
    v <- 9.85
    phi <- -0.39
    beta0 <- -1.05
    beta1 <- -0.85
    sigmatheta <- 0.016
  } else if (step == 4) {
    
    sigmah <- 0.069
    v <- 19.56
    phi <- -0.25
    beta0 <- -1.47
    beta1 <- 0.3
    sigmatheta <- 0.017
  }
  
  theta <- numeric(tot_length)
  theta[1] <- 0 #MY CHOICE
  h <- numeric(tot_length)
  h[1] <- 0 #MY CHOICE
  NFCI <- numeric(tot_length)

  
  file_path <- "DataVulnerabilityAppendix.xls"
  
  # Read the file
  data <- read_excel(file_path)
  data <- data[,1:3]
  NFCI_real <- data$NFCI
  y_real <- data$A191RL1Q225SBEA
  
  #------ Fit an AR(1) model to NFCI data, then generate a new time series using the fitted model
  #fit <- ar(NFCI_real, order.max = 1)
  fit <- Arima(NFCI_real, order=c(1,0,0))  
  
  # Extract  coefficient
  ar_coefficients <- fit$coef["ar1"]  # For AR(1)
  ma_coefficients <- fit$coef["ma1"]  # For MA(1)
  if(is.na(ma_coefficients))
    ma_coefficients <- 0
  c <- fit$coef["intercept"]
  NFCI <- numeric(tot_length)
  NFCI[1] <- -1 #MY CHOICE
  last_error <- rnorm(1, mean = 0, sd = sqrt(fit$sigma2))
  
  # Simulate new values
  for(i in 2:tot_length) {
    new_error <- rnorm(1, mean = 0, sd = sqrt(fit$sigma2))
    NFCI[i] <- c + ar_coefficients * NFCI[i-1] + ma_coefficients * last_error + new_error
    last_error <- new_error
  }

  
  #------ Generate omega(t), nu(t)
  omega <- rnorm(tot_length, 0, sigmatheta)
  nu <- rnorm(tot_length, 0, sigmah)
  
  #------ Generate theta(t) = phi*theta(t-1) + beta0 + beta1*NFCI(t) + omega(t)
  for(t in 2:tot_length) {
    theta[t] <- phi * theta[t-1] + beta0 + beta1*NFCI[t] + omega[t]
  }
  
  #------ Compute c11(v) and E[u(t)]
  c11 <- sqrt(0.5 * v) * gamma(0.5 * (v - 1)) / gamma(0.5 * v)
  Expv <- numeric(tot_length)
  for(t in 1:tot_length) {
    Expv[t] <- c11 * theta[t]
  }
  
  #------ Generate u(t) sampled from NCT (v, theta(t))
  u <- numeric(tot_length)
  for(t in 1:tot_length) {
    #u[t] <- rnct(n = 1, df = v, ncp = theta[t])
    u[t] <- rt(n = 1, df = v, ncp = theta[t])
  }
  
  #------ Generate eps(t) = u(t) - E[u(t)]
  eps <- u - Expv
  
  #------ Generate h(t) = h(t-1) + nu(t)
  for(t in 2:tot_length) {
    h[t] <- h[t-1] + nu[t]
  }
  
  #------ Generate y(t) = exp(h(t)/2) * eps(t), which is my GDP growth
  y <- exp(h/2) * eps
  
  Yh <- stats::filter(y, rep(1/step, step), sides=1) # moving average of the h past values
  
  
  QQ <- seq(0.05, 0.95,by = 0.05) #miscoverage rate
  QQ <- c(0.01,QQ)
    
  indices <- c((n+1):(n+n2))
  
  Z <- cbind(1,Yh[(step+1):n], NFCI[1:(n-step)], y[1:(n-step)])
  ZGDPonly <- cbind(1,Yh[(step+1):n],y[1:(n-step)])
  Z <- as.matrix(Z)
  ZGDPonly <- as.matrix(ZGDPonly)
  
  indices <- c((n+1):(n+n2))
  data <- data.frame(I = 1, Y = Yh[(step+1):n], NFCI = NFCI[1:(n-step)], gdp = y[1:(n-step)])
  data_test <- data.frame(I = 1, Y = Yh[(n+1):(n+n2)], NFCI = NFCI[(n-step+1):(n+n2-step)], gdp = y[(n-step+1):(n+n2-step)])
  
  # DA QUI DEVO SISTMEMAREEEE
  
  Quant_OOS <- matrix(0, n2, length(QQ))
  Quantgdponly_OOS <- matrix(0, n2, length(QQ))
  
  for (jq in 1:length(QQ)) {  
    
    QR <- rq(Y ~ NFCI + gdp, data = data, tau = QQ[jq])     
    Quant_OOS[,jq] <- as.matrix(data_test[,-2]) %*% coef(QR)
    
    QRgdponly <- rq(Y ~ gdp, data = data, tau = QQ[jq])     
    Quantgdponly_OOS[,jq] <- as.matrix(data_test[,-c(2,3)]) %*% coef(QRgdponly)
    
  }
  
  
  full_length <- nrow(data)
  test_length = floor(full_length*50/100)
  data_1 <- data[1:test_length,] #train
  data_2 <- data[(test_length + 1) : full_length,] #calibration
  

    x0 <- data_1[,-2]
    y0 <- data_1[,2]
    x1 <- data_2[,-2]
    y1 <- data_2[,2]
    x_test <- data_test[,-2]
    formula <- "Y ~ NFCI + gdp"
    
    CQuant_OOS <- qCQR(formula,x0,y0,x1,y1,x_test,QQ) 
    
    x0 <- data_1[,-c(2,3)]
    y0 <- data_1[,2]
    x1 <- data_2[,-c(2,3)]
    y1 <- data_2[,2]
    x_test <- data_test[,-c(2,3)]
    formula <- "Y ~ gdp"
    
    
    CQuantgdponly_OOS <- qCQR(formula,x0,y0,x1,y1,x_test,QQ) 

    
    Y_test <- data_test$Y
    coverageQR <- compute_coverage(Y_test, Quant_OOS, QQ)
    coverageQRgdponly <- compute_coverage(Y_test, Quantgdponly_OOS, QQ)
    coverageCQR <- compute_coverage(Y_test, CQuant_OOS, QQ)
    coverageCQRgdponly <- compute_coverage(Y_test, CQuantgdponly_OOS, QQ)
  
    ret <- list(coverageQR, coverageQRgdponly, coverageCQR, coverageCQRgdponly)
    
    return(ret)
    
    #in teoria pronta a essere eseguita
  {
  
  {
  #------------ Now I want to generate n2 = 100 times the next number of the time series, namely n + 1 = 99
  
  new_error <- rnorm(n2, mean = 0, sd = sqrt(fit$sigma2))
  NFCIn2 <- c + ar_coefficients * NFCI[n] + ma_coefficients * last_error + new_error # n2 different values of NFCI[n+1]
  
  omegan2 <- rnorm(n2, 0, sigmatheta)
  nun2 <- rnorm(n2, 0, sigmah)
  
  thetan2 <- phi * theta[n] + beta0 + beta1*NFCIn2 + omegan2 # n2 different values of theta[n+1]
  Expvn2 <- c11 * thetan2
  
  un2 <- numeric(n2)
  for(j in 1:n2) {
    un2[j] <- rt(n = 1, df = v, ncp = thetan2[j])
  }
  
  epsn2 <- un2 - Expvn2
  
  hn2 <- h[n]+nun2
  
  yn2 <- exp(hn2/2) * epsn2
  
  Yh_test = (Yh[n-2]+Yh[n-1]+Yh[n]+yn2)/4 
} #data preparation
  
  #-------------------------------------------------------------------------------------------
  
  
  #QUANTILE REGRESSION 
  
  QQ <- seq(0.05, 0.95,by = 0.05) #miscoverage rates
  QQ <- c(QQ, 0.99)

  QRcov <- matrix(NA,nrow = length(Yh_test),ncol = length(QQ))
  QRGDPonlycov <- matrix(NA,nrow = length(Yh_test),ncol = length(QQ))
  
  
  #OUT OF SAMPLE
    for (jq in 1:length(QQ)) {
      
      #FULL CONDITIONAL
       QRfull <- rq(Yh[(step+1):length(Yh)] ~ Z[1:(length(Yh) - step),-1], tau = (1 - QQ[jq]))
       QRcov[, jq] <- (Yh_test <= as.numeric(Z[(length(Yh) - step + 1), ] %*% coef(QRfull)))
        
      
      #CONDITIONAL, GDP ONLY
       QRgdp <- rq(Yh[(step+1):length(Yh)] ~ ZGDPonly[1:(length(Yh) - step),-1], tau=(1 - QQ[jq]))
       QRGDPonlycov[, jq] <- (Yh_test <= as.numeric(ZGDPonly[(length(Yh) - step + 1), ] %*% coef(QRgdp)))
    }
  
  
  #-------------------------------------------------------------------------------------------
  
  
  # CONFORMALIZED QUANTILE REGRESSION
  
  #OUT OF SAMPLE
  CQRcov <- matrix(NA,nrow = length(Yh_test),ncol = length(QQ))
  CQRGDPonlycov <- matrix(NA,nrow = length(Yh_test),ncol = length(QQ))
  full_length <- length(Yh[(step+1):length(Yh)])
  test_length = full_length*50/100
  Yh1 <- Yh[(step+1):(test_length+step)]
  Yh2 <- Yh[(test_length +step +1):(full_length+step)] #Yh1 and Yh2 correctly have same dimension
  Z1 <- Z[1:test_length,]
  Z2 <- Z[(test_length+1):(length(Yh) - step),] #Z1 and Z2 correctly have same dimension
  ZGDPonly1 <- ZGDPonly[1:test_length,]
  ZGDPonly2 <- ZGDPonly[(test_length+1):(length(Yh) - step),]
    
  for (jq in 1:length(QQ)) {  
      
      #FULL CONDITIONAL
      Y0 <- Yh1
      X0 <- Z1[,-1]
      
      Y1 <- Yh2
      X1 <- Z2[,-1]
      
      Y.test <- Yh_test
      X.test <- matrix(c(rep(Z[(length(Yh) - step + 1), 2], n2),rep(Z[(length(Yh) - step + 1), 3], n2)),nrow = n2, ncol = 2, byrow = FALSE)  
      
      resfull <- cqrmod(Y0,X0,Y1,X1,Y.test,X.test, QQ[jq])
      CQRcov[,jq] <- resfull$cov.o
      
      
      #CONDITIONAL, GDP ONLY
      Y0 <- Yh1
      X0 <- ZGDPonly1[,-1]
      
      Y1 <- Yh2
      X1 <- ZGDPonly2[,-1]
      
      Y.test <- Yh_test
      X.test <- rep(ZGDPonly[(length(Yh) - step + 1),-1 ],n2)
      
      resgdp <- cqrmod(Y0,X0,Y1,X1,Y.test,X.test, QQ[jq])
      CQRGDPonlycov[,jq] <- resgdp$cov.o
    
    }
  

  ret <- list(QRcov, QRGDPonlycov, CQRcov, CQRGDPonlycov)
  
  return (ret)
} #vecchia, non più utile


}



vector_step <- c(1,4)
vector_n <- c(101,201,1001)
n2 <- 1000
n3 <- 1000
n_simul <- n3
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

QQ <- seq(0.05, 0.95,by = 0.05) #miscoverage rate
QQ <- c(0.01,QQ)



count <- 2 #so I leave first row empty



for(step in vector_step){
  for(n in vector_n){
    

# Setup parallel cluster
cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
clusterExport(cl, varlist=c("run_simulation","step","n")) # Export the simulation function to each cluster node
clusterEvalQ(cl, { 
  library(readxl)
  library(quantreg)
  library(forecast)
  source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")
}
  ) # Load required libraries in each cluster node, repeat as necessary for other libraries

# Run simulations in parallel
results <- parLapply(cl, seeds, function(seed) {
  set.seed(seed)
  run_simulation(n,step)
})

# Stop the cluster
stopCluster(cl)

# Extract results
coveragetotQR <- matrix(NA,n3,length(QQ))
coveragetotCQR <- matrix(NA,n3,length(QQ))
coveragetotQRgdponly <- matrix(NA,n3,length(QQ))
coveragetotCQRgdponly <- matrix(NA,n3,length(QQ))




index = 1
for(res in results){
  coveragetotQR[index,] <- res[[1]]
  coveragetotQRgdponly[index,] <- res[[2]]
  coveragetotCQR[index,] <- res[[3]]
  coveragetotCQRgdponly[index,] <- res[[4]]
  
  index <- index + 1
}



# RESULTS EVALUATION

# Calculate average coverage across all iterations and for each quantile
average_coverageQR <- compute_average_coverage(coveragetotQR, QQ) #vector of length = length(QQ)
average_coverageQRgdponly <- compute_average_coverage(coveragetotQRgdponly, QQ) #vector of length = length(QQ)

average_coverageCQR <- compute_average_coverage(coveragetotCQR, QQ) #vector of length = length(QQ)
average_coverageCQRgdponly <- compute_average_coverage(coveragetotCQRgdponly, QQ) #vector of length = length(QQ)



resultsPitSTQR <- compute_results(average_coverageQR, n2*n3, QQ)
resultsPitSTQRgdponly <- compute_results(average_coverageQRgdponly, n2*n3, QQ)

resultsPitSTCQR <- compute_results(average_coverageCQR, n2*n3, QQ)
resultsPitSTCQRgdponly <- compute_results(average_coverageCQRgdponly, n2*n3, QQ)




#------------ Write Results in the excel file


wb <- loadWorkbook(file_path)

# Access the worksheet 
sheet <- "Data"

#In the first column, write the parameters of the simulation:
writeData(wb, sheet = sheet, x = paste("n1 = ", n,", step = ", step, ", QR" ), startCol = 1, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n,", step = ", step, ", QRgdponly" ), startCol = 1, startRow = count+1, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n,", step = ", step, ", CQR" ), startCol = 1, startRow = count+2, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n,", step = ", step, ", CQRgdponly" ), startCol = 1, startRow = count+3, colNames = FALSE)

#In the 2°,3°,4° columns, put how many inside, below and above CI

writeData(wb, sheet = sheet, x = sum(resultsPitSTQR$IsWithinCI) , startCol = 2, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQR$EmpiricalCoverage > resultsPitSTQR$Quantile & !resultsPitSTQR$IsWithinCI) , startCol = 3, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQR$EmpiricalCoverage < resultsPitSTQR$Quantile & !resultsPitSTQR$IsWithinCI) , startCol = 4, startRow = count, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTQRgdponly$IsWithinCI) , startCol = 2, startRow = count + 1, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRgdponly$EmpiricalCoverage > resultsPitSTQRgdponly$Quantile & !resultsPitSTQRgdponly$IsWithinCI) , startCol = 3, startRow = count + 1, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRgdponly$EmpiricalCoverage < resultsPitSTQRgdponly$Quantile & !resultsPitSTQRgdponly$IsWithinCI) , startCol = 4, startRow = count + 1, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQR$IsWithinCI) , startCol = 2, startRow = count + 2, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQR$EmpiricalCoverage > resultsPitSTCQR$Quantile & !resultsPitSTCQR$IsWithinCI) , startCol = 3, startRow = count + 2, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQR$EmpiricalCoverage < resultsPitSTCQR$Quantile & !resultsPitSTCQR$IsWithinCI) , startCol = 4, startRow = count + 2, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRgdponly$IsWithinCI) , startCol = 2, startRow = count + 3, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRgdponly$EmpiricalCoverage > resultsPitSTCQRgdponly$Quantile & !resultsPitSTCQRgdponly$IsWithinCI) , startCol = 3, startRow = count + 3, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRgdponly$EmpiricalCoverage < resultsPitSTCQRgdponly$Quantile & !resultsPitSTCQRgdponly$IsWithinCI) , startCol = 4, startRow = count + 3, colNames = FALSE)

#In the 5° column, the MAE will be placed

writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQR$Quantile-resultsPitSTQR$EmpiricalCoverage)), startCol = 5, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRgdponly$Quantile-resultsPitSTQRgdponly$EmpiricalCoverage)), startCol = 5, startRow = count+1, colNames = FALSE)

writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQR$Quantile-resultsPitSTCQR$EmpiricalCoverage)), startCol = 5, startRow = count+2, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRgdponly$Quantile-resultsPitSTCQRgdponly$EmpiricalCoverage)), startCol = 5, startRow = count+3, colNames = FALSE)

saveWorkbook(wb, file_path, overwrite = TRUE)


count <- count + 5 #lascio una riga vuota



{
  x1 <- resultsPitSTCQR$Quantile
  y1 <- resultsPitSTCQR$EmpiricalCoverage
  x2 <- resultsPitSTQR$Quantile
  y2 <- resultsPitSTQR$EmpiricalCoverage
  
  # Create the plot
  plot(x1, y1, type = "n", xlab = "Quantile", ylab = "Empirical Coverage", main = paste("n1 = ", n-3, "step = ", step))
  
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
  }
}




  
  
  
  
  
  
  
  
  
  
  # 
  # 
  # 
  # 
  # PITtest_env <- new.env()
  # source("PITtest.r",local = PITtest_env)
  # rstestboot_env <- new.env()
  # source("rstestboot.r",local = rstestboot_env)
  # 
  # # PROBABILITY INTEGRAL TRANSFORM (PIT) FOR QUANTILE REGRESSION, FULL CONDITIONAL MODEL
  # rvec <- seq(0, 1, by = 0.001)
  # zST_ecdf1 <- PITtest_env$PITtest(PitST_OOStot, rvec)
  # 
  # 
  # # PROBABILITY INTEGRAL TRANSFORM (PIT) FOR CONFORMAL QUANTILE REGRESSION, FULL CONDITIONAL MODEL
  # 
  # rvec <- seq(0, 1, by = 0.001)
  # CzST_ecdf1 <- PITtest_env$PITtest(CPitST_OOStot, rvec)
  # 
  # 
  # # Plot PIT for quantile regression vs. conformal quantile regression 
  # plot(rvec, CzST_ecdf1, type = 'l', col = 'blue', xlab = 'tau', ylab = 'Empirical CDF', main = 'Full Conditional Model')
  # lines(rvec, zST_ecdf1, type = 'l', col = 'red')
  # lines(rvec, rvec , col = 'black',lty=2)
  # abline(h = seq(0, 1, by = 0.05), col = "lightgray", lty = "dotted")
  # abline(v = seq(0, 1, by = 0.05), col = "lightgray", lty = "dotted")
  # 
  # 
  # legend('bottomright', legend = c('Conformal Quantile Regression', 'Quantile regression'), cex = 1,fill = c('blue', 'red'))
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # 
  # # PROBABILITY INTEGRAL TRANSFORM (PIT) FOR QUANTILE REGRESSION, GDP ONLY MODEL
  # rvec <- seq(0, 1, by = 0.001)
  # zSTGDPonly_ecdf1 <- PITtest_env$PITtest(PitSTGDPonly_OOStot, rvec)
  # 
  # 
  # # PROBABILITY INTEGRAL TRANSFORM (PIT) FOR  CONFORMAL QUANTILE REGRESSION, GDP ONLY MODEL
  # 
  # rvec <- seq(0, 1, by = 0.001)
  # CzSTGDPonly_ecdf1 <- PITtest_env$PITtest(CPitSTGDPonly_OOStot, rvec)
  # 
  # 
  # # Plot PIT for quantile regression vs. conformal quantile regression 
  # plot(rvec, CzSTGDPonly_ecdf1, type = 'l', col = 'blue', xlab = 'tau', ylab = 'Empirical CDF', main = 'GDP ONLY MODEL')
  # lines(rvec, zSTGDPonly_ecdf1, type = 'l', col = 'red')
  # lines(rvec, rvec , col = 'black',lty=2)
  # 
  # 
  # legend('bottomright', legend = c('Conformal Quantile Regression', 'Quantile regression'), cex = 1,fill = c('blue', 'red'))
  # 
  # QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  # QQ <- c(0.01, QQ, 0.99)
  # quant_est_cqr <- unique(CzST_ecdf1)[2: (length(unique(CzST_ecdf1))-1)]
  # quant_est_qr <- unique(zST_ecdf1)[2: (length(unique(zST_ecdf1))-1)]
  # 
  # quant_est_cqrGDP <- unique(CzSTGDPonly_ecdf1)[2: (length(unique(CzSTGDPonly_ecdf1))-1)]
  # quant_est_qrGDP <- unique(zSTGDPonly_ecdf1)[2: (length(unique(zSTGDPonly_ecdf1))-1)]
  # 
  # # Compute RMSE
  # rmse_cqr  <- sqrt(mean((quant_est_cqr  - QQ)^2))
  # rmse_qr  <- sqrt(mean((quant_est_qr  - QQ)^2))
  # 
  # # Compute MAE
  # mae_cqr  <- mean(abs(quant_est_cqr  - QQ))
  # mae_qr  <- mean(abs(quant_est_qr  - QQ))
  # 
  # # Compare RMSE and MAE
  # list(RMSE_cqr  = rmse_cqr , RMSE_qr  = rmse_qr , MAE_cqr  = mae_cqr , MAE_qr  = mae_qr )