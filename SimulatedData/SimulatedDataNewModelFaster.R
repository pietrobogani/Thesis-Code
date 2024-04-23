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


  step <- 4 # Other choice is h = 1
  T <- 3004 #time length
  jtFirstOOS <- 2800 #First index for out-of-sample computations
  
  
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
  
  theta <- numeric(T)
  theta[1] <- 0 #MY CHOICE
  h <- numeric(T)
  h[1] <- 0 #MY CHOICE
  NFCI <- numeric(T)
  #NFCI[1] <- 0 #MY CHOICE
  
  
  file_path <- "DataVulnerabilityAppendix.xls"
  
  # Read the file
  data <- read_excel(file_path)
  data <- data[,1:3]
  NFCI_real <- data$NFCI
  y_real <- data$A191RL1Q225SBEA
  #------ Fit an AR(1) model to NFCI, then generate a new time series using the fitted model
  #fit <- ar(NFCI_real, order.max = 1)
  fit <- Arima(NFCI_real, order=c(1,0,0))  
  
  # Extract  coefficient
  ar_coefficients <- fit$coef["ar1"]  # For AR(1)
  ma_coefficients <- fit$coef["ma1"]  # For MA(1)
  if(is.na(ma_coefficients))
    ma_coefficients <- 0
  c <- fit$coef["intercept"]
  NFCI <- numeric(T)
  NFCI[1] <- -1 #MY CHOICE
  last_error <- rnorm(1, mean = 0, sd = sqrt(fit$sigma2))
  
  # Simulate new values
  for(i in 2:(2*T)) {
    new_error <- rnorm(1, mean = 0, sd = sqrt(fit$sigma2))
    NFCI[i] <- c + ar_coefficients * NFCI[i-1] + ma_coefficients * last_error + new_error
    last_error <- new_error
  }
  NFCI <- NFCI[(T+1):(2*T)] #I implemented a burn-in
  
  #------ Generate ??(t), ??(t)
  omega <- rnorm(T, 0, sigmatheta)
  nu <- rnorm(T, 0, sigmah)
  
  #------ Generate ??(t) = ??*??(t???1) + ??0 + ??1*NFCI(t) + ??(t)
  for(t in 2:T) {
    theta[t] <- phi * theta[t-1] + beta0 + beta1*NFCI[t] + omega[t]
  }
  
  
  #------ Compute c11(v) and E[u(t)]
  c11 <- sqrt(0.5 * v) * gamma(0.5 * (v - 1)) / gamma(0.5 * v)
  Expv <- numeric(T)
  for(t in 1:T) {
    Expv[t] <- c11 * theta[t]
  }
  
  #------ Generate u(t) ??? NCT (??, ??(t))
  u <- numeric(T)
  for(t in 1:T) {
    #u[t] <- rnct(n = 1, df = v, ncp = theta[t])
    u[t] <- rt(n = 1, df = v, ncp = theta[t])
  }
  
  #------ Generate ??(t) = u(t) ??? E[u(t)]
  eps <- u - Expv
  
  #------ Generate h(t) = h(t???1) + ??(t)
  for(t in 2:T) {
    h[t] <- h[t-1] + nu[t]
  }
  
  #------ Generate y(t) = exp(h(t)/2) * ??(t), which is my GDP growth
  y <- exp(h/2) * eps
  
  
  
  
  Yh <- stats::filter(y, rep(1/step, step), sides=1) # moving average of the h past values
  

  Z <- cbind(1, NFCI, y)
  ZGDPonly <- cbind(1, y)
  Z <- as.matrix(Z)
  ZGDPonly <- as.matrix(ZGDPonly)
  
  #-------------------------------------------------------------------------------------------
  
  
  #QUANTILE REGRESSION 
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ, 0.99)

  Quant_OOS <- matrix(NA, nrow(Z), length(QQ))
  QuantGDPonly_OOS <- matrix(NA, nrow(Z), length(QQ))
  Quantunc_OOS <- matrix(NA, nrow(Z), length(QQ))
  
  
  #OUT OF SAMPLE
  for (jt in jtFirstOOS:(length(Yh)- step)) {
    for (jq in 1:length(QQ)) {  
      
      
      #FULL CONDITIONAL
       # QRfull <- rq(Yh[(step+1):jt] ~ Z[1:(jt - step),-1], tau=QQ[jq])
       # Quant_OOS[jt + step, jq] <- Z[jt, ] %*% coef(QRfull)
       # 
      
      #CONDITIONAL, GDP ONLY
      # QR <- rq(Yh[(step+1):jt] ~ ZGDPonly[1:(jt - step),-1], tau=QQ[jq])
      # QuantGDPonly_OOS[jt + step, jq] <- ZGDPonly[jt, ] %*% coef(QR)
    }
  }
  
  #-------------------------------------------------------------------------------------------
  
  
  # CONFORMALIZED QUANTILE REGRESSION
  
  #OUT OF SAMPLE
  CQuant_OOS <- matrix(NA, nrow(Z), length(QQ))
  CQuantGDPonly_OOS <- matrix(NA, nrow(Z), length(QQ))
  Q_low <- matrix(NA, nrow(Z), length(QQ))
  Q_high <- matrix(NA, nrow(Z), length(QQ))
  
  
  for (jt in jtFirstOOS:(length(Yh)- step)) {
    
    full_length <- length(Yh[(step+1):jt])
    test_length = full_length*50/100
    Yh1 <- Yh[(step+1):(step +test_length)]
    Yh2 <- Yh[(test_length +step+1):jt] #Yh1 and Yh2 correctly have same dimension
    Z1 <- Z[1:test_length,]
    Z2 <- Z[(test_length+1):(jt- step),] #Z1 and Z2 correctly have same dimension
    ZGDPonly1 <- ZGDPonly[1:test_length,]
    ZGDPonly2 <- ZGDPonly[(test_length+1):(jt - step),]
    
    for (jq in 1:length(QQ)) {  
      
      #FULL CONDITIONAL
      Q_low[(step+1):(step +length(Yh1)), jq] <- -Inf 
      QR <- rq(Yh1 ~ Z1[,-1], tau= QQ[jq])
      Q_high[(step+1):(step +length(Yh2)), jq] <- as.vector(Z2 %*% coef(QR)) 

      # Initialize a vector for errors
      E_i <- rep(NA, length(Yh2))
      
      # Calculate errors for each point in the test set I2
      for (i in 1:length(E_i)) {
        E_i[i] <- max(Q_low[step +i, jq] - Yh2[i], Yh2[i] - Q_high[step +i, jq])
      }
      
      # Compute Q(QQ[jq])(E, I2) N.B 1 - ?? = QQ[jq]
      quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/length(Yh2)))))
      
      CQuant_OOS[jt+ step,jq] <- Z[jt,] %*% coef(QR) + quantile_E
      
      
      #CONDITIONAL, GDP ONLY
      # Q_low[(step+1):(step +length(Yh1)), jq] <- -Inf 
      # QR <- rq(Yh1 ~ ZGDPonly1[,-1], tau= QQ[jq])
      # Q_high[(step+1):(step +length(Yh2)), jq] <- as.vector(ZGDPonly2 %*% coef(QR)) 
      # 
      # # Initialize a vector for errors
      # E_i <- rep(NA, length(Yh2))
      # 
      # # Calculate errors for each point in the test set I2
      # for (i in 1:length(E_i)) {
      #   E_i[i] <- max(Q_low[step +i, jq] - Yh2[i], Yh2[i] - Q_high[step +i, jq])
      # }
      # 
      # # Compute Q(QQ[jq])(E, I2) N.B 1 - ?? = QQ[jq]
      # quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/length(Yh2)))))
      # 
      # CQuantGDPonly_OOS[jt+ step,jq] <- ZGDPonly[jt,] %*% coef(QR) + quantile_E
    
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
  

  PitST_OOS <- rep(NA, length(Yh))
  PitSTGDPonly_OOS <- rep(NA, length(Yh))
  PitSTunc_OOS <- rep(NA, length(Yh))
  
  
  #OUT OF SAMPLE, FULL CONDITIONAL
  
  for (jt in jtFirstOOS:(length(Yh)- step)){
    
    YhRealized <- Yh[jt + step]
    
     # qqTarg <- Quant_OOS[jt + step, ]
     # PitST_OOS[jt + step] <- cumulative_prob(YhRealized, QQ, qqTarg)
     # 
    # qqTargGDPonly <- QuantGDPonly_OOS[jt + step, ]
    # PitSTGDPonly_OOS[jt + step] <- cumulative_prob(YhRealized, QQ, qqTargGDPonly)
    
  }

  
  
  #-------------------------------------------------------------------------------------------
  #CALIBRATION OF CONFORMAL QUANTILE REGRESSION

  CPitST_OOS <- rep(NA, length(Yh))
  CPitSTGDPonly_OOS <- rep(NA, length(Yh))
  CPitSTunc_OOS <- rep(NA, length(Yh))
  
  
 
  #OUT OF SAMPLE, FULL CONDITIONAL
  
  for (jt in jtFirstOOS:(length(Yh)- step)){
    
    YhRealized <- Yh[jt + step]
    
    qqTarg <- CQuant_OOS[jt + step, ]
    CPitST_OOS[jt + step] <- cumulative_prob(YhRealized, QQ, qqTarg)
    
    # qqTargGDPonly <- CQuantGDPonly_OOS[jt + step, ]
    # CPitSTGDPonly_OOS[jt + step] <- cumulative_prob(YhRealized, QQ, qqTargGDPonly)
    
  }
  
  
  
  
  ret <- list(PitST_OOS, CPitST_OOS, PitSTGDPonly_OOS, CPitSTGDPonly_OOS)
  
  return (ret)
}

n_simul <- 500
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

# Setup parallel cluster
cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
clusterExport(cl, varlist=c("run_simulation")) # Export the simulation function to each cluster node
clusterEvalQ(cl, { 
  library(readxl)
  library(quantreg)
  library(forecast)
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
PitST_OOStot <- numeric(0)
CPitST_OOStot <- numeric(0)
PitSTGDPonly_OOStot <- numeric(0)
CPitSTGDPonly_OOStot <- numeric(0)

for(res in results){
  PitST_OOStot <- c(PitST_OOStot, res[[1]])
  CPitST_OOStot <- c(CPitST_OOStot, res[[2]])
  PitSTGDPonly_OOStot <- c(PitSTGDPonly_OOStot, res[[3]])
  CPitSTGDPonly_OOStot <- c(CPitSTGDPonly_OOStot, res[[4]])
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
#   
  
# #-------------------------------- per caricare dati salvati -----------------------
# 
 # filename <- paste("Sim500_1", ".RData", sep="")
# cat(paste("Loading results from file", filename, "\n"))
# 
 # load(filename)
# 
# #-----------------------------------------------------------------------------------
  
  
  
#-------------------------------------------------- WILSON SCORES
#Posso usare anche PitST_OOStot e simili, infatti sono i valori arrotondati al quantile inferiore della cumulativa di Yh_realized
#Se metto dentro "success" < e non <=, sono a posto


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


#--------------- FULL MODEL, QUANTILE REGRESSION
resultsPitST <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

PitST_OOStot <- na.omit(PitST_OOStot)
for (quantile in c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)) {
  successes <- sum(PitST_OOStot < quantile)
  n <- length(PitST_OOStot)
  
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

# Print the counts
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))
print(paste("Percentage of quantiles within the confidence interval:", percentage_within_ci))



#--------------- FULL MODEL, CONFORMAL QUANTILE REGRESSION  
resultsCPitST <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

CPitST_OOStot <- na.omit(CPitST_OOStot)
n <- length(CPitST_OOStot)

for (quantile in c(0.01, seq(0.05, 0.95, by = 0.05), 0.99)) {
  successes <- sum(CPitST_OOStot < quantile)

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

# Print the counts
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))
print(paste("Percentage of quantiles within the confidence interval:", percentage_within_ci))

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  PITtest_env <- new.env()
  source("PITtest.r",local = PITtest_env)
  rstestboot_env <- new.env()
  source("rstestboot.r",local = rstestboot_env)
  
  # PROBABILITY INTEGRAL TRANSFORM (PIT) FOR QUANTILE REGRESSION, FULL CONDITIONAL MODEL
  rvec <- seq(0, 1, by = 0.001)
  zST_ecdf1 <- PITtest_env$PITtest(PitST_OOStot, rvec)
  
  
  # PROBABILITY INTEGRAL TRANSFORM (PIT) FOR CONFORMAL QUANTILE REGRESSION, FULL CONDITIONAL MODEL
  
  rvec <- seq(0, 1, by = 0.001)
  CzST_ecdf1 <- PITtest_env$PITtest(CPitST_OOStot, rvec)
  
  
  # Plot PIT for quantile regression vs. conformal quantile regression 
  plot(rvec, CzST_ecdf1, type = 'l', col = 'blue', xlab = 'tau', ylab = 'Empirical CDF', main = 'Full Conditional Model')
  lines(rvec, zST_ecdf1, type = 'l', col = 'red')
  lines(rvec, rvec , col = 'black',lty=2)
  abline(h = seq(0, 1, by = 0.05), col = "lightgray", lty = "dotted")
  abline(v = seq(0, 1, by = 0.05), col = "lightgray", lty = "dotted")

  
  legend('bottomright', legend = c('Conformal Quantile Regression', 'Quantile regression'), cex = 1,fill = c('blue', 'red'))
  
  
  
  
  
  
  
  
  # PROBABILITY INTEGRAL TRANSFORM (PIT) FOR QUANTILE REGRESSION, GDP ONLY MODEL
  rvec <- seq(0, 1, by = 0.001)
  zSTGDPonly_ecdf1 <- PITtest_env$PITtest(PitSTGDPonly_OOStot, rvec)
  
  
  # PROBABILITY INTEGRAL TRANSFORM (PIT) FOR  CONFORMAL QUANTILE REGRESSION, GDP ONLY MODEL
  
  rvec <- seq(0, 1, by = 0.001)
  CzSTGDPonly_ecdf1 <- PITtest_env$PITtest(CPitSTGDPonly_OOStot, rvec)
  
  
  # Plot PIT for quantile regression vs. conformal quantile regression 
  plot(rvec, CzSTGDPonly_ecdf1, type = 'l', col = 'blue', xlab = 'tau', ylab = 'Empirical CDF', main = 'GDP ONLY MODEL')
  lines(rvec, zSTGDPonly_ecdf1, type = 'l', col = 'red')
  lines(rvec, rvec , col = 'black',lty=2)

  
  legend('bottomright', legend = c('Conformal Quantile Regression', 'Quantile regression'), cex = 1,fill = c('blue', 'red'))

  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ, 0.99)
  quant_est_cqr <- unique(CzST_ecdf1)[2: (length(unique(CzST_ecdf1))-1)]
  quant_est_qr <- unique(zST_ecdf1)[2: (length(unique(zST_ecdf1))-1)]

  quant_est_cqrGDP <- unique(CzSTGDPonly_ecdf1)[2: (length(unique(CzSTGDPonly_ecdf1))-1)]
  quant_est_qrGDP <- unique(zSTGDPonly_ecdf1)[2: (length(unique(zSTGDPonly_ecdf1))-1)]

  # Compute RMSE
  rmse_cqr  <- sqrt(mean((quant_est_cqr  - QQ)^2))
  rmse_qr  <- sqrt(mean((quant_est_qr  - QQ)^2))
  
  # Compute MAE
  mae_cqr  <- mean(abs(quant_est_cqr  - QQ))
  mae_qr  <- mean(abs(quant_est_qr  - QQ))
  
  # Compare RMSE and MAE
  list(RMSE_cqr  = rmse_cqr , RMSE_qr  = rmse_qr , MAE_cqr  = mae_cqr , MAE_qr  = mae_qr )