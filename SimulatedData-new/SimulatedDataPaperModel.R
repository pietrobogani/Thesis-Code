cqrmod <- function(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig){
  
  #beta.lo <- rq.fit.br(cbind(1,X0),Y0,tau=(0.001))$coefficients
  beta.hi <- rq.fit.br(cbind(1,X0),Y0,tau=(1-alpha.sig + 0.001))$coefficients
  #beta.50 <- rq.fit.br(cbind(1,X0),Y0,tau=0.5)$coefficients
  
  #tq.lo  <- as.matrix(cbind(1,X1))%*%beta.lo
  tq.hi  <- as.matrix(cbind(1,X1))%*%beta.hi
  #tq.50  <- as.matrix(cbind(1,X1))%*%beta.50
  
  #qsr <- t(apply(cbind(tq.lo,tq.50,tq.hi),1,FUN=sort))
  
  q.lo <- rep(-Inf,length(Y1))
  #q.50 <- qsr[,2]
  q.hi <- tq.hi 
  
  Eo.vec <- Em.vec <- Er.vec <- rep(NA,length(Y1))
  for (t in 1:(length(Y1))){
    Eo.vec[t]   <- max(q.lo[t]-Y1[t],Y1[t]-q.hi[t])
    #Em.vec[t]   <- max((q.lo[t]-Y1[t])/(q.50[t]-q.lo[t]),(Y1[t]-q.hi[t])/(q.hi[t]-q.50[t]))
    #Er.vec[t]   <- max((q.lo[t]-Y1[t])/(q.hi[t]-q.lo[t]),(Y1[t]-q.hi[t])/(q.hi[t]-q.lo[t]))
  }
  
  k     <- ceiling((1-alpha.sig)*(1+length(Y1)))
  Q.Eo  <- sort(Eo.vec)[k]
  #Q.Em  <- sort(Em.vec)[k]
  #Q.Er  <- sort(Er.vec)[k]
  
  #tq.test.lo <- as.matrix(cbind(1,X.test))%*%beta.lo
  #tq.test.50 <- as.matrix(cbind(1,X.test))%*%beta.50
  tq.test.hi <- as.matrix(cbind(1,X.test))%*%beta.hi
  
  #qs.test <- t(apply(cbind(tq.test.lo,tq.test.50,tq.test.hi),1,sort))
  
  #q.test.lo <- qs.test[,1]
  #q.test.50 <- qs.test[,2]
  q.test.hi <-tq.test.hi 
  
  # lb.o  <- q.test.lo - Q.Eo
  ub.o  <- q.test.hi + Q.Eo
  # lb.m  <- q.test.lo - Q.Em * (q.test.50-q.test.lo) 
  # ub.m  <- q.test.hi + Q.Em * (q.test.hi-q.test.50)  
  # lb.r  <- q.test.lo - Q.Er * (q.test.hi-q.test.lo)  
  # ub.r  <- q.test.hi + Q.Er * (q.test.hi-q.test.lo) 
  # 
  cov.o <- (Y.test<=ub.o )
  # cov.m <- (Y.test<=ub.m & Y.test>=lb.m)
  # cov.r <- (Y.test<=ub.r & Y.test>=lb.r)
  # 
  # leng.o <- ub.o-lb.o
  # leng.m <- ub.m-lb.m
  # leng.r <- ub.r-lb.r
  
  return(list(cov.o=cov.o,ub.o=ub.o))
  
} #this one only computes qhigh^

cqr <- function(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig){
  
  beta.lo <- rq.fit.br(cbind(1,X0),Y0,tau=(0.001))$coefficients
  beta.hi <- rq.fit.br(cbind(1,X0),Y0,tau=(1-alpha.sig + 0.001))$coefficients
  beta.50 <- rq.fit.br(cbind(1,X0),Y0,tau=0.5)$coefficients
  
  tq.lo  <- as.matrix(cbind(1,X1))%*%beta.lo
  tq.hi  <- as.matrix(cbind(1,X1))%*%beta.hi
  tq.50  <- as.matrix(cbind(1,X1))%*%beta.50
  
  qsr <- t(apply(cbind(tq.lo,tq.50,tq.hi),1,FUN=sort))
  
  q.lo <- qsr[,1]
  q.50 <- qsr[,2]
  q.hi <- qsr[,3]  
  
  Eo.vec <- Em.vec <- Er.vec <- rep(NA,length(Y1))
  for (t in 1:(length(Y1))){
    Eo.vec[t]   <- max(q.lo[t]-Y1[t],Y1[t]-q.hi[t])
    Em.vec[t]   <- max((q.lo[t]-Y1[t])/(q.50[t]-q.lo[t]),(Y1[t]-q.hi[t])/(q.hi[t]-q.50[t]))
    Er.vec[t]   <- max((q.lo[t]-Y1[t])/(q.hi[t]-q.lo[t]),(Y1[t]-q.hi[t])/(q.hi[t]-q.lo[t]))
  }
  
  k     <- ceiling((1-alpha.sig)*(1+length(Y1)))
  Q.Eo  <- sort(Eo.vec)[k]
  Q.Em  <- sort(Em.vec)[k]
  Q.Er  <- sort(Er.vec)[k]
  
  tq.test.lo <- as.matrix(cbind(1,X.test))%*%beta.lo
  tq.test.50 <- as.matrix(cbind(1,X.test))%*%beta.50
  tq.test.hi <- as.matrix(cbind(1,X.test))%*%beta.hi
  
  qs.test <- t(apply(cbind(tq.test.lo,tq.test.50,tq.test.hi),1,sort))
  
  q.test.lo <- qs.test[,1]
  q.test.50 <- qs.test[,2]
  q.test.hi <- qs.test[,3] 
  
  lb.o  <- q.test.lo - Q.Eo
  ub.o  <- q.test.hi + Q.Eo
  lb.m  <- q.test.lo - Q.Em * (q.test.50-q.test.lo) 
  ub.m  <- q.test.hi + Q.Em * (q.test.hi-q.test.50)  
  lb.r  <- q.test.lo - Q.Er * (q.test.hi-q.test.lo)  
  ub.r  <- q.test.hi + Q.Er * (q.test.hi-q.test.lo) 
  
  cov.o <- (Y.test<=ub.o & Y.test>=lb.o)
  cov.m <- (Y.test<=ub.m & Y.test>=lb.m)
  cov.r <- (Y.test<=ub.r & Y.test>=lb.r)
  
  leng.o <- ub.o-lb.o
  leng.m <- ub.m-lb.m
  leng.r <- ub.r-lb.r
  
  return(list(cov.o=cov.o,cov.m=cov.m,cov.r=cov.r,leng.o=leng.o,leng.m=leng.m,leng.r=leng.r))
  
} #this one uses also qlow, as 0.001


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

{
  step <- 4 # Other choice is h = 1
  n1 <- 198 #time length
  n2 <- 100
  
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
  
  theta <- numeric(n1)
  theta[1] <- 0 #MY CHOICE
  h <- numeric(n1)
  h[1] <- 0 #MY CHOICE
  NFCI <- numeric(n1)
  #NFCI[1] <- 0 #MY CHOICE
  
  
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
  NFCI <- numeric(n1)
  NFCI[1] <- -1 #MY CHOICE
  last_error <- rnorm(1, mean = 0, sd = sqrt(fit$sigma2))
  
  # Simulate new values
  for(i in 2:(2*n1)) {
    new_error <- rnorm(1, mean = 0, sd = sqrt(fit$sigma2))
    NFCI[i] <- c + ar_coefficients * NFCI[i-1] + ma_coefficients * last_error + new_error
    last_error <- new_error
  }
  NFCI <- NFCI[(n1+1):(2*n1)] #I implemented a burn-in
  
  #------ Generate omega(t), nu(t)
  omega <- rnorm(n1, 0, sigmatheta)
  nu <- rnorm(n1, 0, sigmah)
  
  #------ Generate theta(t) = phi*theta(t-1) + beta0 + beta1*NFCI(t) + omega(t)
  for(t in 2:n1) {
    theta[t] <- phi * theta[t-1] + beta0 + beta1*NFCI[t] + omega[t]
  }
  
  #------ Compute c11(v) and E[u(t)]
  c11 <- sqrt(0.5 * v) * gamma(0.5 * (v - 1)) / gamma(0.5 * v)
  Expv <- numeric(n1)
  for(t in 1:n1) {
    Expv[t] <- c11 * theta[t]
  }
  
  #------ Generate u(t) sampled from NCT (v, theta(t))
  u <- numeric(n1)
  for(t in 1:n1) {
    #u[t] <- rnct(n = 1, df = v, ncp = theta[t])
    u[t] <- rt(n = 1, df = v, ncp = theta[t])
  }
  
  #------ Generate eps(t) = u(t) - E[u(t)]
  eps <- u - Expv
  
  #------ Generate h(t) = h(t-1) + nu(t)
  for(t in 2:n1) {
    h[t] <- h[t-1] + nu[t]
  }
  
  #------ Generate y(t) = exp(h(t)/2) * eps(t), which is my GDP growth
  y <- exp(h/2) * eps
  
  Yh <- stats::filter(y, rep(1/step, step), sides=1) # moving average of the h past values
  

  Z <- cbind(1, NFCI, y)
  ZGDPonly <- cbind(1, y)
  Z <- as.matrix(Z)
  ZGDPonly <- as.matrix(ZGDPonly)
  
  
  
  #------------ Now I want to generate n2 = 100 times the next number of the time series, namely n1 + 1 = 99
  
  new_error <- rnorm(n2, mean = 0, sd = sqrt(fit$sigma2))
  NFCIn2 <- c + ar_coefficients * NFCI[n1] + ma_coefficients * last_error + new_error # n2 different values of NFCI[n1+1]
  
  omegan2 <- rnorm(n2, 0, sigmatheta)
  nun2 <- rnorm(n2, 0, sigmah)
  
  thetan2 <- phi * theta[n1] + beta0 + beta1*NFCIn2 + omegan2 # n2 different values of theta[n1+1]
  Expvn2 <- c11 * thetan2
  
  un2 <- numeric(n2)
  for(j in 1:n2) {
    un2[j] <- rt(n = 1, df = v, ncp = thetan2[j])
  }
  
  epsn2 <- un2 - Expvn2
  
  hn2 <- h[n1]+nun2
  
  yn2 <- exp(hn2/2) * epsn2
  
  Yh_test = (Yh[n1-2]+Yh[n1-1]+Yh[n1]+yn2)/4 
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
}

n1 <- 198
n2 <- 100
n3 <- 100
n_simul <- n3
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

# QQ <- seq(0.95, 0.05,by = -0.05) #miscoverage rate
# QQ <- c(0.99, QQ)
QQ <- seq(0.05, 0.95,by = 0.05) #miscoverage rate
QQ <- c(QQ,0.99)
# QQ <-0.1

# Setup parallel cluster
cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
clusterExport(cl, varlist=c("run_simulation","cqr","cqrmod")) # Export the simulation function to each cluster node
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
covQRfulltot <- rep(0,length(QQ))
covQRgdptot <- rep(0,length(QQ))
covCQRfulltot <- rep(0,length(QQ))
covCQRgdptot <- rep(0,length(QQ))




for(res in results){
  for(jq in 1:length(QQ)){
    covQRfulltot[jq] <- covQRfulltot[jq] + sum(res[[1]][,jq])
    covQRgdptot[jq] <- covQRgdptot[jq] +  sum(res[[2]][,jq])
    covCQRfulltot[jq] <- covCQRfulltot[jq] + sum(res[[3]][,jq])
    covCQRgdptot[jq] <- covCQRgdptot[jq] + sum(res[[4]][,jq])
  }
  
}



SuccessRateQRfulltot <- rep(NA,length(QQ))
SuccessRateQRgdptot <- rep(NA,length(QQ))
SuccessRateCQRfulltot <- rep(NA,length(QQ))
SuccessRateCQRgdptot <- rep(NA,length(QQ))


for(jq in 1:length(QQ)) {
  SuccessRateQRfulltot[jq] <- covQRfulltot[jq]/(n2*n3)
  SuccessRateQRgdptot[jq] <- covQRgdptot[jq]/(n2*n3)
  SuccessRateCQRfulltot[jq] <- covCQRfulltot[jq]/(n2*n3)
  SuccessRateCQRgdptot[jq] <- covCQRgdptot[jq]/(n2*n3)
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
  
  


#------------------- WILSON SCORES

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


#------------------------- QR full

{
  resultsPitSTQRfull <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
  
  for (jq in 1:length(QQ)) {
    n <- n2*n3
    
    # Calculate success rate
    success_rate <- SuccessRateQRfulltot[jq]
    
    # Calculate Wilson score interval
    ci <- wilson_score_interval(covQRfulltot[jq], n)
    
    # Add to results data frame
    resultsPitSTQRfull <- rbind(resultsPitSTQRfull, data.frame(Quantile = (1-QQ[jq]), SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
  }
  
  
  
  # Apply the function to each row and add the result as a new column
  resultsPitSTQRfull$IsWithinCI <- mapply(is_within_ci, resultsPitSTQRfull$Quantile, resultsPitSTQRfull$CI_Lower, resultsPitSTQRfull$CI_Upper)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTQRfull$IsWithinCI) * 100
  
  # Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTQRfull$SuccessRate > resultsPitSTQRfull$Quantile & !resultsPitSTQRfull$IsWithinCI)
  
  # Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTQRfull$SuccessRate < resultsPitSTQRfull$Quantile & !resultsPitSTQRfull$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for QR full:", percentage_within_ci))
  print(paste("Number of times SuccessRate is above the Quantile:", count_above))
  print(paste("Number of times SuccessRate is below the Quantile:", count_below))
  
}
#cqr senza qlow: 2 - 2 - 16

#------------------------- QR gdp

{
  resultsPitSTQRgdp <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
  
  for (jq in 1:length(QQ)) {
    n <- n2*n3
    
    # Calculate success rate
    success_rate <- SuccessRateQRgdptot[jq]
    
    # Calculate Wilson score interval
    ci <- wilson_score_interval(covQRgdptot[jq], n)
    
    # Add to results data frame
    resultsPitSTQRgdp <- rbind(resultsPitSTQRgdp, data.frame(Quantile = (1-QQ[jq]), SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
  }
  
  
  
  # Apply the function to each row and add the result as a new column
  resultsPitSTQRgdp$IsWithinCI <- mapply(is_within_ci, resultsPitSTQRgdp$Quantile, resultsPitSTQRgdp$CI_Lower, resultsPitSTQRgdp$CI_Upper)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTQRgdp$IsWithinCI) * 100
  
  # Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTQRgdp$SuccessRate > resultsPitSTQRgdp$Quantile & !resultsPitSTQRgdp$IsWithinCI)
  
  # Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTQRgdp$SuccessRate < resultsPitSTQRgdp$Quantile & !resultsPitSTQRgdp$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for QR gdp:", percentage_within_ci))
  print(paste("Number of times SuccessRate is above the Quantile:", count_above))
  print(paste("Number of times SuccessRate is below the Quantile:", count_below))
}
#cqr senza qlow: 3 - 6 - 11

#------------------------- CQR full

{
  
  resultsPitSTCQRfull <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
  
  for (jq in 1:length(QQ)) {
    n <- n2*n3
    
    # Calculate success rate
    success_rate <- SuccessRateCQRfulltot[jq]
    
    # Calculate Wilson score interval
    ci <- wilson_score_interval(covCQRfulltot[jq], n)
    
    # Add to results data frame
    resultsPitSTCQRfull <- rbind(resultsPitSTCQRfull, data.frame(Quantile = (1-QQ[jq]), SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
  }
  
  
  
  # Apply the function to each row and add the result as a new column
  resultsPitSTCQRfull$IsWithinCI <- mapply(is_within_ci, resultsPitSTCQRfull$Quantile, resultsPitSTCQRfull$CI_Lower, resultsPitSTCQRfull$CI_Upper)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTCQRfull$IsWithinCI) * 100
  
  # Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTCQRfull$SuccessRate > resultsPitSTCQRfull$Quantile & !resultsPitSTCQRfull$IsWithinCI)
  
  # Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTCQRfull$SuccessRate < resultsPitSTCQRfull$Quantile & !resultsPitSTCQRfull$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for CQR full:", percentage_within_ci))
  print(paste("Number of times SuccessRate is above the Quantile:", count_above))
  print(paste("Number of times SuccessRate is below the Quantile:", count_below))
  
}
#cqr con qlow: 0 - 20 - 0
#cqr senza qlow: 9 - 2 - 9

#------------------------- CQR gdp

{ resultsPitSTCQRgdp <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
  
  for (jq in 1:length(QQ)) {
    n <- n2*n3
    
    # Calculate success rate
    success_rate <- SuccessRateCQRgdptot[jq]
    
    # Calculate Wilson score interval
    ci <- wilson_score_interval(covCQRgdptot[jq], n)
    
    # Add to results data frame
    resultsPitSTCQRgdp <- rbind(resultsPitSTCQRgdp, data.frame(Quantile = (1-QQ[jq]), SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
  }
  
  
  
  # Apply the function to each row and add the result as a new column
  resultsPitSTCQRgdp$IsWithinCI <- mapply(is_within_ci, resultsPitSTCQRgdp$Quantile, resultsPitSTCQRgdp$CI_Lower, resultsPitSTCQRgdp$CI_Upper)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTCQRgdp$IsWithinCI) * 100
  
  # Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTCQRgdp$SuccessRate > resultsPitSTCQRgdp$Quantile & !resultsPitSTCQRgdp$IsWithinCI)
  
  # Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTCQRgdp$SuccessRate < resultsPitSTCQRgdp$Quantile & !resultsPitSTCQRgdp$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for CQR gdp):", percentage_within_ci))
  print(paste("Number of times SuccessRate is above the Quantile:", count_above))
  print(paste("Number of times SuccessRate is below the Quantile:", count_below))
}
#cqr con qlow: 4 - 12 - 4
#cqr senza qlow: 2 - 9 - 9


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

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