#Using ONLINE code I estimate quantiles from 0.01 to 0.99 of "bio" dataset. I repeat the operation with 20 diff seeds and do a final average

#RESULTS: CQR IMPROVES ALWAYS BUT ONCE.

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

library(quantreg)
library(dplyr)

CASP <- read.csv("CASP.csv")



#Standardizzo
CASP[, -1] <- scale(CASP[, -1])  # scale() defaults to zero mean and unit variance

# Rescale the response variable (assuming it's in the first column)
CASP[, 1] <- CASP[, 1] / mean(abs(CASP[, 1]))

seed_length <- 3
QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(0.01, QQ, 0.99)

# Number of rows in the dataset
norig <- nrow(CASP)
CASPorig <- CASP
coverageCQR <- matrix(NA,seed_length, length(QQ))
coverageQR <- matrix(NA,seed_length, length(QQ))


for (seed in 1:seed_length) {
  
  set.seed(seed)
  
  indices <- sample(1:(norig/50))
  CASP <- CASPorig[indices,]
  n <- nrow(CASP)
  
  # Shuffling indices to randomly select 80% of the data for training
  train_indices <- sample(1:n, size = floor(0.8 * n))
  
  # Create training and test datasets
  train_data <- CASP[train_indices, ]
  test_data <- CASP[-train_indices, ]
  
  full_length <- nrow(train_data)
  test_length = floor(full_length*66/100)
  train2 <- train_data[1:test_length,]
  calibration <- train_data[(test_length + 1) : full_length,]
  
  
  
  # -------------------------------- CQR
  
  

  for (jq in 1:length(QQ)) {  
    
    Y0 <- train2[,1]
    X0 <- train2[,-1]
    
    Y1 <- calibration[,1]
    X1 <- calibration[,-1]
    Y.test <- test_data[,1]
    X.test <- test_data[,-1]
    
    res.cqr <- cqrmod(Y0,X0,Y1,X1,Y.test,X.test,1-QQ[jq])
    
    
    coverageCQR[seed,jq] <- sum(res.cqr$cov.o) / length(res.cqr$cov.o)
    
    
  }
  
  #---------------------- QR
  
  
  Quant <- matrix(NA, nrow(test_data),length(QQ))
  for (jq in 1:length(QQ)) {
    QR <- rq(RMSD ~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9, data = train_data, tau = QQ[jq]) 
    Quant[, jq] <- as.matrix(cbind(1,test_data[,-1])) %*% coef(QR) 
    
    #CALIBRATION OF QR
    
    PitQROOS <- rep(NA, nrow(test_data))
    
    for (i in 1:nrow(test_data)){
      if (test_data[i,1] <= Quant[i,jq]) {
        PitQROOS[i]  <- 1 
      }
      else 
        PitQROOS[i]  <- 0
    }
    
    coverageQR[seed,jq] <- sum(PitQROOS) / length(PitQROOS)
  }



}



# RESULTS EVALUATION

  # Calculate average coverage across all iterations and for each quantile
  average_coverageCQR <- rep(NA,length(QQ))
  average_coverageQR <- rep(NA,length(QQ))
  
  for (jq in 1:length(QQ)) { 
    average_coverageCQR[jq] <- mean(coverageCQR[,jq])
    average_coverageQR[jq] <- mean(coverageQR[,jq])
  }
  
  # Print results
  cat("Average coverage CQR: ", average_coverageCQR, "\n")  
  cat("Average coverage QR: ", average_coverageQR, "\n")
  
  
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
  
  
  
  resultsPitST <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
  
  for (jq in 1:length(QQ)) {
    
    success_rate <- average_coverageQR[jq]

    n <- nrow(test_data)*seed_length
    
    successes <- success_rate*n
    
    # Calculate Wilson score interval
    ci <- wilson_score_interval(successes, n)
    
    # Add to results data frame
    resultsPitST <- rbind(resultsPitST, data.frame(Quantile = QQ[jq], SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
  }
  
  
  # Apply the function to each row and add the result as a new column
  resultsPitST$IsWithinCI <- mapply(is_within_ci, resultsPitST$Quantile, resultsPitST$CI_Lower, resultsPitST$CI_Upper)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitST$IsWithinCI) * 100
  
  # Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitST$SuccessRate > resultsPitST$Quantile & !resultsPitST$IsWithinCI)
  
  # Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitST$SuccessRate < resultsPitST$Quantile & !resultsPitST$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
  print(paste("Number of times SuccessRate is above the Quantile:", count_above))
  print(paste("Number of times SuccessRate is below the Quantile:", count_below))
  print(paste("MAE of QR:", mean(abs(resultsPitST$Quantile-resultsPitST$SuccessRate))))
  
  
  
  
  resultsCPitST <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
  
  for (jq in 1:length(QQ)) {
    
    success_rate <- average_coverageCQR[jq]
    
    n <- nrow(test_data)*seed_length
    
    successes <- success_rate*n
    
    # Calculate Wilson score interval
    ci <- wilson_score_interval(successes, n)
    
    # Add to results data frame
    resultsCPitST <- rbind(resultsCPitST, data.frame(Quantile = QQ[jq], SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
  }
  
  
  
  # Apply the function to each row and add the result as a new column
  resultsCPitST$IsWithinCI <- mapply(is_within_ci, resultsCPitST$Quantile, resultsCPitST$CI_Lower, resultsCPitST$CI_Upper)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsCPitST$IsWithinCI) * 100
  
  # Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsCPitST$SuccessRate > resultsCPitST$Quantile & !resultsCPitST$IsWithinCI)
  
  # Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsCPitST$SuccessRate < resultsCPitST$Quantile & !resultsCPitST$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for CQR :", percentage_within_ci))
  print(paste("Number of times SuccessRate is above the Quantile:", count_above))
  print(paste("Number of times SuccessRate is below the Quantile:", count_below))
  print(paste("MAE of CQR:", mean(abs(resultsCPitST$Quantile-resultsCPitST$SuccessRate))))
  
  
  
  
  
  #66-33 split, 20 seeds, 2% of dataset:
  #qr:  19-1-1 , MAE 0.00508456934686443
  #cqr: 21-0-0 , MAE 0.00429612282071297
  
  #66-33 split, 10 seeds, 2% of dataset:
  #qr:  18-1-2 , MAE 0.00625552953421804
  #cqr: 21-0-0 , MAE 0.00542805100182145
  
  #66-33 split, 3 seeds, 2% of dataset:
  #qr:  18-0-3 , MAE 0.0160048573163328
  #cqr: 18-0-3 , MAE 0.0221597710122301
  
  
  
  #66-33 split, 20 seeds, 10% of dataset:
  #qr:  21-0-0 , MAE 0.00126723913609161
  #cqr: 20-1-0 , MAE 0.00144678636481914
  
  #66-33 split, 10 seeds, 10% of dataset:
  #qr:  20-0-1 , MAE 0.00252406973718452
  #cqr: 21-0-0 , MAE 0.00206088992974239 
  
  #66-33 split, 3 seeds, 10% of dataset:
  #qr:  20-0-1 , MAE 0.00663544106167061
  #cqr: 21-0-0 , MAE 0.00508283459103132
  