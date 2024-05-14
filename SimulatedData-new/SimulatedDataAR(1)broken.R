#---------------------------------------------------------------------------------------------------------------------------------------------------------------

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



# da qui comincia la simulazione

run_simulation <- function(){
  
  n <- 1001 #3 elements will be lost in the DGP
  n2 <- 100

#-------------------Generate n1 points for the AR(2) model, this is our DGP
  
  phi_ar2 <- 1  # AR coefficients for AR(2)
  Y_ar2 <- numeric(n)
  Y_ar2[1] <- 1
  for (i in 2:n) {
    Y_ar2[i] <- phi_ar2 * Y_ar2[i-1] + rt(n = 1, df = 2) #heavy tails
  }
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])

  Y_test <- phi_ar2 * Y_ar2[n] + rt(n = n2, df = 2) #I create n2 = 100 points for testing
  
  # Prepare the dataset for AR(3) (excluding the first three NA values)
  data_ar2 <- data.frame(I = 1, Y = Y_ar2[-c(1:3)], Y_lag1 = Y_lag1_ar2[-c(1:3)])
  

  
#--------------- QUANTILE REGRESSION 
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)
  
  qrf_model <- quantregForest(x = as.matrix(data_ar2[, 3]), y = as.matrix(data_ar2[,2]))
  

  QuantAR1_OOS <- matrix(NA, length(QQ))

    for (jq in 1:length(QQ)) {  
      
       # QR1 <- rq(Y ~ Y_lag1, data = data_ar2, tau=QQ[jq])
       # QuantAR1_OOS[jq] <- (c(1, Y_ar2[n])) %*% coef(QR1)
      
      QuantAR1_OOS[jq] <- predict(qrf_model, newdata = as.matrix(rep(Y_ar2[n],n2)), what = QQ[jq])[1]
    
      
    }
  
  

  
  
#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  CQuantAR1_OOS <- matrix(NA, length(QQ))
  full_length <- nrow(data_ar2)
  train_indices = sample(1:full_length,full_length*66/100)
  data_ar2_1 <- data_ar2[train_indices,] #train
  data_ar2_2 <- data_ar2[-train_indices,] #calibration
  Q_low <- matrix(NA, nrow(data_ar2_2), length(QQ))
  Q_high <- matrix(NA, nrow(data_ar2_2), length(QQ))
  
  coverageCQRAR1 <- rep(NA,length(QQ))
  qrf_model <- quantregForest(x = as.matrix(data_ar2_1[, 3]), y = as.matrix(data_ar2_1[,2]))
  
  


    for (jq in 1:length(QQ)) {  
      
      #AR(1)
       Q_low[1:nrow(data_ar2_2), jq] <- -Inf 

       #QR1 <- rq(Y ~ Y_lag1, data = data_ar2_1, tau=QQ[jq])
       #Q_high[,jq] <- as.matrix(data_ar2_2[, c(1,3)]) %*% coef(QR1)
       Q_high[,jq] <- predict(qrf_model, newdata = as.matrix(data_ar2_2[, 3]), what = QQ[jq])
       
       # Initialize a vector for errors
       E_i <- rep(NA, nrow(data_ar2_2))
       
       # Calculate errors for each point in the test set I2
       for (i in 1:length(E_i)) {
         E_i[i] <- max(Q_low[i, jq] - data_ar2_2[i,2], data_ar2_2[i,2] - Q_high[i, jq])
       }
       
       # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
       quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(data_ar2_2)))
       
       #CQuantAR1_OOS[jq] <-  (c(1, Y_ar2[n])) %*% coef(QR1) + quantile_E
       CQuantAR1_OOS[jq] <- predict(qrf_model, newdata = as.matrix(rep(Y_ar2[n],n2)), what = QQ[jq])[1] + quantile_E
       
      # Y0 <- data_ar2_1[,2]
      # X0 <- data_ar2_1[,3]
      # 
      # Y1 <- data_ar2_2[,2]
      # X1 <- data_ar2_2[,3]
      # 
      # Y.test <- Y_test
      # X.test <- rep(Y_ar2[n],n2)
      # temp <- cqrmod(Y0,X0,Y1,X1,Y.test,X.test,1-QQ[jq])
      # coverageCQRAR1[jq] <- sum(temp$cov.o) / length(temp$cov.o)
      
    }
  
  
  
  #-------------------------------------------------------------------------------------------
  
  #CALIBRATION OF QR AND CQR
  
  # This function returns the cumulative probability for a given value X, assuming a step function as CDF
  coverageQRAR1 <- rep(NA,length(QQ))

  
  for (jq in 1:length(QQ)) {
  PitQROOSAR1 <- rep(NA, length(Y_test))

  PitCQROOSAR1 <- rep(NA, length(Y_test))

  for (i in 1:length(Y_test)){
    if (Y_test[i] <= QuantAR1_OOS[jq]) {
      PitQROOSAR1[i]  <- 1 
    }
    else 
      PitQROOSAR1[i]  <- 0
  }
  
  coverageQRAR1[jq] <- sum(PitQROOSAR1) / length(PitQROOSAR1)


   for (i in 1:length(Y_test)){
     if (Y_test[i] <= CQuantAR1_OOS[jq]) {
       PitCQROOSAR1[i]  <- 1 
     }
     else 
       PitCQROOSAR1[i]  <- 0
   }
  
  coverageCQRAR1[jq] <- sum(PitCQROOSAR1) / length(PitCQROOSAR1)

  }
  
    
  ret <- list(coverageQRAR1, coverageCQRAR1)
  
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
  library(readxl)
  library(quantreg)
  library(quantregForest)}
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
coveragetotCQRAR1 <- matrix(NA,n3,length(QQ))



index = 1
for(res in results){
  coveragetotQRAR1[index,] <- res[[1]]
  coveragetotCQRAR1[index,] <- res[[2]]
  index <- index + 1
}



# RESULTS EVALUATION

# Calculate average coverage across all iterations and for each quantile
average_coverageQRAR1 <- rep(NA,length(QQ))
average_coverageCQRAR1 <- rep(NA,length(QQ))

for (jq in 1:length(QQ)) { 
  average_coverageQRAR1[jq] <- mean(coveragetotQRAR1[,jq])

  average_coverageCQRAR1[jq] <- mean(coveragetotCQRAR1[,jq])
}


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



#----------------------- QR AR(1)
{
resultsPitSTQRAR1 <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

for (jq in 1:length(QQ)) {
  
  success_rate <- average_coverageQRAR1[jq]
  
  n <- n2*n3
  
  successes <- success_rate*n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsPitSTQRAR1 <- rbind(resultsPitSTQRAR1, data.frame(Quantile = QQ[jq], SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsPitSTQRAR1$IsWithinCI <- mapply(is_within_ci, resultsPitSTQRAR1$Quantile, resultsPitSTQRAR1$CI_Lower, resultsPitSTQRAR1$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsPitSTQRAR1$IsWithinCI) * 100

# Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsPitSTQRAR1$SuccessRate > resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI)

# Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsPitSTQRAR1$SuccessRate < resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))
print(paste("MAE of QR:", mean(abs(resultsPitSTQRAR1$Quantile-resultsPitSTQRAR1$SuccessRate))))
}


#----------------------- CQR AR(1)
{
resultsPitSTCQRAR1 <- data.frame(Quantile = numeric(), SuccessRate = numeric(), CI_Lower = numeric(), CI_Upper = numeric())

for (jq in 1:length(QQ)) {
  
  success_rate <- average_coverageCQRAR1[jq]
  
  n <- n2*n3
  
  successes <- success_rate*n
  
  # Calculate Wilson score interval
  ci <- wilson_score_interval(successes, n)
  
  # Add to results data frame
  resultsPitSTCQRAR1 <- rbind(resultsPitSTCQRAR1, data.frame(Quantile = QQ[jq], SuccessRate = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
}


# Apply the function to each row and add the result as a new column
resultsPitSTCQRAR1$IsWithinCI <- mapply(is_within_ci, resultsPitSTCQRAR1$Quantile, resultsPitSTCQRAR1$CI_Lower, resultsPitSTCQRAR1$CI_Upper)

# Calculate the percentage of quantiles within the CI
percentage_within_ci <- mean(resultsPitSTCQRAR1$IsWithinCI) * 100

# Count the number of times SuccessRate is above the Quantile when IsWithinCI is FALSE
count_above <- sum(resultsPitSTCQRAR1$SuccessRate > resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI)

# Count the number of times SuccessRate is below the Quantile when IsWithinCI is FALSE
count_below <- sum(resultsPitSTCQRAR1$SuccessRate < resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI)

# Print the result
print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
print(paste("Number of times SuccessRate is above the Quantile:", count_above))
print(paste("Number of times SuccessRate is below the Quantile:", count_below))
print(paste("MAE of CQR:", mean(abs(resultsPitSTCQRAR1$Quantile-resultsPitSTCQRAR1$SuccessRate))))
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





