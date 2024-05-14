#FA QUELLO CHE DEVE. QUI I PREDICTION INTERVALS SONO MEGLIO CALIBRATI DI QUELLI DI QR

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
  
}

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
  
  n <- 101 #3 elements will be lost in the DGP
  n2 <- 100

#-------------------Generate n1 points for the AR(2) model, this is our DGP
  
  phi_ar2 <- c(0.5, -0.2)  # AR coefficients for AR(2)
  Y_ar2 <- numeric(2*n)
  Y_ar2[1] <- 1
  Y_ar2[2] <- -2
  for (i in 3:(2*n)) {
    Y_ar2[i] <- phi_ar2[1] * Y_ar2[i-1] + phi_ar2[2] * Y_ar2[i-2] + rnorm(n=1,0,1) #rt(n = 1, df = 2) #heavy tails
  }
  Y_ar2 <- Y_ar2[(n + 1): (2*n)] #I implemented a burn-in

  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])
  Y_lag3_ar2 <- c(NA, NA, NA, Y_ar2[1:(length(Y_ar2)-3)])
  
  Y_test <- phi_ar2[1] * Y_ar2[n] + phi_ar2[2] * Y_ar2[n-1] + rnorm(n=n2,0,1)#rt(n = n2, df = 2) #I create n2 = 100 points for testing
  
  # Prepare the dataset for AR(3) (excluding the first three NA values)
  data_ar2 <- data.frame(Y = Y_ar2[-c(1:3)], Y_lag1 = Y_lag1_ar2[-c(1:3)], Y_lag2 = Y_lag2_ar2[-c(1:3)], Y_lag3 = Y_lag3_ar2[-c(1:3)])
  

  
#--------------- QUANTILE REGRESSION 
  
  QQ <- 0.7 #miscoverage rate

  
      QR1low <- rq(Y ~ Y_lag1, data = data_ar2, tau=0.001)
      QR1high <- rq(Y ~ Y_lag1, data = data_ar2, tau=(1-QQ+0.001))
      
      QR1cov <- (Y_test <= as.numeric((c(1, Y_ar2[n])) %*% coef(QR1high)) & Y_test >= as.numeric((c(1, Y_ar2[n])) %*% coef(QR1low)))
      
      QR2low <- rq(Y ~ Y_lag1 + Y_lag2, data = data_ar2, tau=0.001)
      QR2high <- rq(Y ~ Y_lag1 + Y_lag2, data = data_ar2, tau=(1-QQ+0.001))
      
      QR2cov <- (Y_test <= as.numeric((c(1, Y_ar2[n], Y_ar2[n-1])) %*% coef(QR2high)) & Y_test >= as.numeric((c(1, Y_ar2[n], Y_ar2[n-1])) %*% coef(QR2low)))
       
  
#--------------- CONFORMALIZED QUANTILE REGRESSION

  full_length <- nrow(data_ar2)
  test_length = floor(full_length*50/100)
  data_ar2_1 <- data_ar2[1:test_length,] #train
  data_ar2_2 <- data_ar2[(test_length + 1) : full_length,] #calibration
  

      #AR(1)
      Y0 <- data_ar2_1[,1]
      X0 <- data_ar2_1[,2]
      
      Y1 <- data_ar2_2[,1]
      X1 <- data_ar2_2[,2]
      
      Y.test <- Y_test
      X.test <- rep(Y_ar2[n],n2)
      
      resAR1 <- cqrmod(Y0,X0,Y1,X1,Y.test,X.test, QQ)
      CQR1cov <- resAR1$cov.o
      
      #AR(2)
      Y0 <- data_ar2_1[,1]
      X0 <- data_ar2_1[,c(2,3)]
      
      Y1 <- data_ar2_2[,1]
      X1 <- data_ar2_2[,c(2,3)]
      
      Y.test <- Y_test
      X.test <- matrix(c(rep(Y_ar2[n], n2),rep(Y_ar2[n-1], n2)),nrow = n2, ncol = 2, byrow = FALSE)          
      
      resAR2 <- cqrmod(Y0,X0,Y1,X1,Y.test,X.test, QQ)
      CQR2cov <- resAR2$cov.o
      
   
  ret <- list(QR1cov, QR2cov, CQR1cov, CQR2cov)
  
  return(ret)
}





n3 = 100 
n_simul <- n3
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

# Setup parallel cluster
cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
clusterExport(cl, varlist=c("run_simulation","cqr","cqrmod")) # Export the simulation function to each cluster node
clusterEvalQ(cl, { 
  library(readxl)
  library(quantreg)}
) # Load required libraries in each cluster node, repeat as necessary for other libraries

# Run simulations in parallel
results <- parLapply(cl, seeds, function(seed) {
  set.seed(seed)
  run_simulation()
})


# Stop the cluster
stopCluster(cl)

# Extract results
covAR1QRtot <- numeric(0)
covAR2QRtot <- numeric(0)
covAR1CQRtot <- numeric(0)
covAR2CQRtot <- numeric(0)




for(res in results){
  covAR1QRtot <- c(covAR1QRtot, na.omit(res[[1]]))
  covAR2QRtot <- c(covAR2QRtot,  na.omit(res[[2]]))
  covAR1CQRtot <- c(covAR1CQRtot, na.omit(res[[3]]))
  covAR2CQRtot <- c(covAR2CQRtot, na.omit(res[[4]]))


}
paste( "coverage AR(1) QR:", sum(covAR1QRtot)/length(covAR1QRtot))
paste( "coverage AR(2) QR:", sum(covAR2QRtot)/length(covAR2QRtot))
paste( "coverage AR(1) CQR:", sum(covAR1CQRtot)/length(covAR1CQRtot))
paste( "coverage AR(2) CQR:", sum(covAR2CQRtot)/length(covAR2CQRtot))

#con quant low = 0.01

# > paste( "coverage AR(1) QR:", sum(covAR1QRtot)/length(covAR1QRtot))
# [1] "coverage AR(1) QR: 0.8851"
# > paste( "coverage AR(2) QR:", sum(covAR2QRtot)/length(covAR2QRtot))
# [1] "coverage AR(2) QR: 0.8727"
# > paste( "coverage AR(1) CQR:", sum(covAR1CQRtot)/length(covAR1CQRtot))
# [1] "coverage AR(1) CQR: 0.9031"
# > paste( "coverage AR(2) CQR:", sum(covAR2CQRtot)/length(covAR2CQRtot))
# [1] "coverage AR(2) CQR: 0.9072"

#senza quant loq, usando "cqrmod":

# > paste( "coverage AR(1) QR:", sum(covAR1QRtot)/length(covAR1QRtot))
# [1] "coverage AR(1) QR: 0.8851"
# > paste( "coverage AR(2) QR:", sum(covAR2QRtot)/length(covAR2QRtot))
# [1] "coverage AR(2) QR: 0.8727"
# > paste( "coverage AR(1) CQR:", sum(covAR1CQRtot)/length(covAR1CQRtot))
# [1] "coverage AR(1) CQR: 0.9032"
# > paste( "coverage AR(2) CQR:", sum(covAR2CQRtot)/length(covAR2CQRtot))
# [1] "coverage AR(2) CQR: 0.9045"



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
  
}
