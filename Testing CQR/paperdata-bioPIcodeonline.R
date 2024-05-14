# I calculate centered PI on two datasets from CQR paper using the online version of CQR. Good results, equal to paper


#NOT A LOT OF ANALYSIS TO DO HERE, JUST TO PLAY WITH

#I'll do the same with code found online
library(quantreg)
library(dplyr)
library(lubridate)
library(tidyr)

cqr <- function(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig){
  
  beta.lo <- rq.fit.br(cbind(1,X0),Y0,tau=(alpha.sig/2))$coefficients
  beta.hi <- rq.fit.br(cbind(1,X0),Y0,tau=(1-alpha.sig/2))$coefficients
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


alpha.sig <- 0.1 #so produces 0.9 intervals



CASP <- read.csv("CASP.csv")
Concrete_data <- read.csv("Concrete_data.csv")
#STAR <- read.csv("STAR.csv")


#Standardizzo
CASP[, -1] <- scale(CASP[, -1])  # scale() defaults to zero mean and unit variance
# Rescale the response variable (assuming it's in the first column)
CASP[, 1] <- CASP[, 1] / mean(abs(CASP[, 1]))

#Standardizzo
Concrete_data[, -1] <- scale(Concrete_data[, -1])  # scale() defaults to zero mean and unit variance
# Rescale the response variable (assuming it's in the first column)
Concrete_data[, 1] <- Concrete_data[, 1] / mean(abs(Concrete_data[, 1]))

#Standardizzo
# STAR[, -1] <- scale(STAR[, -1])  # scale() defaults to zero mean and unit variance
# # Rescale the response variable (assuming it's in the first column)
# STAR[, 1] <- STAR[, 1] / mean(abs(STAR[, 1]))

norig <- nrow(CASP)
CASPorig <- CASP
#----------------------------------- CASP dataset

length_seed <- 8
# Initialize vectors to store the coverage results for each method
coverageCQR_casp <- numeric(length_seed)
coverageQR_casp <- numeric(length_seed)

for (seed in 1:length_seed) {
  set.seed(seed)
  
  indices <- sample(1:(norig/10))
  CASP <- CASPorig[indices,]
  n <- nrow(CASP)

# Number of rows in the dataset
n <- nrow(CASP)

# Shuffling indices to randomly select 80% of the data for training
train_indices <- sample(1:n, size = floor(0.8 * n))

# Create training and test datasets
train_data <- CASP[train_indices, ]
test_data <- CASP[-train_indices, ]



#--------------------------- CQR
full_length <- nrow(train_data)
test_length = floor(full_length*66/100)
train2 <- train_data[1:test_length,]
calibration <- train_data[(test_length + 1) : full_length,]


Y0 <- train2[,1]
X0 <- train2[,-1]

Y1 <- calibration[,1]
X1 <- calibration[,-1]
Y.test <- test_data[,1]
X.test <- test_data[,-1]

res.cqr <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)


coverageCQR <- sum(res.cqr$cov.o) / length(res.cqr$cov.o)


#-------------------- QR


QRlow <- rq(RMSD ~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9, data = train_data, tau = (alpha.sig/2))
Quantdown <- as.matrix(cbind(1,test_data[,-1])) %*% coef(QRlow) 


QRhigh <- rq(RMSD ~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9, data = train_data, tau = (1-alpha.sig/2))
Quantup <- as.matrix(cbind(1,test_data[,-1])) %*% coef(QRhigh) 
  



#CALIBRATION OF QR


PitQROOS <- rep(NA, nrow(test_data))
for (i in 1:nrow(test_data)){
  if ( test_data[i,1] >= Quantdown[i,1] && test_data[i,1] <= Quantup[i,1]) {
    PitQROOS[i]  <- 1 
  }
  else 
    PitQROOS[i]  <- 0
}

coverageQR <- sum(PitQROOS) / length(PitQROOS)

coverageCQR_casp[seed] <- coverageCQR
coverageQR_casp[seed] <- coverageQR

}
# Print results
cat("Average coverage CQR: ", mean(coverageCQR_casp), "\n") #0.8964778 
cat("Average coverage QR: ", mean(coverageQR_casp), "\n")









#----------------------------------- Concrete_data dataset

# Initialize vectors to store the coverage results for each method
coverageCQR_concrete <- numeric(20)
coverageQR_concrete <- numeric(20)

# for (seed in 1:20) {
#   set.seed(seed)
#   
#   # Number of rows in the dataset
#   n <- nrow(Concrete_data)
#   
#   # Shuffling indices to randomly select 80% of the data for training
#   train_indices <- sample(1:n, size = floor(0.8 * n))
#   
#   # Create training and test datasets
#   train_data <- Concrete_data[train_indices, ]
#   test_data <- Concrete_data[-train_indices, ]
#   
#   
#   
#   #--------------------------- CQR
#   full_length <- nrow(train_data)
#   test_length = floor(full_length*66/100)
#   train2 <- train_data[1:test_length,]
#   calibration <- train_data[(test_length + 1) : full_length,]
#   
#   
#   Y0 <- train2[,ncol(train2)]
#   X0 <- train2[,-ncol(train2)]
#   
#   Y1 <- calibration[,ncol(train2)]
#   X1 <- calibration[,-ncol(train2)]
#   Y.test <- test_data[,ncol(train2)]
#   X.test <- test_data[,-ncol(train2)]
#   
#   res.cqr <- cqr(Y0,X0,Y1,X1,Y.test,X.test,alpha.sig)
#   
#   
#   coverageCQR <- sum(res.cqr$cov.o) / length(res.cqr$cov.o)
#   
#   
#   #-------------------- QR
#   
#   
#   QRlow <- rq(as.formula(paste(names(train_data)[ncol(train_data)], "~", paste(names(train_data)[-ncol(train_data)], collapse=" + "))), data = train_data, tau = (alpha.sig/2))
#   Quantdown <- as.matrix(cbind(1,test_data[,-ncol(test_data)])) %*% coef(QRlow) 
#   
#   
#   QRhigh <- rq(as.formula(paste(names(train_data)[ncol(train_data)], "~", paste(names(train_data)[-ncol(train_data)], collapse=" + "))), data = train_data, tau = (1-alpha.sig/2))
#   Quantup <- as.matrix(cbind(1,test_data[,-ncol(test_data)])) %*% coef(QRhigh) 
#   
#   
#   
#   
#   #CALIBRATION OF QR
#   
#   
#   PitQROOS <- rep(NA, nrow(test_data))
#   for (i in 1:nrow(test_data)){
#     if ( test_data[i,ncol(test_data)] >= Quantdown[i,1] && test_data[i,ncol(test_data)] <= Quantup[i,1]) {
#       PitQROOS[i]  <- 1 
#     }
#     else 
#       PitQROOS[i]  <- 0
#   }
#   
#   coverageQR <- sum(PitQROOS) / length(PitQROOS)
#   
#   coverageCQR_concrete[seed] <- coverageCQR
#   coverageQR_concrete[seed] <- coverageQR
#   
# }
# 











# Calculate average coverage across all iterations
average_coverageCQR <- mean(c(coverageCQR_casp,coverageCQR_concrete))
average_coverageQR <- mean(c(coverageQR_casp,coverageQR_concrete))

# Print results
cat("Average coverage CQR: ", average_coverageCQR, "\n") #0.8964778 
cat("Average coverage QR: ", average_coverageQR, "\n")


