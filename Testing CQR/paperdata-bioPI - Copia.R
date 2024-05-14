#calcolo tanti PI per vedere se almeno qui CQR funziona meglio



library(quantreg)
library(quantregForest)



CASP <- read.csv("CASP.csv")
Concrete_data <- read.csv("Concrete_data.csv")


#Standardizzo
CASP[, -1] <- scale(CASP[, -1])  # scale() defaults to zero mean and unit variance
# Rescale the response variable (assuming it's in the first column)
CASP[, 1] <- CASP[, 1] / mean(abs(CASP[, 1]))

#Standardizzo
Concrete_data[, -1] <- scale(Concrete_data[, -1])  # scale() defaults to zero mean and unit variance
# Rescale the response variable (assuming it's in the first column)
Concrete_data[, 1] <- Concrete_data[, 1] / mean(abs(Concrete_data[, 1]))
QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(QQ, 0.99)
CASPorig <- CASP
norig <- nrow(CASPorig)

#----------------------------------- CASP dataset

# Initialize vectors to store the coverage results for each method
coverageCQR_casp <- matrix(NA,5,length(QQ))
coverageQR_casp <- matrix(NA,5,length(QQ))

for (seed in 1:5) {
  set.seed(seed)
  
  indices <- sample(1:(norig/50))
  CASP <- CASPorig[indices,]
  n <- nrow(CASP)
  
  # Number of rows in the dataset
  n <- nrow(CASP)

# Number of rows in the dataset
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


PitCQROOS <- matrix(NA, nrow(test_data),length(QQ))
PitQROOS <- matrix(NA, nrow(test_data),length(QQ))
coverageCQR <- rep(NA,length(QQ))
coverageQR <- rep(NA,length(QQ))

for(jq in 1:length(QQ)) {
   QRlow <- rq(RMSD ~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9, data = train2, tau = (QQ[jq]/2)) 
   Q_low <-  as.matrix(cbind(1,calibration[,-1])) %*% coef(QRlow) 
   
   QRhigh <- rq(RMSD ~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9, data = train2, tau = (1 - (QQ[jq]/2)))
   Q_high <- as.matrix(cbind(1,calibration[,-1])) %*% coef(QRhigh) 

# Fit the quantile random forest model
  # qrf_model <- quantregForest(train2[, c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9")], train2$RMSD)#, alpha = QQ[jq]/2)
  # predictions <- predict(qrf_model, newdata = calibration[,-1], what = c(QQ[jq]/2, 1 - QQ[jq]/2))
  # Q_low <- predictions[,1]
  # Q_high <- predictions[,2]

  
  # Initialize a vector for errors
  E_i <- rep(NA, nrow(calibration))
  
  # Calculate errors for each point in the test set I2
  for (i in 1:length(E_i)) {
    E_i[i] <- max(Q_low[i] - calibration[i,1], calibration[i,1] - Q_high[i])
  }
  
  
  # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
  quantile_E <- quantile(E_i, (1-QQ[jq]) * (1 + 1/nrow(calibration)))
  
  #CQuantup <- predict(qrf_model, newdata = test_data[,-1], what = 1 - QQ[jq]/2) + quantile_E
  #CQuantdown <- predict(qrf_model, newdata = test_data[,-1], what = QQ[jq]/2) - quantile_E
  CQuantup <- as.matrix(cbind(1,test_data[,-1])) %*% coef(QRhigh) + quantile_E
  CQuantdown <- as.matrix(cbind(1,test_data[,-1])) %*% coef(QRlow) - quantile_E



for (i in 1:nrow(test_data)){
  if ( test_data[i,1] >= CQuantdown[i] && test_data[i,1] <= CQuantup[i]) {
      PitCQROOS[i,jq]  <- 1 
  }
  else 
    PitCQROOS[i,jq]  <- 0
}

coverageCQR[jq] <- sum(PitCQROOS[,jq]) / length(PitCQROOS[,jq])


#---------------------- QR


QRlow <- rq(RMSD ~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9, data = train_data, tau = QQ[jq]/2) 
Quantdown <- as.matrix(cbind(1,test_data[,-1])) %*% coef(QRlow) 

QRhigh <- rq(RMSD ~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9, data = train_data, tau = (1 - (QQ[jq]/2))) 
Quantup <- as.matrix(cbind(1,test_data[,-1])) %*% coef(QRhigh) 





for (i in 1:nrow(test_data)){
  if ( test_data[i,1] >= Quantdown[i] && test_data[i,1] <= Quantup[i]) {
    PitQROOS[i,jq]  <- 1 
  }
  else 
    PitQROOS[i,jq]  <- 0
}

coverageQR[jq] <- sum(PitQROOS[,jq]) / length(PitQROOS[,jq])
}

paste("coverage CQR: ", coverageCQR)
paste("coverage QR: ", coverageQR)

coverageCQR_casp[seed,] <- coverageCQR
coverageQR_casp[seed,] <- coverageQR

}



# Calculate average coverage across all iterations and for each quantile
average_coverageCQR <- rep(NA,length(QQ))
average_coverageQR <- rep(NA,length(QQ))

for (jq in 1:length(QQ)) { 
  average_coverageCQR[jq] <- mean(coverageCQR_casp[,jq])
  average_coverageQR[jq] <- mean(coverageQR_casp[,jq])
}

# Print results
cat("Average coverage CQR: ", average_coverageCQR, "\n")  
cat("Average coverage QR: ", average_coverageQR, "\n")


print(paste("MAE of QR:", mean(abs(average_coverageQR - 1 + QQ))))
print(paste("MAE of CQR:", mean(abs(average_coverageCQR - 1 + QQ))))































#----------------------------------- Concrete_data dataset

# Initialize vectors to store the coverage results for each method
coverageCQR_concrete <- numeric(20)
coverageQR_concrete <- numeric(20)

for (seed in 1:20) {
  set.seed(seed)
  
  # Number of rows in the dataset
  n <- nrow(Concrete_data)
  
  # Shuffling indices to randomly select 80% of the data for training
  train_indices <- sample(1:n, size = floor(0.8 * n))
  
  # Create training and test datasets
  train_data <- Concrete_data[train_indices, ]
  test_data <- Concrete_data[-train_indices, ]
  
  
  
  #--------------------------- CQR
  full_length <- nrow(train_data)
  test_length = floor(full_length*66/100)
  train2 <- train_data[1:test_length,]
  calibration <- train_data[(test_length + 1) : full_length,]
  
  QRlow <- rq(as.formula(paste(names(train_data)[ncol(train_data)], "~", paste(names(train_data)[-ncol(train_data)], collapse=" + "))), data = train2, tau = (alpha.sig/2)) 
  Q_low <-  as.matrix(cbind(1,calibration[,-ncol(train_data)])) %*% coef(QRlow) 
  
  QRhigh <- rq(as.formula(paste(names(train_data)[ncol(train_data)], "~", paste(names(train_data)[-ncol(train_data)], collapse=" + "))), data = train2, tau = (1 - (alpha.sig/2)))
  Q_high <- as.matrix(cbind(1,calibration[,-ncol(train_data)])) %*% coef(QRhigh) 
  
  # Fit the quantile random forest model
  # qrf_model <- quantregForest(train2[, c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9")], train2$RMSD)#, alpha = alpha.sig/2)
  # predictions <- predict(qrf_model, newdata = calibration[,-1], what = c(alpha.sig/2, 1 - alpha.sig/2))
  # Q_low <- predictions[,1]
  # Q_high <- predictions[,2]
  # 
  res <- t(apply(cbind(Q_low,Q_high),1,FUN=sort))
  Q_low <- res[,1]
  Q_high <- res[,2]
  
  # Initialize a vector for errors
  E_i <- rep(NA, nrow(calibration))
  
  # Calculate errors for each point in the test set I2
  for (i in 1:length(E_i)) {
    E_i[i] <- max(Q_low[i] - calibration[i,1], calibration[i,1] - Q_high[i])
  }
  
  
  # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
  quantile_E <- quantile(E_i, (1-alpha.sig) * (1 + 1/nrow(calibration)))
  
  #CQuantup <- predict(qrf_model, newdata = test_data[,-1], what = 1 - alpha.sig/2) + quantile_E
  #CQuantdown <- predict(qrf_model, newdata = test_data[,-1], what = alpha.sig/2) - quantile_E
  CQuantup <- as.matrix(cbind(1,test_data[,-ncol(train_data)])) %*% coef(QRhigh) + quantile_E
  CQuantdown <- as.matrix(cbind(1,test_data[,-ncol(train_data)])) %*% coef(QRlow) - quantile_E
  
  
  
  PitCQROOS <- rep(NA, nrow(test_data))
  for (i in 1:nrow(test_data)){
    if ( test_data[i,1] >= CQuantdown[i] && test_data[i,1] <= CQuantup[i]) {
      PitCQROOS[i]  <- 1 
    }
    else 
      PitCQROOS[i]  <- 0
  }
  
  coverageCQR <- sum(PitCQROOS) / length(PitCQROOS)
  
  
  #-------------------- QR
  
  
  QRlow <- rq(as.formula(paste(names(train_data)[ncol(train_data)], "~", paste(names(train_data)[-ncol(train_data)], collapse=" + "))), data = train_data, tau = (alpha.sig/2))
  Quantdown <- as.matrix(cbind(1,test_data[,-ncol(test_data)])) %*% coef(QRlow) 
  
  
  QRhigh <- rq(as.formula(paste(names(train_data)[ncol(train_data)], "~", paste(names(train_data)[-ncol(train_data)], collapse=" + "))), data = train_data, tau = (1-alpha.sig/2))
  Quantup <- as.matrix(cbind(1,test_data[,-ncol(test_data)])) %*% coef(QRhigh) 
  
  
  
  
  #CALIBRATION OF QR
  
  
  PitQROOS <- rep(NA, nrow(test_data))
  for (i in 1:nrow(test_data)){
    if ( test_data[i,ncol(test_data)] >= Quantdown[i,1] && test_data[i,ncol(test_data)] <= Quantup[i,1]) {
      PitQROOS[i]  <- 1 
    }
    else 
      PitQROOS[i]  <- 0
  }
  
  coverageQR <- sum(PitQROOS) / length(PitQROOS)
  
  coverageCQR_concrete[seed] <- coverageCQR
  coverageQR_concrete[seed] <- coverageQR
  
}






# Calculate average coverage across all iterations
average_coverageCQR <- mean(c(coverageCQR_casp,coverageCQR_concrete))
average_coverageQR <- mean(c(coverageQR_casp,coverageQR_concrete))

# Print results
cat("Average coverage CQR: ", average_coverageCQR, "\n") #0.8964778 
cat("Average coverage QR: ", average_coverageQR, "\n")


