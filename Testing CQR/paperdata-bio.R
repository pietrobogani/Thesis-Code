#Using MY code I estimate quantiles from 0.01 to 0.99 of "bio" dataset. I repeat the operation with 20 diff seeds and do a final average

# SI CONFERMA CHE CQR FUNZIONA MEGLIO DI QR. PARE ANCHE CHE LA MIA VERSIONE DEL CODICE FUNZIONI MEGLIO DI QUELLA TROVATA ONLINE (DA PRENDERE CON LE PINZE)

library(quantreg)

CASP <- read.csv("CASP.csv")

#Standardizzo
 CASP[, -1] <- scale(CASP[, -1])  # scale() defaults to zero mean and unit variance
 
# Rescale the response variable (assuming it's in the first column)
 CASP[, 1] <- CASP[, 1] / mean(abs(CASP[, 1]))

# Set seed for reproducibility
 seed_length <- 3

# Number of rows in the dataset
norig <- nrow(CASP)
CASPorig <- CASP

QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(0.01, QQ, 0.99)
coverageCQR <- matrix(NA,seed_length, length(QQ))
coverageQR <- matrix(NA,seed_length, length(QQ))


for (seed in 1:seed_length) {
set.seed(seed)

  indices <- sample(1:(norig/20))
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

QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(0.01, QQ, 0.99)


# -------------------------------- CQR


Q_low <- matrix(NA, nrow(calibration),length(QQ))
Q_high <- matrix(NA, nrow(calibration),length(QQ))
CQuant <- matrix(NA, nrow(test_data),length(QQ))

for (jq in 1:length(QQ)) {  
  
Q_low[, jq] <- -Inf 

QR <- rq(RMSD ~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9, data = train2, tau = QQ[jq]) 
Q_high[,jq] <- as.matrix(cbind(1,calibration[,-1])) %*% coef(QR) 

# Initialize a vector for errors
E_i <- rep(NA, nrow(calibration))

# Calculate errors for each point in the test set I2
for (i in 1:length(E_i)) {
  E_i[i] <- max(Q_low[i,jq] - calibration[i,1], calibration[i,1] - Q_high[i,jq])
}


# Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(calibration)))

CQuant[,jq] <- as.matrix((cbind(1, test_data[,-1]))) %*% coef(QR) + quantile_E




#CALIBRATION OF CQR

PitCQROOS <- rep(NA, nrow(test_data))

for (i in 1:nrow(test_data)){
  if (test_data[i,1] <= CQuant[i,jq]) {
    PitCQROOS[i]  <- 1 
  }
  else 
    PitCQROOS[i]  <- 0
}

coverageCQR[seed,jq] <- sum(PitCQROOS) / length(PitCQROOS)
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


#QR
{
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
}

#CQR
{
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
}



#66-33 split, 20 seeds, 2% of dataset:
#qr:  19-1-1 , MAE 0.00508456934686443
#cqr: 21-0-0 , MAE 0.0030575071558678

#66-33 split, 10 seeds, 2% of dataset:
#qr:  18-1-2 , MAE 0.00625552953421804
#cqr: 21-0-0 , MAE 0.00430392922196199 

#66-33 split, 3 seeds, 2% of dataset:
#qr:  18-0-3 , MAE 0.0160048573163328
#cqr: 18-0-3 , MAE 0.0220747679764073



#66-33 split, 20 seeds, 5% of dataset:
#qr:  16-0-5 , MAE 0.00704824287793724
#cqr: 9-0-12 , MAE 0.0089831565814099

#66-33 split, 10 seeds, 5% of dataset:
#qr:  18-0-3 , MAE 0.00667914327302978
#cqr: 18-0-3 , MAE 0.00779164067373677

#66-33 split, 3 seeds, 5% of dataset:
#qr:  15-1-5 , MAE 0.016514867955916
#cqr: 19-1-1 , MAE 0.0108865321965759



#66-33 split, 20 seeds, 10% of dataset:
#qr:  21-0-0 , MAE 0.00126723913609161
#cqr: 21-0-0 , MAE 0.00136091595107989

#66-33 split, 10 seeds, 10% of dataset:
#qr:  20-0-1 , MAE 0.00252406973718452
#cqr: 21-0-0 , MAE 0.00204527712724435 

#66-33 split, 3 seeds, 10% of dataset:
#qr:  20-0-1 , MAE 0.00663544106167061
#cqr: 20-1-0 , MAE 0.00537774308266113 



#66-33 split, 20 seeds, 50% of dataset:
#qr:  21-0-0 , MAE 0.000792956587839625
#cqr: 18-0-3 , MAE 0.00157654139722807

#66-33 split, 10 seeds, 50% of dataset:
#qr:  18-0-0 , MAE 0.000832005664719447
#cqr: 19-0-2 , MAE 0.00187956223381549

#66-33 split, 3 seeds, 50% of dataset:
#qr:  21-0-0 , MAE 0.00168171357762432
#cqr: 21-0-0 , MAE 0.00203679985005156

