# Generate data
NFCI <- rnorm(10000, 0, 1)
A191RL1Q225SBEA <- 2 * NFCI + rnorm(10000, 0, 1)
data <- data.frame(NFCI, A191RL1Q225SBEA)
alpha <- 0.2

# Randomly shuffle indices
indices <- sample(1:nrow(data))

# Define index ranges for splitting
I1 <- indices[1:3333]
I2 <- indices[3334:6666]
I3 <- indices[6668:10000]

# Create training, test, and validation sets
train <- data[I1,]
test <- data[I2,]
validation <- data[I3,]

# Fit quantile regression on the first split (I1)
fit_lo <- rq(train$A191RL1Q225SBEA ~ train$NFCI, tau = 0.3) 
fit_hi <- rq(train$A191RL1Q225SBEA ~ train$NFCI, tau = alpha + 0.3) 

# Extract quantile functions evaluated on test
X_test <- data.frame(NFCI = test$NFCI)
q_alpha_lo <- predict(fit_lo, newdata = X_test)
q_alpha_hi <- predict(fit_hi, newdata = X_test)

# Initialize a vector for errors
E_i <- rep(NA, length(I2))

# Calculate errors for each point in the test set
for (i in seq_along(I2)) {
  E_i[i] <- max(q_alpha_lo[i] - test$A191RL1Q225SBEA[i], test$A191RL1Q225SBEA[i] - q_alpha_hi[i])
}

# Compute Q(1-??)(E, I2)
quantile_E <- quantile(E_i, pmin(1, pmax(0, (1 - alpha) * (1 + 1/length(I2)))))


# Evaluate adjusted quantile regression lines on validation set
X_validation <- data.frame(NFCI = validation$NFCI)  # Create a data frame with the predictor variable
q_alpha_lo_eval <- predict(fit_lo, newdata = X_validation)
q_alpha_hi_eval <- predict(fit_hi, newdata = X_validation)

# Adjusted quantile regression lines with conformal regression adjustment
q_alpha_lo_adjusted_eval <- q_alpha_lo_eval - quantile_E
q_alpha_hi_adjusted_eval <- q_alpha_hi_eval + quantile_E

# Check if A191RL1Q225SBEA falls between adjusted quantiles
in_quantile_range_adjusted <- (validation$A191RL1Q225SBEA >= q_alpha_lo_adjusted_eval) & (validation$A191RL1Q225SBEA <= q_alpha_hi_adjusted_eval)

# Calculate the percentage of observations within the adjusted quantile range
percentage_in_quantile_range_adjusted <- mean(in_quantile_range_adjusted) * 100

# Display result
cat("Percentage of observations within the adjusted quantile range:", percentage_in_quantile_range_adjusted, "%\n")


