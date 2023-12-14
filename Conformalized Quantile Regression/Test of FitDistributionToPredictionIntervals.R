# Load required libraries
library(sn)
library(moments)

# Function to construct confidence intervals
construct_confidence_intervals <- function(levels) {
  lower_bounds <- numeric(length(levels))
  upper_bounds <- numeric(length(levels))
  
  for (i in seq_along(levels)) {
    lower_bounds[i] <- qst((1 - levels[i]) / 2, true_mu, true_sigma, true_alpha, true_nu)
    upper_bounds[i] <- qst(1 - (1 - levels[i]) / 2, true_mu, true_sigma, true_alpha, true_nu)
  }
  
  return(list(lower_bounds = lower_bounds, upper_bounds = upper_bounds))
}

# Set true parameters
true_mu <- 0
true_sigma <- 1
true_alpha <- 2
true_nu <- 5


# Construct confidence intervals
confidence_levels <- c(0.05, 0.2, 0.4, 0.7,0.9)
confidence_intervals <- construct_confidence_intervals(confidence_levels)

# Test your optimization function
result <- FitDistributionToPredictionIntervals(
  lower_bounds = confidence_intervals$lower_bounds,
  upper_bounds = confidence_intervals$upper_bounds,
  levels = confidence_levels
)

# Display the true and estimated parameters
cat("True Parameters:\n")
cat("mu:", true_mu, "\n")
cat("sigma:", true_sigma, "\n")
cat("alpha:", true_alpha, "\n")
cat("nu:", true_nu, "\n")

cat("\nEstimated Parameters:\n")
cat("mu:", result$lc, "\n")
cat("sigma:", result$sc, "\n")
cat("alpha:", result$sh, "\n")
cat("nu:", result$df, "\n")

