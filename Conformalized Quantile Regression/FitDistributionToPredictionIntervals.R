FitDistributionToPredictionIntervals <- function(lower_bounds, upper_bounds, levels) {
  
  
  {
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
  lower_bounds = confidence_intervals$lower_bounds
  upper_bounds = confidence_intervals$upper_bounds
  levels = confidence_levels
  
  }
  
  
  
  
  # Set bounds and options for optimization
  LB <- c(-20, 0, -30, 1)  # Adjust the bounds as needed
  UB <- c(20, 50, 30, 30)  # Adjust the bounds as needed
  
  # Set initial conditions for optimization
  lc0 <- NULL  # Provide initial value if available
  sc0 <- NULL  # Provide initial value if available
  sh0 <- NULL  # Provide initial value if available
  df0 <- NULL  # Provide initial value if available
  
  # Set initial conditions for optimization (if not provided)
  if(is.null(lc0)) {
    lc0 <- mean(upper_bounds)
    sc0 <- sd(upper_bounds)
    sh0 <- 3
    df0 <- 1  # Default value, adjust as needed
  }
  
  # Initial conditions
  X0 <- c(lc0, sc0, sh0, df0)
  
  # Estimate approximation error for each possible value of the degrees of freedom parameter, 
  # optimizing over the other three continuous-valued parameters.
  par <- matrix(NA, 30, 4)
  ssq <- numeric(30)
  for (df in 1:30) {
    result <- optim(X0, 
                    fn=function(p) sum((upper_bounds - 
                                          get_upper_bounds(lower_bounds, levels, mu=p[1], sigma=p[2], alpha=p[3], nu=df))^2),
                    lower=LB, 
                    upper=UB, 
                    method="L-BFGS-B",
                    control = list(maxit = 10000))
    
    # Print the parameter values during each iteration
    cat("Iteration for df =", df, "Parameters:", result$par, "\n")
    par[df,] <- result$par
    ssq[df] <- result$value
  }
  
  # Find degree of freedom value that provides the best fit, along with the other three parameters.
  X <- numeric(4)
  X[4] <- which.min(ssq)
  X[1:3] <- par[X[4],]
  
  list(lc = X[1], sc = X[2], sh = X[3], df = X[4])
}




get_upper_bounds <- function(lower_bounds, levels, mu, sigma, alpha, nu) {
  # Initialize vector to store upper bounds
  
  mu<-1.410193 
  sigma<- 0.7405623 
  alpha<- 1 
  nu <- 5
  upper_bounds <- numeric(length(levels))
  
  # Calculate the quantile at the specified levels for each lower bound
  for (i in seq_along(lower_bounds)) {
    # Calculate the CDF at the current lower bound
    lower_cdf <- pst(lower_bounds[i], mu = mu, sigma = sigma, alpha = alpha, nu = nu)
    
    # Calculate the quantile at the specified level
    upper_bounds[i] <- qst(levels[i] + lower_cdf, mu = mu, sigma = sigma, alpha = alpha, nu = nu)
  }
  
  return(upper_bounds)
}
