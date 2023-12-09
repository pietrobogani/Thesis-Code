QuantilesInterpolation <- function(qqTarg, QQ, lc0=NULL, sc0=NULL, sh0=NULL) {
  
  #correttooo!
  
  

  # Set bounds and options for optimization
  LB <- c(-20, 0, -30)
  UB <- c(20, 50, 30)
  
  # Locate target quantiles
  jq50 <- which.min(abs(QQ - 0.50))
  jq25 <- which.min(abs(QQ - 0.25))
  jq75 <- which.min(abs(QQ - 0.75))
  jq05 <- which.min(abs(QQ - 0.05))
  jq95 <- which.min(abs(QQ - 0.95))
  
  # Set initial conditions for optimization (if not provided)
  if(is.null(lc0)) {
    iqn <- qnorm(0.75) - qnorm(0.25)
    lc0 <- qqTarg[jq50]
    sc0 <- (qqTarg[jq75] - qqTarg[jq25]) / iqn
    sh0 <- 0
  }
  
  # Initial conditions
  X0 <- c(lc0, sc0, sh0)
  
  # Indices of target quantiles
  Select <- c(jq05, jq25, jq75, jq95)
  
  # Estimate approximation error for each possible value of the degrees of freedom parameter, 
  # optimizing over the other three continuous-valued parameters.
  par <- matrix(NA, 30, 3)
  ssq <- numeric(30)
  for (df in 1:30) {
    result <- optim(X0, 
                    fn=function(p) sum((qqTarg[Select] - qst(QQ[Select], xi=p[1], omega=p[2], alpha=p[3], nu=df))^2),
                    lower=LB, 
                    upper=UB, 
                    method="L-BFGS-B")
    par[df,] <- result$par
    ssq[df] <- result$value
  }
  
  # Find degree of freedom value that provides the best fit, along with the three other parameters.
  X <- numeric(4)
  X[4] <- which.min(ssq)
  X[1:3] <- par[X[4],]
  
  list(lc = X[1], sc = X[2], sh = X[3], df = X[4])
}
