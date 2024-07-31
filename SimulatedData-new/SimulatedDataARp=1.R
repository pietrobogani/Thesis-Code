#Si vuole testare cosa accade quando c'è una sola covariata. Punto di racconto tra parte empirica e pratica (dove abbiamo solo
# NFCI come covariata). 

# Purtroppo i risultati non sono quelli sperati! Non è vero che "n1" piccoli "siano come" "p" molto alti (In termini di MAE classico)

#Un aspetto interessante chje si nota in tutte le curve di calibrazione è che CQR tende a fare meglio per i quantili estremi.
#La stima conformal è quasi sempre conservativa e questa è una buona cosa per i quantili medio-alti,
#mentre è invece meno conservativa di QRsui quantili bassi ma anche questa è una buona cosa


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
library(mvtnorm)
library(quantregForest)
source("functions.R")



run_simulation <- function(){
  
  # Set parameters
  n <- 121       # Number of data points for training. 3 will be lost!!
  n2 <- 100      # Number of data points for testing
  p <- 1
  phi_ar2 <- c(0.5, -0.2)   # AR(2) coefficients
  num_exog_vars <-  p     # Parameter to control the number of exogenous variables
  beta <- runif(num_exog_vars, 0.1, 1) # Randomly generated coefficients for exogenous variables
  means <- rnorm(p)
  
  if (p == 1) {
    sigma <- as.matrix(runif(p))
  }  else {
    sigma <- diag(runif(p))
  }
  
  # Initialize variables
  Y_ar2 <- numeric(n)  # Series for AR(2) data
  exog_vars <- matrix(nrow = n + 1, ncol = num_exog_vars)  # Matrix for exogenous variables
  
  exog_vars <- rmvnorm( n = n + 1, mean = means, sigma = sigma)

  # Initialize Y values
  Y_ar2[1:2] <- rnorm(2)  # Initial values for the AR(2) process
  
  # Simulate AR(2) with exogenous variables for training data
  for (t in 3:n) {
    Y_ar2[t] <- phi_ar2[1] * Y_ar2[t - 1] + phi_ar2[2] * Y_ar2[t - 2] +
      sum(beta * exog_vars[t, ]) + rt(n = 1, df = 2)
  }
  
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])
  Y_lag3_ar2 <- c(NA, NA, NA, Y_ar2[1:(length(Y_ar2)-3)])
  
  
  # Prepare dataset for model excluding the first two NA values
  data_ar2 <- data.frame(I = 1, Y = Y_ar2[-c(1:3)], Y_lag1 = Y_lag1_ar2[-c(1:3)], Y_lag2 = Y_lag2_ar2[-c(1:3)], Y_lag3 = Y_lag3_ar2[-c(1:3)])
  for (j in 1:num_exog_vars) {
    data_ar2[paste0("X", j)] <- exog_vars[-c(1:3, n+1), j]
  }
  
  # Extract exogenous variables at time n+1 for testing
  exog_vars_n_plus_1 <- exog_vars[n+1, ]
  
  # Simulate test data using the same exogenous variables at n+1
  Y_test <- phi_ar2[1] * Y_ar2[n] + phi_ar2[2] * Y_ar2[n -1] +
    sum(beta * exog_vars_n_plus_1) + rt(n = n2, df = 2)

  QQ <- seq(0.05, 0.85, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)
  
  #--------------- QUANTILE REGRESSION 
  
  QuantAR2_OOS <- matrix(0, length(QQ))

  { for (jq in 1:length(QQ)) {  

    QR2 <- rq(make_formula(2, num_exog_vars), data = data_ar2, tau=QQ[jq])
    QuantAR2_OOS[jq] <- c(1, Y_ar2[n], Y_ar2[n-1], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR2)
  }
  }
  
  
  #--------------- CONFORMAL QUANTILE REGRESSION
  
  full_length <- nrow(data_ar2)
  test_length = floor(full_length*50/100)
  data_ar2_1 <- data_ar2[1:test_length,] #train
  data_ar2_2 <- data_ar2[(test_length + 1) : full_length,] #calibration
  
  {  

    x0 <- data_ar2_1[,-c(2,5)]
    y0 <- data_ar2_1[,2]
    x1 <- data_ar2_2[,-c(2,5)]
    y1 <- data_ar2_2[,2]
    x_test <- c(1,Y_ar2[n],Y_ar2[n-1],exog_vars_n_plus_1)
    
    CQuantAR2_OOS <- qCQR(make_formula(2, num_exog_vars),x0,y0,x1,y1,x_test,QQ) 
    
  }

  
  # CALIBRATION
  
  coverageQRAR2 <- compute_coverage(Y_test, QuantAR2_OOS, QQ)
  coverageCQRAR2 <- compute_coverage(Y_test, CQuantAR2_OOS, QQ)


  
    
  ret <- list(coverageQRAR2, coverageCQRAR2)
  
  return(ret)
}




n2 = 100
n3 = 1000
n_simul <- n3
seeds <- 1:n_simul # Creates a vector of seeds, one for each simulation

QQ <- seq(0.05, 0.85, by = 0.05) #vector of quantiles I'll do quantile regression on
QQ <- c(0.01, QQ)

source("functions.R")


# Setup parallel cluster
cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
clusterExport(cl, varlist=c("run_simulation")) # Export the simulation function to each cluster node
clusterEvalQ(cl, { 
  library(quantreg)
  library(quantregForest)
  library(mvtnorm)
  source("functions.R")
  
}
) # Load required libraries in each cluster node, repeat as necessary for other libraries

# Run simulations in parallel
results <- parLapply(cl, seeds, function(seed) {
  set.seed(seed)
  run_simulation()
})


# Stop the cluster
stopCluster(cl)

# Extract results
coveragetotQRAR2 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR2 <- matrix(NA,n3,length(QQ))





index = 1
for(res in results){
  coveragetotQRAR2[index,] <- res[[1]]
  coveragetotCQRAR2[index,] <- res[[2]]

  index <- index + 1
}



# RESULTS EVALUATION


# Calculate average coverage across all iterations and for each quantile
average_coverageQRAR2 <- compute_average_coverage(coveragetotQRAR2, QQ)
average_coverageCQRAR2 <- compute_average_coverage(coveragetotCQRAR2, QQ)






#----------------------- QR AR(2)
{
  resultsPitSTQRAR2 <- compute_results(average_coverageQRAR2, n2*n3, QQ)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTQRAR2$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTQRAR2$EmpiricalCoverage > resultsPitSTQRAR2$Quantile & !resultsPitSTQRAR2$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTQRAR2$EmpiricalCoverage < resultsPitSTQRAR2$Quantile & !resultsPitSTQRAR2$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of QR:", mean(abs(resultsPitSTQRAR2$Quantile-resultsPitSTQRAR2$EmpiricalCoverage))))

}


#----------------------- CQR AR(2)
{
  resultsPitSTCQRAR2 <- compute_results(average_coverageCQRAR2, n2*n3, QQ)
  
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTCQRAR2$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTCQRAR2$EmpiricalCoverage > resultsPitSTCQRAR2$Quantile & !resultsPitSTCQRAR2$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTCQRAR2$EmpiricalCoverage < resultsPitSTCQRAR2$Quantile & !resultsPitSTCQRAR2$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of CQR:", mean(abs(resultsPitSTCQRAR2$Quantile-resultsPitSTCQRAR2$EmpiricalCoverage))))

}










#------------- PLOT CALIBRATION CURVES OF QR and CQR and DCP with AR(2) MODEL

{
  x1 <- resultsPitSTCQRAR2$Quantile
  y1 <- resultsPitSTCQRAR2$EmpiricalCoverage
  x2 <- resultsPitSTQRAR2$Quantile
  y2 <- resultsPitSTQRAR2$EmpiricalCoverage

  
  # Create the plot
  plot(x1, y1, type = "n", xlab = "Quantile", ylab = "Empirical Coverage", main = "n1 = 18, p = 1", ylim = c(0,1))
  
  # Add the staircase lines for the first dataset
  for (i in 2:length(x1)) {
    segments(x1[i-1], y1[i-1], x1[i-1], y1[i], col = "blue") # Vertical segment
    segments(x1[i-1], y1[i], x1[i], y1[i], col = "blue")     # Horizontal segment
  }
  
  # Add the points for the first dataset
  points(x1, y1, col = "blue")
  
  # Add the staircase lines for the second dataset
  for (i in 2:length(x2)) {
    segments(x2[i-1], y2[i-1], x2[i-1], y2[i], col = "red") # Vertical segment
    segments(x2[i-1], y2[i], x2[i], y2[i], col = "red")     # Horizontal segment
  }
  
  # Add the points for the second dataset
  points(x2, y2, col = "red")
  

  # Add the diagonal line
  abline(a = 0, b = 1, col = "black", lty = 2)
  
  # Optional: Add a legend
  legend("bottomright", legend = c("CQR AR2", "QR AR2", "DCP AR2"), col = c("blue", "red", "green"), pch = c(16, 17, 18), bty = "n", cex = 1.5, pt.cex = 1.5)
}











 filename <- paste("", ".RData",sep="")
 cat(paste("Saving results to file", filename, "\n"))
 
# Save all the variables to the .RData file
 save(
resultsPitSTDCPAR2,
resultsPitSTCQRAR2,
resultsPitSTQRAR2,   
   file=filename
 )
 


# #-------------------------------- per caricare dati salvati -----------------------
# 
# filename <- paste("Sim500_1", ".RData", sep="")
# cat(paste("Loading results from file", filename, "\n"))
# 
# load(filename)
# 
# #-----------------------------------------------------------------------------------



#Plot the performance for p = 1 and n1 moving
 
 library(ggplot2)
 library(dplyr)
 
 # Define data vectors
 numerosities <- c(198, 178, 158, 138, 118, 98, 78, 58, 48, 38, 28, 18, 15, 12, 9)
 qr_values <- c(0.006478, 0.008654, 0.006674, 0.007256, 0.0100125, 0.0099925, 0.015586, 0.019365, 0.0235905, 0.029447, 0.03406947, 0.05467437, 0.064342778, 0.079756111, 0.104554706)
 cqr_values <- c(0.0048655, 0.010363, 0.005653, 0.005991, 0.0062955, 0.005524, 0.0135575, 0.021302, 0.023358, 0.030311, 0.035098947, 0.045901053, 0.062988333, 0.08124, 0.093165294)
 
 # Create a data frame
 data <- data.frame(
   Num = numerosities,
   QR = qr_values,
   CQR = cqr_values,
   Difference = qr_values - cqr_values  # Calculate the difference
 )
 
 # Melt the data for plotting with ggplot2
 library(tidyr)
 data_long <- pivot_longer(data, cols = c("QR", "CQR", "Difference"), names_to = "Type", values_to = "Value")
 
 # Plotting
 ggplot(data_long, aes(x = Num, y = Value, color = Type, group = Type)) +
   geom_point() +
   geom_line() +
   labs(title = "Number of observations vs. MAE for QR, CQR for 1 covariate",
        x = "Number of observations", y = "MAE", color = "Type") +
   theme_minimal() +
   scale_x_reverse()  # Adjust x-axis to display numerosity in reverse order
 
 # Run this plot code in an R environment
 