#quando predico il 101-esimo, con alpha piccolo 0.01 ar1 fa peggio!!! di parecchio
#quando invece predico i futuri 100, ar(1) non ha più problemi...

#ancora con 101-esimo valore!

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
library(ggplot2)


# Clear workspace 
rm(list = ls())

source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")



run_simulation <- function(){
  
  # Set parameters
  n <- 201      # Number of data points for training. 3 will always be lost!
  n2 <- 100      # Number of data points for testing
  p <- 0.3 * n
  phi_ar2 <- c(0.5, -0.2)   # AR(2) coefficients
  num_exog_vars <-  p     # Parameter to control the number of exogenous variables
  beta <- runif(num_exog_vars, 0.1, 1) # Randomly generated coefficients for exogenous variables
  means <- rnorm(p)
  sigma <- diag(runif(p))
  
  # Initialize variables
  Y_ar2 <- numeric(n)  # Series for AR(2) data
  exog_vars <- matrix(nrow = n + 1, ncol = num_exog_vars)  # Matrix for exogenous variables
  
  exog_vars <- rmvnorm( n = n + 1, mean = means, sigma = sigma)

  # Initialize Y values
  Y_ar2[1:2] <- rnorm(2)  
  
  # Simulate AR(2) with exogenous variables for training data
  for (t in 3:n) {
    Y_ar2[t] <- phi_ar2[1] * Y_ar2[t - 1] + phi_ar2[2] * Y_ar2[t - 2] +
      sum(beta * exog_vars[t, ]) + 0.0001*rt(n = 1, df = 2)
  }
  
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])
  Y_lag3_ar2 <- c(NA, NA, NA, Y_ar2[1:(length(Y_ar2)-3)])
  
  
  # Prepare dataset for model excluding the first three NA values
  data_ar2 <- data.frame(I = 1, Y = Y_ar2[-c(1:3)], Y_lag1 = Y_lag1_ar2[-c(1:3)], Y_lag2 = Y_lag2_ar2[-c(1:3)], Y_lag3 = Y_lag3_ar2[-c(1:3)])
  for (j in 1:num_exog_vars) {
    data_ar2[paste0("X", j)] <- exog_vars[-c(1:3, n+1), j]
  }
  
  # Extract exogenous variables at time n+1 for testing
  exog_vars_n_plus_1 <- exog_vars[n+1, ]
  
  # Simulate test data using the same exogenous variables at n+1
  Y_test <- phi_ar2[1] * Y_ar2[n] + phi_ar2[2] * Y_ar2[n -1] +
    sum(beta * exog_vars_n_plus_1) + 0.0001*rt(n = n2, df = 2)
  
  QQ <- seq(0.05, 0.95, by = 0.05) #vector of quantiles I'll do quantile regression on
  QQ <- c(0.01, QQ)
  
  #--------------- QUANTILE REGRESSION 
  
  QuantAR1_OOS <- matrix(0, length(QQ))
  QuantAR2_OOS <- matrix(NA, length(QQ))
  QuantAR3_OOS <- matrix(NA, length(QQ))

  { for (jq in 1:length(QQ)) {  
       
      QR1 <- rq(make_formula(1, num_exog_vars), data = data_ar2, tau=QQ[jq])     
      QuantAR1_OOS[jq] <- c(1, Y_ar2[n], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR1)
      
      QR2 <- rq(make_formula(2, num_exog_vars), data = data_ar2, tau=QQ[jq])
      QuantAR2_OOS[jq] <- c(1, Y_ar2[n], Y_ar2[n-1], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR2)
      
      QR3 <- rq(make_formula(3, num_exog_vars), data = data_ar2, tau=QQ[jq])
      QuantAR3_OOS[jq] <- c(1, Y_ar2[n], Y_ar2[n-1], Y_ar2[n-2], exog_vars[n+1, 1:num_exog_vars]) %*% coef(QR3)

  }
  }
  
  #--------------- QUANTILE RANDOM FOREST 
  
  {
  # qrf_model1 <- quantregForest(x = (data_ar2[,-c(1,2,4,5)]), y = (data_ar2[,2]))
  # qrf_model2 <- quantregForest(x = (data_ar2[,-c(1,2,5)]), y = (data_ar2[,2]))
  # qrf_model3 <- quantregForest(x = (data_ar2[,-c(1,2)]), y = (data_ar2[,2]))
  
  #for (jq in 1:length(QQ)) {  
    
    #QuantAR1_OOS[jq] <- predict(qrf_model1, newdata = c(Y_ar2[n],exog_vars_n_plus_1), what = QQ[jq])
    #QuantAR2_OOS[jq] <- predict(qrf_model2, newdata = c(Y_ar2[n],Y_ar2[n-1],exog_vars_n_plus_1), what = QQ[jq])
    #QuantAR3_OOS[jq] <- predict(qrf_model3, newdata = c(Y_ar2[n],Y_ar2[n-1],Y_ar2[n-2],exog_vars_n_plus_1), what = QQ[jq])
    
  #}
  }
  
  #--------------- CONFORMAL QUANTILE REGRESSION
  
  full_length <- nrow(data_ar2)
  test_length = floor(full_length*50/100)
  data_ar2_1 <- data_ar2[1:test_length,] #train
  data_ar2_2 <- data_ar2[(test_length + 1) : full_length,] #calibration
  
  {  
    x0 <- data_ar2_1[,-c(2,4,5)]
    y0 <- data_ar2_1[,2]
    x1 <- data_ar2_2[,-c(2,4,5)]
    y1 <- data_ar2_2[,2]
    x_test <- c(1,Y_ar2[n],exog_vars_n_plus_1)
    
    CQuantAR1_OOS <- qCQR(make_formula(1, num_exog_vars),x0,y0,x1,y1,x_test,QQ) 
    
    
    x0 <- data_ar2_1[,-c(2,5)]
    y0 <- data_ar2_1[,2]
    x1 <- data_ar2_2[,-c(2,5)]
    y1 <- data_ar2_2[,2]
    x_test <- c(1,Y_ar2[n],Y_ar2[n-1],exog_vars_n_plus_1)
    
    CQuantAR2_OOS <- qCQR(make_formula(2, num_exog_vars),x0,y0,x1,y1,x_test,QQ) 
    
    x0 <- data_ar2_1[,-2]
    y0 <- data_ar2_1[,2]
    x1 <- data_ar2_2[,-2]
    y1 <- data_ar2_2[,2]
    x_test <- c(1,Y_ar2[n],Y_ar2[n-1],Y_ar2[n-2],exog_vars_n_plus_1)
    
    CQuantAR3_OOS <- qCQR(make_formula(3, num_exog_vars),x0,y0,x1,y1,x_test,QQ) 
  }
  
  #--------------- CONFORMAL QUANTILE RANDOM FOREST
  
  {
  #  
  #  x0 <- data_ar2_1[,-c(1,2,4,5)]
  #  y0 <- data_ar2_1[,2]
  #  x1 <- data_ar2_2[,-c(1,2,4,5)]
  #  y1 <- data_ar2_2[,2]
  #  x_test <- c(Y_ar2[n],exog_vars_n_plus_1)
  #  
  #  CQuantAR1_OOS <- qCQRF(x0,y0,x1,y1,x_test,QQ) 
  #  
  #  x0 <- data_ar2_1[,-c(1,2,5)]
  #  y0 <- data_ar2_1[,2]
  #  x1 <- data_ar2_2[,-c(1,2,5)]
  #  y1 <- data_ar2_2[,2]
  #  x_test <- c(Y_ar2[n],Y_ar2[n-1],exog_vars_n_plus_1)
  #  
  #  CQuantAR2_OOS <- qCQRF(x0,y0,x1,y1,x_test,QQ) 
  #  
  #  x0 <- data_ar2_1[,-c(1,2)]
  #  y0 <- data_ar2_1[,2]
  #  x1 <- data_ar2_2[,-c(1,2)]
  #  y1 <- data_ar2_2[,2]
  #  x_test <- c(Y_ar2[n],Y_ar2[n-1],Y_ar2[n-2],exog_vars_n_plus_1)
  #  
  #  CQuantAR3_OOS <- qCQRF(x0,y0,x1,y1,x_test,QQ) 
   }

  #--------------- DISTRIBUTIONAL CONFORMAL PREDICTION
  
  {
    
    
    # source("temp1.R")
    # source("full.R")
    # funs=rq.funs()
    # 
    # predict.fun=funs$predict.fun
    # train.fun=funs$train.fun
    # 
     DQuantAR2_OOS <- matrix(0, length(QQ))
    # 
    #  y <- data_ar2[,2]
    #  x <- data_ar2[,-c(1,2,5)]
    #  x0 <- c(Y_ar2[n],Y_ar2[n-1],exog_vars[n+1, 1:num_exog_vars])
    # 
    #  for (jq in 1:length(QQ)) {
    #    test = dist.conformal.pred(x,y,x0,train.fun = train.fun, predict.fun = predict.fun, verbose = T, alpha = 1 - QQ[jq])
    #    DQuantAR2_OOS[jq] <- test$up[1,1]
    #    cat("Upper level of CI:", test$up[1, 1], "- Lower level of CI:", test$lo[1, 1], "\n")
    #    cat("Minimum of vector Y:", min(y), "\n")
    #  }
    #  
  } 
  
  #-------------------------------------------------------------------------------------------
  
  # CALIBRATION
  
  coverageQRAR1 <- compute_coverage(Y_test, QuantAR1_OOS, QQ)
  coverageQRAR2 <- compute_coverage(Y_test, QuantAR2_OOS, QQ)
  coverageQRAR3 <- compute_coverage(Y_test, QuantAR3_OOS, QQ)
  
  coverageCQRAR1 <- compute_coverage(Y_test, CQuantAR1_OOS, QQ)
  coverageCQRAR2 <- compute_coverage(Y_test, CQuantAR2_OOS, QQ)
  coverageCQRAR3 <- compute_coverage(Y_test, CQuantAR3_OOS, QQ)
  
  coverageDCPAR2 <- compute_coverage(Y_test, DQuantAR2_OOS, QQ)
  
  #------------------- Calculation of CRPS
  
   # observed <- Y_test[1]
   # quantile_values <- QuantAR1_OOS
   # quantile_levels <- QQ
   # mean(abs(coverageQRAR1-QQ))
   # mean(abs(coverageCQRAR1-QQ))
  
  CRPSQRAR1 <- calculate_crps_from_quantiles(-9, QuantAR1_OOS, QQ)
  CRPSQRAR2 <- calculate_crps_from_quantiles(Y_test, QuantAR2_OOS, QQ)
  CRPSQRAR3 <- calculate_crps_from_quantiles(Y_test, QuantAR3_OOS, QQ)
  
  CRPSCQRAR1 <- calculate_crps_from_quantiles(-9, CQuantAR1_OOS, QQ)
  CRPSCQRAR2 <- calculate_crps_from_quantiles(Y_test, CQuantAR2_OOS, QQ)
  CRPSCQRAR3 <- calculate_crps_from_quantiles(Y_test, CQuantAR3_OOS, QQ)
  
  CRPSDCPAR2 <- calculate_crps_from_quantiles(Y_test, DQuantAR2_OOS, QQ)
  
    
  
  # 
  # plot(sort(QuantAR1_OOS[,1]),QQ, type = "s", xlim = c(min(QuantAR1_OOS[1,], Y_test[1]) - 1, max(QuantAR1_OOS[1,], Y_test[1]) + 1), ylim = c(-0.1,1.1),
  #      ylab = "Quantile Levels", xlab = "Values", main = "Heaviside Function and Quantiles")
  # 
  # # Adding a horizontal line for the first observed value
  # abline(v = Y_test[1], col = "red", lty = 2)
  # 
  # # Adding the Heaviside function
  # # This function is 1 if quantile value >= observed value (Y_test[1]), otherwise 0
  # heaviside_values <- as.numeric(QuantAR1_OOS[,1] >= Y_test[1])
  # points(QuantAR1_OOS[,1], heaviside_values, col = "blue", pch = 19)
  # 
  # 
  # plot(sort(CQuantAR1_OOS),QQ, type = "s", xlim = c(min(CQuantAR1_OOS, Y_test[1]) - 1, max(CQuantAR1_OOS, Y_test[1]) + 1), ylim = c(-0.1,1.1),
  #      ylab = "Quantile Levels", xlab = "Values", main = "Heaviside Function and Quantiles")
  # 
  # # Adding a horizontal line for the first observed value
  # abline(v = Y_test[1], col = "red", lty = 2)
  # 
  # # Adding the Heaviside function
  # # This function is 1 if quantile value >= observed value (Y_test[1]), otherwise 0
  # heaviside_valuesC <- as.numeric(CQuantAR1_OOS >= Y_test[1])
  # points(CQuantAR1_OOS, heaviside_valuesC, col = "blue", pch = 19)
  # 
  # 
  
  
  # #ANALYZE BOXPLOTS OF AR1,AR2,AR3 ON QUANTILE VALUES! NOT LEVELS
  # data_df <- data.frame(
  #   Value = c(QuantAR1_OOS, QuantAR2_OOS, QuantAR3_OOS, CQuantAR1_OOS, CQuantAR2_OOS, CQuantAR3_OOS),
  #   Group = factor(rep(c("QuantAR1_OOS", "QuantAR2_OOS", "QuantAR3_OOS", "CQuantAR1_OOS", "CQuantAR2_OOS", "CQuantAR3_OOS"), times = c(length(QuantAR1_OOS), length(QuantAR2_OOS), length(QuantAR3_OOS), length(CQuantAR1_OOS), length(CQuantAR2_OOS), length(CQuantAR3_OOS))))
  # )
  # 
  # # Simplify the group labels
  # data_df$Group <- factor(data_df$Group, levels = c("QuantAR1_OOS", "QuantAR2_OOS", "QuantAR3_OOS", "CQuantAR1_OOS", "CQuantAR2_OOS", "CQuantAR3_OOS"),
  #                         labels = c("AR1", "AR2", "AR3", "CAR1", "CAR2", "CAR3"))
  # 
  # # Create the boxplot
  # ggplot(data_df, aes(x = Group, y = Value, fill = Group)) +
  #   geom_boxplot() +
  #   labs(title = "Comparison of Multiple Vectors", y = "Values", x = "Group") +
  #   theme_minimal() +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
  #     axis.title = element_text(size = 12),
  #     legend.title = element_text(size = 12),
  #     legend.text = element_text(size = 10)
  #   ) +
  #   scale_fill_manual(values = c("AR1" = "lightblue", "AR2" = "lightgreen", "AR3" = "lightcoral", "CAR1" = "lightpink", "CAR2" = "lightgoldenrod", "CAR3" = "lightgray"))
  # 
  # 
  # 
  
  ret <- list(coverageQRAR1, coverageQRAR2, coverageQRAR3, coverageCQRAR1, coverageCQRAR2, coverageCQRAR3, coverageDCPAR2, 
              CRPSQRAR1, CRPSQRAR2, CRPSQRAR3, CRPSCQRAR1, CRPSCQRAR2, CRPSCQRAR3, CRPSDCPAR2)
  
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
  library(quantregForest)
  library(mvtnorm)
  source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")
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
coveragetotQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotQRAR2 <- matrix(NA,n3,length(QQ))
coveragetotQRAR3 <- matrix(NA,n3,length(QQ))

coveragetotCQRAR1 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR2 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR3 <- matrix(NA,n3,length(QQ))

coveragetotDCPAR2 <- matrix(NA,n3,length(QQ))

CRPStotQRAR1 <- rep(NA,n3)
CRPStotQRAR2 <- rep(NA,n3)
CRPStotQRAR3 <- rep(NA,n3)

CRPStotCQRAR1 <- rep(NA,n3)
CRPStotCQRAR2 <- rep(NA,n3)
CRPStotCQRAR3 <- rep(NA,n3)

CRPStotDCPAR2 <- rep(NA,n3)




index = 1
for(res in results){
  coveragetotQRAR1[index,] <- res[[1]]
  coveragetotQRAR2[index,] <- res[[2]]
  coveragetotQRAR3[index,] <- res[[3]]
  
  coveragetotCQRAR1[index,] <- res[[4]]
  coveragetotCQRAR2[index,] <- res[[5]]
  coveragetotCQRAR3[index,] <- res[[6]]
  
  coveragetotDCPAR2[index,] <- res[[7]]
  
  CRPStotQRAR1[index] <- res[[8]]
  CRPStotQRAR2[index] <- res[[9]]
  CRPStotQRAR3[index] <- res[[10]]
  
  CRPStotCQRAR1[index] <- res[[11]]
  CRPStotCQRAR2[index] <- res[[12]]
  CRPStotCQRAR3[index] <- res[[13]]
  
  CRPStotDCPAR2[index] <- res[[14]]
  
  
  index <- index + 1
}


# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)

# Function to randomly select column indices
select_random_columns <- function(mat, num_cols = 5) {
  sample(ncol(mat), num_cols)
}

# Select 5 random columns
set.seed(123)  # Setting seed for reproducibility
selected_columns <- select_random_columns(coveragetotQRAR1, 5)

# Function to extract and reshape data for a given column
extract_and_reshape <- function(df, col_index, type_label) {
  df_selected <- df[, col_index, drop = FALSE]
  colnames(df_selected) <- paste0("V", col_index)
  df_selected %>%
    mutate(Row = 1:nrow(df_selected)) %>%
    pivot_longer(cols = -Row, names_to = "Column", values_to = "Value") %>%
    mutate(Type = type_label)
}

# List of matrices and their labels
matrices <- list(coveragetotQRAR1, coveragetotQRAR2, coveragetotQRAR3, 
                 coveragetotCQRAR1, coveragetotCQRAR2, coveragetotCQRAR3)
labels <- c("QRAR1", "QRAR2", "QRAR3", "CQRAR1", "CQRAR2", "CQRAR3")

# Loop over selected columns and create plots
plots <- list()
for (col_index in selected_columns) {
  combined_long <- bind_rows(
    lapply(seq_along(matrices), function(i) {
      extract_and_reshape(as.data.frame(matrices[[i]]), col_index, labels[i])
    })
  )
  
  # Ensure the columns are ordered correctly
  combined_long$Column <- factor(combined_long$Column, levels = paste0("V", col_index))
  
  # Create an interleaved column for plotting
  combined_long <- combined_long %>%
    mutate(Interaction = paste0(Column, "_", Type)) %>%
    arrange(Column, Type)
  
  # Plot the data
  plot <- ggplot(combined_long, aes(x = Interaction, y = Value, fill = Type)) +
    geom_boxplot() +
    labs(x = 'Column and Type', y = 'Value', title = paste("Boxplot for Quantile", QQ[col_index])) +
    scale_x_discrete(labels = function(x) sub("_.*", "", x)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  plots[[paste0("Plot_", col_index)]] <- plot
}

# Print the plots
for (plot_name in names(plots)) {
  print(plots[[plot_name]])
}


# RESULTS EVALUATION


# Calculate average coverage across all iterations and for each quantile
average_coverageQRAR1 <- compute_average_coverage(coveragetotQRAR1, QQ) #vector of length = length(QQ)
average_coverageQRAR2 <- compute_average_coverage(coveragetotQRAR2, QQ)
average_coverageQRAR3 <- compute_average_coverage(coveragetotQRAR3, QQ)

average_coverageCQRAR1 <- compute_average_coverage(coveragetotCQRAR1, QQ)
average_coverageCQRAR2 <- compute_average_coverage(coveragetotCQRAR2, QQ)
average_coverageCQRAR3 <- compute_average_coverage(coveragetotCQRAR3, QQ)

average_coverageDCPAR2 <- compute_average_coverage(coveragetotDCPAR2, QQ)




#----------------------- QR AR(1)
{
  resultsPitSTQRAR1 <- compute_results(average_coverageQRAR1, n2*n3, QQ)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTQRAR1$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTQRAR1$EmpiricalCoverage > resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTQRAR1$EmpiricalCoverage < resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of QR:", mean(abs(resultsPitSTQRAR1$Quantile-resultsPitSTQRAR1$EmpiricalCoverage))))
  print(paste("CRPS of QR:", mean(CRPStotQRAR1)))
}


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
  print(paste("CRPS of QR:", mean(CRPStotQRAR2)))
  
}


#----------------------- QR AR(3)
{
  resultsPitSTQRAR3 <- compute_results(average_coverageQRAR3, n2*n3, QQ)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTQRAR3$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTQRAR3$EmpiricalCoverage > resultsPitSTQRAR3$Quantile & !resultsPitSTQRAR3$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTQRAR3$EmpiricalCoverage < resultsPitSTQRAR3$Quantile & !resultsPitSTQRAR3$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for QR:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of QR:", mean(abs(resultsPitSTQRAR3$Quantile-resultsPitSTQRAR3$EmpiricalCoverage))))
  print(paste("CRPS of QR:", mean(CRPStotQRAR3)))
  
}


#----------------------- CQR AR(1)
{
  resultsPitSTCQRAR1 <- compute_results(average_coverageCQRAR1, n2*n3, QQ)
  
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTCQRAR1$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTCQRAR1$EmpiricalCoverage > resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTCQRAR1$EmpiricalCoverage < resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of CQR:", mean(abs(resultsPitSTCQRAR1$Quantile-resultsPitSTCQRAR1$EmpiricalCoverage))))
  print(paste("CRPS of CQR:", mean(CRPStotCQRAR1)))
  
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
  print(paste("CRPS of CQR:", mean(CRPStotCQRAR2)))
  
}


#----------------------- CQR AR(3)
{
  resultsPitSTCQRAR3 <- compute_results(average_coverageCQRAR3, n2*n3, QQ)
  
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTCQRAR3$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTCQRAR3$EmpiricalCoverage > resultsPitSTCQRAR3$Quantile & !resultsPitSTCQRAR3$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTCQRAR3$EmpiricalCoverage < resultsPitSTCQRAR3$Quantile & !resultsPitSTCQRAR3$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for CQR:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of CQR:", mean(abs(resultsPitSTCQRAR3$Quantile-resultsPitSTCQRAR3$EmpiricalCoverage))))
  print(paste("CRPS of CQR:", mean(CRPStotCQRAR3)))
  
}



#----------------------- DCP AR(2)

{
  resultsPitSTDCPAR2 <- compute_results(average_coverageDCPAR2, n2*n3, QQ)
  
  # Calculate the percentage of quantiles within the CI
  percentage_within_ci <- mean(resultsPitSTDCPAR2$IsWithinCI) * 100
  
  # Count the number of times EmpiricalCoverage is above the Quantile when IsWithinCI is FALSE
  count_above <- sum(resultsPitSTDCPAR2$EmpiricalCoverage > resultsPitSTDCPAR2$Quantile & !resultsPitSTDCPAR2$IsWithinCI)
  
  # Count the number of times EmpiricalCoverage is below the Quantile when IsWithinCI is FALSE
  count_below <- sum(resultsPitSTDCPAR2$EmpiricalCoverage < resultsPitSTDCPAR2$Quantile & !resultsPitSTDCPAR2$IsWithinCI)
  
  # Print the result
  print(paste("Percentage of quantiles within the confidence interval for DCP:", percentage_within_ci))
  print(paste("Number of times EmpiricalCoverage is above the Quantile:", count_above))
  print(paste("Number of times EmpiricalCoverage is below the Quantile:", count_below))
  print(paste("MAE of DCP:", mean(abs(resultsPitSTDCPAR2$Quantile-resultsPitSTDCPAR2$EmpiricalCoverage))))
  print(paste("CRPS of DCP:", mean(CRPStotDCPAR2)))
  
}









#------------- PLOT CALIBRATION CURVES OF QR and CQR with AR(1) MODEL

{
  x1 <- resultsPitSTCQRAR1$Quantile
  y1 <- resultsPitSTCQRAR1$EmpiricalCoverage
  x2 <- resultsPitSTQRAR1$Quantile
  y2 <- resultsPitSTQRAR1$EmpiricalCoverage
  
  # Create the plot
  plot(x1, y1, type = "n", xlab = "Quantile", ylab = "Empirical Coverage", main = "Staircase Plots with Diagonal Line")
  
  # Add the staircase lines for the first dataset
  for (i in 2:length(x1)) {
    segments(x1[i-1], y1[i-1], x1[i-1], y1[i], col = "blue") # Vertical segment
    segments(x1[i-1], y1[i], x1[i], y1[i], col = "blue")     # Horizontal segment
  }
  
  # Add the points for the first dataset
  points(x1, y1)
  
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
  legend("bottomright", legend = c("CQR AR1", "QR AR1"), col = c("blue", "red"), pch = c(16, 17), bty = "n", cex = 1.5, pt.cex = 1.5)
}


#------------- PLOT CALIBRATION CURVES OF QR and CQR and DCP with AR(2) MODEL

{
  x1 <- resultsPitSTCQRAR2$Quantile
  y1 <- resultsPitSTCQRAR2$EmpiricalCoverage
  x2 <- resultsPitSTQRAR2$Quantile
  y2 <- resultsPitSTQRAR2$EmpiricalCoverage
  

  
  

  
  # Create the plot
  plot(x1, y1, type = "n", xlab = "Quantile", ylab = "Empirical Coverage", main = "n1 = 998, p/n = [0.1,0.2,0.3,0.4]", ylim = c(0,1))
  
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
  
  # Add the staircase lines for the third dataset
  # for (i in 2:length(x3)) {
  #   segments(x3[i-1], y3[i-1], x3[i-1], y3[i], col = "green") # Vertical segment
  #   segments(x3[i-1], y3[i], x3[i], y3[i], col = "green")     # Horizontal segment
  # }
  
  # Add the points for the third dataset
  # points(x3, y3, col = "green")
  
  # Add the diagonal line
  abline(a = 0, b = 1, col = "black", lty = 2)
  
  # Optional: Add a legend
  legend("bottomright", legend = c("CQR AR2", "QR AR2", "DCP AR2"), col = c("blue", "red", "green"), pch = c(16, 17, 18), bty = "n", cex = 1.5, pt.cex = 1.5)
}


#------------- PLOT CALIBRATION CURVES OF QR and CQR with AR(3) MODEL
{
  x1 <- resultsPitSTCQRAR3$Quantile
  y1 <- resultsPitSTCQRAR3$EmpiricalCoverage
  x2 <- resultsPitSTQRAR3$Quantile
  y2 <- resultsPitSTQRAR3$EmpiricalCoverage
  
  # Create the plot
  plot(x1, y1, type = "n", xlab = "Quantile", ylab = "Empirical Coverage", main = "Staircase Plots with Diagonal Line")
  
  # Add the staircase lines for the first dataset
  for (i in 2:length(x1)) {
    segments(x1[i-1], y1[i-1], x1[i-1], y1[i], col = "blue") # Vertical segment
    segments(x1[i-1], y1[i], x1[i], y1[i], col = "blue")     # Horizontal segment
  }
  
  # Add the points for the first dataset
  points(x1, y1)
  
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
  legend("bottomright", legend = c("CQR AR3", "QR AR3"), col = c("blue", "red"), pch = c(16, 17), bty = "n", cex = 1.5, pt.cex = 1.5)
}









#  filename <- paste("pdivison=0,3cp98", ".RData",sep="")
#  cat(paste("Saving results to file", filename, "\n"))
#  
# # Save all the variables to the .RData file
#  save(
# resultsPitSTDCPAR2,
# resultsPitSTCQRAR2,
# resultsPitSTQRAR2, 
# coveragetotDCPAR2, #only for some this is added
#    file=filename
#  )
#  


# #-------------------------------- per caricare dati salvati -----------------------
# 
# filename <- paste("Sim500_1", ".RData", sep="")
# cat(paste("Loading results from file", filename, "\n"))
# 
# load(filename)
# 
# #-----------------------------------------------------------------------------------


 #------------ PLOT TOGETHER 0.1,0.2,0.3,0.4 on AR(2) (sfumato)
 
 

 # Define colors
 cqr_colors <- colorRampPalette(c("#ADD8E6", "#0000FF"))(4)
 qr_colors <- colorRampPalette(c("#F08080", "#FF0000"))(4)
 
 cqr_colors <- c("#E0F7FA", "#81D4FA", "#0288D1", "#0000FF")  # Very light blue to bright blue
 qr_colors <- c("#FFEBEE", "#FFCDD2", "#E57373", "#D32F2F")
 
 x1 <- resultsPitSTCQR0.1$Quantile
 y1 <- resultsPitSTCQR0.1$EmpiricalCoverage
 x2 <- resultsPitSTQR0.1$Quantile
 y2 <- resultsPitSTQR0.1$EmpiricalCoverage
 
 x3 <- resultsPitSTCQR0.2$Quantile
 y3 <- resultsPitSTCQR0.2$EmpiricalCoverage
 x4 <- resultsPitSTQR0.2$Quantile
 y4 <- resultsPitSTQR0.2$EmpiricalCoverage
 
 x5 <- resultsPitSTCQR0.3$Quantile
 y5 <- resultsPitSTCQR0.3$EmpiricalCoverage
 x6 <- resultsPitSTQR0.3$Quantile
 y6 <- resultsPitSTQR0.3$EmpiricalCoverage
 
 x7 <- resultsPitSTCQR0.4$Quantile
 y7 <- resultsPitSTCQR0.4$EmpiricalCoverage
 x8 <- resultsPitSTQR0.4$Quantile
 y8 <- resultsPitSTQR0.4$EmpiricalCoverage
 
 # Create the plot
 plot(x2, y2, type = "n", xlab = "Quantile", ylab = "Empirical Coverage", main = "n1 = 198, p/n = [0.1,0.2,0.3,0.4]", ylim = c(0, 1))
 
 # Function to add staircase lines
 add_staircase <- function(x, y, color) {
   for (i in 2:length(x)) {
     segments(x[i-1], y[i-1], x[i-1], y[i], col = color) # Vertical segment
     segments(x[i-1], y[i], x[i], y[i], col = color)     # Horizontal segment
   }
   points(x, y, col = color)
 }
 
 # Add the staircases
 #add_staircase(x1, y1, cqr_colors[1])
 add_staircase(x2, y2, qr_colors[1])
 #add_staircase(x3, y3, cqr_colors[2])
 add_staircase(x4, y4, qr_colors[2])
 #add_staircase(x5, y5, cqr_colors[3])
 add_staircase(x6, y6, qr_colors[3])
 add_staircase(x7, y7, cqr_colors[4])
 add_staircase(x8, y8, qr_colors[4])
 
 # Add the diagonal line
 abline(a = 0, b = 1, col = "black", lty = 2)
 
 
 
 
 
 # Define colors
 cqr_colors <- c("#E0F7FA", "#81D4FA", "#0288D1", "#0000FF")  # Very light blue to bright blue
 qr_colors <- c("#FFEBEE", "#FFCDD2", "#E57373", "#D32F2F")
 
 # Prepare data
 results_list <- list(
   data.frame(Quantile = resultsPitSTQR0.1$Quantile, EmpiricalCoverage = resultsPitSTQR0.1$EmpiricalCoverage, Group = "QR 0.1", Color = qr_colors[1]),
   data.frame(Quantile = resultsPitSTQR0.2$Quantile, EmpiricalCoverage = resultsPitSTQR0.2$EmpiricalCoverage, Group = "QR 0.2", Color = qr_colors[2]),
   data.frame(Quantile = resultsPitSTQR0.3$Quantile, EmpiricalCoverage = resultsPitSTQR0.3$EmpiricalCoverage, Group = "QR 0.3", Color = qr_colors[3]),
   data.frame(Quantile = resultsPitSTQR0.4$Quantile, EmpiricalCoverage = resultsPitSTQR0.4$EmpiricalCoverage, Group = "QR 0.4", Color = qr_colors[4]),
   data.frame(Quantile = resultsPitSTCQR0.1$Quantile, EmpiricalCoverage = resultsPitSTCQR0.1$EmpiricalCoverage, Group = "CQR 0.1", Color = cqr_colors[1]),
   data.frame(Quantile = resultsPitSTCQR0.2$Quantile, EmpiricalCoverage = resultsPitSTCQR0.2$EmpiricalCoverage, Group = "CQR 0.2", Color = cqr_colors[2]),
   data.frame(Quantile = resultsPitSTCQR0.3$Quantile, EmpiricalCoverage = resultsPitSTCQR0.3$EmpiricalCoverage, Group = "CQR 0.3", Color = cqr_colors[3]),
   data.frame(Quantile = resultsPitSTCQR0.4$Quantile, EmpiricalCoverage = resultsPitSTCQR0.4$EmpiricalCoverage, Group = "CQR 0.4", Color = cqr_colors[4])
 )
 
 # Combine data into a single data frame
 plot_data <- bind_rows(results_list)
 
 # Create staircase plot
 ggplot() +
   geom_step(data = filter(plot_data, Group %in% c("QR 0.1", "QR 0.2", "QR 0.3", "QR 0.4")), aes(x = Quantile, y = EmpiricalCoverage, color = Group), direction = "vh", size = 1) +
   geom_step(data = filter(plot_data, Group %in% c("CQR 0.4")), aes(x = Quantile, y = EmpiricalCoverage, color = Group), direction = "vh", size = 1) +
   scale_color_manual(values = c("QR 0.1" = qr_colors[1], "QR 0.2" = qr_colors[2], "QR 0.3" = qr_colors[3], "QR 0.4" = qr_colors[4], "CQR 0.4" = cqr_colors[4])) +
   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) +
   labs(title = "n1 = 198, p/n = [0.1,0.2,0.3,0.4]", x = "Quantile", y = "Empirical Coverage", color = "Legend") +
   theme_minimal() +
   theme(
     plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
     axis.title = element_text(size = 9),
     legend.title = element_text(size = 7),
     legend.text = element_text(size = 6)
   )
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 # Add a legend
 # legend("bottomright", legend = c("CQR 0.1", "QR 0.1", "CQR 0.2", "QR 0.2", "CQR 0.3", "QR 0.3", "CQR 0.4", "QR 0.4"), 
 #        col = c(cqr_colors, qr_colors), pch = 16, bty = "n", cex = 1.5, pt.cex = 1.5)
 
 

 
 filename <- paste("plotsfumato198", ".RData",sep="")
 cat(paste("Saving results to file", filename, "\n"))
 
 # Save all the variables to the .RData file
 save(
   resultsPitSTCQR0.1,
   resultsPitSTCQR0.2,
   resultsPitSTCQR0.3, 
   resultsPitSTCQR0.4,
   resultsPitSTQR0.1,
   resultsPitSTQR0.2,
   resultsPitSTQR0.3, 
   resultsPitSTQR0.4,
   file=filename
 )
 
  filename <- paste("plotsfumato198",".RData", sep="")
  cat(paste("Loading results from file", filename, "\n"))
  
  load(filename)
 
 
 
 
 
 
 
 
#------------ PLOT RESULTS WITH n1 < 98
 library(ggplot2)
 

 # Manually creating the data frame
 data <- data.frame(
   Type = rep(c("QR", "CQR"), each = 14),
   Num = rep(c(198, 148, 98, 88, 78, 68, 58, 48, 38, 28, 25, 21, 18, 17), 2),  # Repeat this sequence for each QR and CQR
   E = c(0.029678421, 0.03048579, 0.03671421, 0.032993684, 0.03562316, 0.03934, 0.0417137, 0.04438684, 0.043571579, 0.05242316, 0.05449211, 0.06234, 0.07019316, 0.07412526,  # QR values
         0.01013421, 0.01192368, 0.00635947, 0.0147, 0.00917105, 0.01271421, 0.01647632, 0.02045842, 0.02559211, 0.03197684, 0.04342895, 0.04490579, 0.04297947, 0.04790789)  # CQR values
 )
 data <- data[!(data$Num %in% c(198, 148)),]
 qr_data <- filter(data, Type == "QR")
 cqr_data <- filter(data, Type == "CQR")
 
 # Calculate the difference
 difference_data <- data.frame(
   Type = "Difference",
   Num = qr_data$Num,
   E = qr_data$E - cqr_data$E
 )
 
 # Combine all data
 final_data <- bind_rows(data, difference_data)
 
 # Plotting using ggplot2
 ggplot(final_data, aes(x = Num, y = E, color = Type, group = Type)) +
   geom_point() +
   geom_line() +
   labs(title = "Number of Observations vs. MAE for p/n = 0.1",
        x = "Number of Observations", y = "MAE", color = "Type") +
   theme_minimal() +
   scale_x_reverse(breaks = c(98, 88, 78, 68, 58, 48, 38, 28, 17))  # Adjust breaks to exclude 198 and 148
 