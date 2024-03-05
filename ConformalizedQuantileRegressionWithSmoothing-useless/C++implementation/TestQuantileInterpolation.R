# Load necessary libraries
library(Rcpp)
library(sn)  # Only if you need to use the sn package in R as well

# Source the C++ code
Sys.setenv(PKG_CXXFLAGS="-IC:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include")

Rcpp::sourceCpp("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/QuantilesInterpolationC++.cpp")

# Prepare test data
qqTarg <- seq(-2, 2, length.out = 100)
QQ <- pnorm(qqTarg)  # Example cumulative probabilities

# Call the C++ function
result <- QuantilesInterpolationCpp(qqTarg, QQ)

# Print the result
print(result)
