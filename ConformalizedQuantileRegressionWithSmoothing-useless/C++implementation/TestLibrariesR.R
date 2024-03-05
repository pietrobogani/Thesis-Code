library(RcppEigen)
library(Rcpp)

Sys.setenv(PKG_CXXFLAGS="-IC:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include")

Rcpp::sourceCpp("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/testlibraries.cpp")

# Call the function from R
result <- optimizeRosenbrock()
print(result)
# Assuming you've sourced the C++ code already
result <- optimizeQuadratic()

# Print the optimization result
print(result)

# Optimize Himmelblau's function
result_himmelblau <- optimizeFunction("Himmelblau")
print(result_himmelblau)

# Optimize Beale's function
result_beale <- optimizeFunction("Beale")
print(result_beale)




