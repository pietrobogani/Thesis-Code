#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]


// Adjust the include paths according to where you placed the LBFGS++ and Eigen libraries
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGS.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSB.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/BFGSMat.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/BKLDLT.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/Cauchy.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/LineSearchBacktracking.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/LineSearchBracketing.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/LineSearchMoreThuente.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/LineSearchNocedalWright.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/Param.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/SubspaceMin.h"


// TEST TO VERIFY IF LIBRARIES ARE CORRECTLY ADDED
{
// Define Himmelblau's function class
class HimmelblausFunction {
private:
  Eigen::VectorXd Vector;
  
public:
  double operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    const double x1 = x(0);
    const double x2 = x(1);
    grad(0) = 4 * x1 * (x1 * x1 + x2 - 11) + 2 * (x1 + x2 * x2 - 7);
    grad(1) = 2 * (x1 * x1 + x2 - 11) + 4 * x2 * (x1 + x2 * x2 - 7);
    return (x1 * x1 + x2 - 11) * (x1 * x1 + x2 - 11) + (x1 + x2 * x2 - 7) * (x1 + x2 * x2 - 7);
  }
};

// Define Beale's function class
class BealesFunction {
private:
  Eigen::VectorXd Vector;
  
public:
  double operator()(const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
    const double x1 = x(0);
    const double x2 = x(1);
    grad(0) = 2 * (1.5 - x1 + x1 * x2) * (-1 + x2) +
      2 * (2.25 - x1 + x1 * x2 * x2) * (-1 + x2 * x2) +
      2 * (2.625 - x1 + x1 * x2 * x2 * x2) * (-1 + x2 * x2 * x2);
    grad(1) = 2 * (1.5 - x1 + x1 * x2) * x1 +
      4 * (2.25 - x1 + x1 * x2 * x2) * x1 * x2 +
      6 * (2.625 - x1 + x1 * x2 * x2 * x2) * x1 * x2 * x2;
    return (1.5 - x1 + x1 * x2) * (1.5 - x1 + x1 * x2) +
      (2.25 - x1 + x1 * x2 * x2) * (2.25 - x1 + x1 * x2 * x2) +
      (2.625 - x1 + x1 * x2 * x2 * x2) * (2.625 - x1 + x1 * x2 * x2 * x2);
  }
};

// [[Rcpp::export]]
Rcpp::List optimizeFunction(std::string functionName) {
  Eigen::VectorXd x(2);
  x << -3, -4; // Initial guess
  
  LBFGSpp::LBFGSBParam<double> param;
  param.max_iterations = 200; // Adjust as needed
  param.epsilon = 1e-6; // Convergence tolerance
  
  LBFGSpp::LBFGSBSolver<double> solver(param);
  
  Eigen::VectorXd lb(2);
  Eigen::VectorXd ub(2);
  lb << -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(); // No lower bounds
  ub << std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity();  // No upper bounds
  
  double fx;
  int niter;
  if (functionName == "Himmelblau") {
    HimmelblausFunction f;
    niter = solver.minimize(f, x, fx, lb, ub);
  } else if (functionName == "Beale") {
    BealesFunction f;
    niter = solver.minimize(f, x, fx, lb, ub);
  } else {
    return Rcpp::List::create(Rcpp::Named("error") = "Unknown function name");
  }
  
  return Rcpp::List::create(Rcpp::Named("x") = Rcpp::wrap(x),
                            Rcpp::Named("f(x)") = fx,
                            Rcpp::Named("iterations") = niter);
}

}