#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

// pst has to be implemented still!!!!!



// [[Rcpp::export]]
NumericMatrix qst_bounds(double p, double alpha, double nu) {
 
  double s = std::signbit(alpha) ? -1 : 1;
  NumericVector lowerVec = Rcpp::qt(NumericVector::create(p), nu, true, false);
  double lower = lowerVec[0]; // Extract the value correctly
  NumericVector upperVec = Rcpp::qf(NumericVector::create(p), 1, nu, true, false);
  double upper = std::sqrt(upperVec[0]); // Extract and calculate correctly
  bool wide = (upper - lower) > 5;
  
  if (wide) {
    double step = 5;
    int m = 0;
    while (true) {
      lower = upper - step;
      double p0 = pst(lower, alpha, nu, 2); // Assuming this call is correct for your pst function.
      if (p0 < p) break;
      step *= std::pow(2, 2.0 / (m + 2));
      ++m;
    }
  }
  
  NumericMatrix result(1, 2);
  if (s > 0) {
    result(0, 0) = lower;
    result(0, 1) = upper;
  } else {
    result(0, 0) = -upper;
    result(0, 1) = -lower;
  }
  
  return result;
}