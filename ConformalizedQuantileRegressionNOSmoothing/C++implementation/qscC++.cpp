//This should work.


#include <Rcpp.h>



#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <algorithm>

// Helper functions
double sign(double value) {
  if (value > 0) return 1.0;
  if (value < 0) return -1.0;
  return 0.0;
}

std::vector<double> qsc(std::vector<double> p, double xi = 0, double omega = 1, double alpha = 0, std::vector<double> dp = {}) {
  
  for (size_t i = 0; i < p.size(); ++i) {
    
    
    // Main calculation
    double u = (p[i] - 0.5) * M_PI;
    double delta;
    if (std::abs(alpha) == std::numeric_limits<double>::infinity()) {
      delta = sign(alpha);
    } else {
      delta = alpha / std::sqrt(1 + alpha * alpha);
    }
    
    double z = delta / std::cos(u) + std::tan(u);
    
    p[i] = xi + omega * z;
  }
  
  return p;
}




//originale R da tradurre
// function (p, xi = 0, omega = 1, alpha = 0, dp = NULL) 
// {
//   if (!is.null(dp)) {
//     if (!missing(alpha)) 
//       stop("You cannot set both 'dp' and component parameters")
//       xi <- dp[1]
//     omega <- dp[2]
//     alpha <- dp[3]
//   }
//   na <- is.na(p) | (p < 0) | (p > 1)
//     zero <- (p == 0)
//     one <- (p == 1)
//     p <- replace(p, (na | zero | one), 0.5)
//     u <- (p - 0.5) * pi
//   delta <- if (abs(alpha) == Inf) 
//     sign(alpha)
//     else alpha/sqrt(1 + alpha^2)
//       z <- delta/cos(u) + tan(u)
//       z <- replace(z, na, NA)
//       z <- replace(z, zero, -Inf)
//       z <- replace(z, one, Inf)
//       q <- (xi + omega * z)
//       names(q) <- names(p)
//       return(q)
// }