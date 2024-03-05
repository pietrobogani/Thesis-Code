//Dovrebbe essere pronta, devo implementare 3 funzioni sotto e si può usare
  
// #Funzioni aggiuntive da implementare:

// # - qsc DONE
// # - qst_bounds DONE
// # - pst

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>



// Main function to compute quantile
std::vector<double> qst(const std::vector<double>& p, double xi = 0, double omega = 1, double alpha = 0, double nu = std::numeric_limits<double>::infinity(), double tol = 1e-08, std::vector<double> dp = {}, int method = 0) {
  
  
  std::vector<double> q(p.size(), std::numeric_limits<double>::quiet_NaN());
  
  for (size_t i = 0; i < p.size(); ++i) {
   
    
    if (nu == 1) {
      q[i] = qsc(p[i], xi, omega, alpha);
      continue;
    }
    
    if (alpha == std::numeric_limits<double>::infinity()) {
      q[i] = xi + omega * std::sqrt(p[i]); // Placeholder for the correct implementation involving qf
      continue;
    }
    
    if (alpha == -std::numeric_limits<double>::infinity()) {
      q[i] = xi - omega * std::sqrt(1 - p[i]); // Placeholder for the correct implementation involving qf
      continue;
    }
    
    // Handle NA, Inf, -Inf with your logic here, if applicable
    
    double adjustedP = (alpha < 0) ? 1 - p[i] : p[i];
    auto bounds = qst_bounds(adjustedP, std::abs(alpha), nu);
    double lower = bounds[0];
    double upper = bounds[1];
    
    double midpoint, pstValue;
    while ((upper - lower) > tol) {
      midpoint = (lower + upper) / 2;
      pstValue = pst(midpoint, xi, omega, alpha, nu, method) - adjustedP;
      
      if (pstValue > 0) {
        upper = midpoint;
      } else {
        lower = midpoint;
      }
    }
    
    q[i] = xi + omega * ((alpha < 0) ? -1 : 1) * midpoint;
  }
  
  return q;
}




//#ORIGINAL R CODE:


// function (p, xi = 0, omega = 1, alpha = 0, nu = Inf, tol = 1e-08, 
//           dp = NULL, method = 0, ...) 
// {
//   if (!is.null(dp)) {                         I don't need this if
//     if (!missing(alpha)) 
//       stop("You cannot set both component parameters and 'dp'")
//       xi <- dp[1]
//     omega <- dp[2]
//     alpha <- dp[3]
//     nu <- dp[4]
//   }
//   if (length(alpha) > 1) 
//     stop("'alpha' must be a single value")
//     if (length(nu) > 1) 
//       stop("'nu' must be a single value")
//       if (nu <= 0) 
//         stop("'nu' must be non-negative")
//         if (nu > 10000) 
//           return(qsn(p, xi, omega, alpha))
//           if (nu == 1) 
//             return(qsc(p, xi, omega, alpha))
//             if (alpha == Inf) 
//               return(xi + omega * sqrt(qf(p, 1, nu)))
//               if (alpha == -Inf) 
//                 return(xi - omega * sqrt(qf(1 - p, 1, nu)))
//                 na <- is.na(p) | (p < 0) | (p > 1)
//                 abs.alpha <- abs(alpha)
//                 if (alpha < 0) 
//                   p <- (1 - p)
//                   zero <- (p == 0)
//                   one <- (p == 1)
//                   x <- xa <- xb <- xc <- fa <- fb <- fc <- rep(NA, length(p))
//                   nc <- rep(TRUE, length(p))
//                   nc[(na | zero | one)] <- FALSE
//                 fc[!nc] <- 0
//                 bounds <- qst_bounds(p[nc], abs.alpha, nu)
//                   xa[nc] <- bounds[, "lower"]
//                 xb[nc] <- bounds[, "upper"]
//                 fa[nc] <- pst(xa[nc], 0, 1, abs.alpha, nu, method = method, 
//                                                ...) - p[nc]
//                 fb[nc] <- pst(xb[nc], 0, 1, abs.alpha, nu, method = method, 
//                                                ...) - p[nc]
//                 regula.falsi <- FALSE
//                 while (sum(nc) > 0) {
//                   xc[nc] <- if (regula.falsi) 
//                     xb[nc] - fb[nc] * (xb[nc] - xa[nc])/(fb[nc] - fa[nc])
//                     else (xb[nc] + xa[nc])/2
//                     fc[nc] <- pst(xc[nc], 0, 1, abs.alpha, nu, method = method) - 
//                     p[nc]
//                     pos <- (fc[nc] > 0)
//                       xa[nc][!pos] <- xc[nc][!pos]
//                     fa[nc][!pos] <- fc[nc][!pos]
//                     xb[nc][pos] <- xc[nc][pos]
//                     fb[nc][pos] <- fc[nc][pos]
//                     fail <- ((xc[nc] - xa[nc]) * (xc[nc] - xb[nc])) > 0
//                     fail[is.na(fail)] <- TRUE
//                     xc[fail] <- NA
//                     x[nc] <- xc[nc]
//                     nc[fail] <- FALSE
//                     nc[(abs(fc) < tol)] <- FALSE
//                     regula.falsi <- !regula.falsi
//                 }
//                 x <- replace(x, zero, -Inf)
//                   x <- replace(x, one, Inf)
//                   Sign <- function(x) sign(x) + as.numeric(x == 0)
//                   q <- as.numeric(xi + omega * Sign(alpha) * x)
//                   names(q) <- names(p)
//                   return(q)
// }
// <bytecode: 0x00000200b313f250>
//   <environment: namespace:sn>