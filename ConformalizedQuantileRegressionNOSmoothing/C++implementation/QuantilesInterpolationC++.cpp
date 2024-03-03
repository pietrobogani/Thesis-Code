
// CALL FROM R THE FUNCTION qst
#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/Param.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGS.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSB.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/BFGSMat.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/BKLDLT.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/Cauchy.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/LineSearchBacktracking.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/LineSearchBracketing.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/LineSearchMoreThuente.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/LineSearchNocedalWright.h"
#include "C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/ConformalizedQuantileRegression/C++implementation/LBFGSpp-master/LBFGSpp-master/include/LBFGSpp/SubspaceMin.h"

int findClosestIndex(const Rcpp::NumericVector& vec, double target) {
  Rcpp::NumericVector absDiffs = Rcpp::abs(vec - target);
  return Rcpp::which_min(absDiffs);
}

class QuantilesObjective {
public:
  QuantilesObjective(const Rcpp::NumericVector& QQ, const Rcpp::NumericVector& qqTarg, const Rcpp::IntegerVector& Select, int df)
    : QQ(QQ), qqTarg(qqTarg), Select(Select), current_df(df) {}
  
  double operator()(const Eigen::VectorXd& p, Eigen::VectorXd& grad) {
    Rcpp::Environment sn("package:sn");
    Rcpp::Function qstFunc = sn["qst"];
    
    double sum_sq_error = 0.0;
    for (int i = 0; i < Select.size(); ++i) {
      int selectIndex = Select[i]; 
      
      Rcpp::NumericVector qq = qstFunc(QQ[selectIndex], 
                                       Rcpp::Named("xi", p[0]), 
                                       Rcpp::Named("omega", p[1]), 
                                       Rcpp::Named("alpha", p[2]), 
                                       Rcpp::Named("nu", current_df));
      
      double diff = qqTarg[selectIndex] - qq[0];
      sum_sq_error += diff * diff;
    }
    return sum_sq_error;
  }
  
private:
  Rcpp::NumericVector QQ;
  Rcpp::NumericVector qqTarg;
  Rcpp::IntegerVector Select;
  int current_df;
};

// [[Rcpp::export]]
Rcpp::List QuantilesInterpolationCpp(const Rcpp::NumericVector& qqTarg, const Rcpp::NumericVector& QQ, 
                                     double lc0 = NA_REAL, double sc0 = NA_REAL, double sh0 = NA_REAL) {
  Rcpp::NumericVector LB = Rcpp::NumericVector::create(-20, 0, -30);
  Rcpp::NumericVector UB = Rcpp::NumericVector::create(20, 50, 30);
  
  Rcpp::IntegerVector Select = Rcpp::IntegerVector::create(
    findClosestIndex(QQ, 0.05),
    findClosestIndex(QQ, 0.25),
    findClosestIndex(QQ, 0.50),
    findClosestIndex(QQ, 0.75),
    findClosestIndex(QQ, 0.95)
  );
  
  if (Rcpp::NumericVector::is_na(lc0)) {
    double iqn = R::qnorm(0.75, 0, 1, 1, 0) - R::qnorm(0.25, 0, 1, 1, 0);
    lc0 = qqTarg[findClosestIndex(QQ, 0.50)];
    sc0 = (qqTarg[findClosestIndex(QQ, 0.75)] - qqTarg[findClosestIndex(QQ, 0.25)]) / iqn;
    sh0 = 0;
  }
  
  Rcpp::NumericVector X0 = Rcpp::NumericVector::create(lc0, sc0, sh0);
  Rcpp::NumericMatrix par(30, 3);
  Rcpp::NumericVector ssq(30);
  
  LBFGSpp::LBFGSBParam<double> param;
  param.epsilon = 1e-6;
  param.max_iterations = 100;
  
  for (int df = 1; df <= 30; ++df) {
    LBFGSpp::LBFGSBSolver<double> solver(param);
    QuantilesObjective obj(QQ, qqTarg, Select, df);
    
    Eigen::Map<Eigen::VectorXd> X0_eigen(X0.begin(), X0.size());
    Eigen::Map<Eigen::VectorXd> LB_eigen(LB.begin(), LB.size());
    Eigen::Map<Eigen::VectorXd> UB_eigen(UB.begin(), UB.size());
    
    double fx;
    solver.minimize(obj, X0_eigen, fx, LB_eigen, UB_eigen);
    
    par.row(df - 1) = X0;
    ssq[df - 1] = fx;
  }
  int best_df = Rcpp::which_min(ssq) + 1;
  Rcpp::NumericVector best_params = par.row(best_df - 1);
  
  return Rcpp::List::create(Rcpp::Named("lc") = best_params[0],
                            Rcpp::Named("sc") = best_params[1],
                                                           Rcpp::Named("sh") = best_params[2],
                                                                                          Rcpp::Named("df") = best_df);
}
