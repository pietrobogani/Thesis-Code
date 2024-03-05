#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double objectiveFun(NumericVector params, NumericVector qqTarg, NumericVector QQ, IntegerVector Select, double df) {
  Environment sn = Environment::namespace_env("sn");
  Function qst = sn["qst"];
  double sumSq = 0.0;
  for(int i = 0; i < Select.size(); i++) {
    NumericVector qstResult = qst(Named("p") = QQ[Select[i]], Named("xi") = params[0], Named("omega") = params[1], Named("alpha") = params[2], Named("nu") = df);
    sumSq += pow(qqTarg[Select[i]] - qstResult[0], 2);
  }
  return sumSq;
}

// [[Rcpp::export]]
List findOptimalParams(NumericVector qqTarg, NumericVector QQ, NumericVector initialParams, IntegerVector Select) {
  double bestSsq = R_PosInf;
  NumericVector bestParams = clone(initialParams);
  double bestDf = 0;
  for(int df = 3; df <= 30; df++) {
    double ssq = objectiveFun(initialParams, qqTarg, QQ, Select, df);
    if(ssq < bestSsq) {
      bestSsq = ssq;
      bestDf = df;
    }
  }
  return List::create(Named("lc") = bestParams[0], Named("sc") = bestParams[1], Named("sh") = bestParams[2], Named("df") = bestDf, Named("ssq") = bestSsq);
}
