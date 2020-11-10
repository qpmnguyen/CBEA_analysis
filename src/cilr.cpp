#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

NumericVector cilr(NumericMatrix X, NumericMatrix A){
  
}



NumericVector evaluate(NumericVector score, 
                       NumericMatrix X, 
                       NumericMatrix A, float thresh, Int32 nperm){
  int ncolumns = X.ncol();
  IntegerVector sequence = seq_len(ncolumns);
  for (int i=0;i<=nperm;i=i+1){
    IntegerVector samp = sample(sequence, ncolumns, false);
    NumericMatrix X_samp = X(_, sequence);
    
  }
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
