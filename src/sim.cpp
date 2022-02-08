#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
NumericMatrix sim_data(int ntax, int nsamp,
                       int setsize, int nsets, 
                       float prop_set_inflate, 
                       float prop_samp_inflate, 
                       float eff_size, 
                       float b_rho, 
                       char method, 
                       bool vary_params) {
    return(0);
}

// [[Rcpp::export]]
NumericMatrix get_sigma(int ntax, float b_rho){
    NumericMatrix sigma(ntax);
    std::fill(sigma.begin(), sigma.end(), b_rho);
    return(sigma);
}