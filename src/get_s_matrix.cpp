# include <Rcpp.h>

using namespace Rcpp;

// Getting the S matrix
// [[Rcpp::export]]

NumericMatrix generate_s_matrix(NumericMatrix A){
  int nr = A.nrow();
  int nc = A.ncol();
  NumericMatrix out(nr, nc);
  for (int i=0;i<nc;i++){
    double size = 0;
    NumericVector column = A(_,i);
    for(int j=0;j<nr;j++){
      if(column[j] > 0){
        size++;
      }
    }
    for(int j=0;j<nr;j++){
      if(column[j] > 0){
        column[j] = 1/size;
      }
      else{ // if column[j] < 0
        column[j] = -1/(nr - size);
      }
    }
    out(_,i) = column;
  }
  return(out);
}

// Getting vector sizes 
// [[Rcpp::export]]
NumericVector get_sizes(NumericMatrix A){
  int nr = A.nrow();
  int nc = A.ncol();
  NumericVector out(nc);
  for (int i=0;i<nc;i++){
    NumericVector column = A(_,i);
    double size = 0;
    for (int j=0;j<nr;j++){
      if (column[j] > 0){
        size++;
      }
    }
    // std::cout << size << std::endl;
    float value = sqrt((size*nr - size*size)/nr);
    out[i] = value;
  }
  return(out);
}
