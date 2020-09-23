library(compositions)
library(zCompositions)

#' @title Function to perform simple cilr transformation for a set
#' @param X Matrix of n x p dimensions 
#' @param A Matrix of p x m dimensions
simple_cilr <- function(X, A){
  R <- matrix(0, ncol = ncol(A), nrow = nrow(X))
  p <- ncol(X)
  for (i in seq(ncol(A))){
    size <- length(which(A[,i] == 1))
    scale <- sqrt((size*(p - size))/p)
    if(scale == 0){
      warning("scale is 0 here")
    }
    R[,i] <- scale * geometricmeanRow(x = X[,A[,i] == 1])
  }
  colnames(R) <- colnames(A)
  return(R)
}
