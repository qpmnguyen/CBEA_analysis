library(compositions)
library(zCompositions)

#' @title Function to perform simple cilr transformation for a set
#' @param X Matrix of n x p dimensions 
#' @param A Matrix of p x m dimensions
simple_cilr <- function(X, A, process = T){
  #TODO: Add zero imputation using zCompositions 
  if (process == T){
    X = process(X)
  }
  R <- matrix(0, ncol = ncol(A), nrow = nrow(X))
  p <- ncol(X)
  #TODO: Make loop faster
  for (i in seq(ncol(A))){
    size <- length(which(A[,i] == 1))
    scale <- sqrt((size*(p - size))/p)
    if(scale == 0){
      warning("scale is 0 here")
    }
    num <- geometricmeanRow(x = X[,A[,i] == 1])
    denom <- geometricmeanRow(x = X[,A[,i] == 0])
    R[,i] <- scale * log(num/denom)
  }
  colnames(R) <- colnames(A)
  rownames(R) <- rownames(X)
  return(R)
}


process <- function(X){
  X[X == 0] <- 1
  return(unclass(acomp(X)))
}




