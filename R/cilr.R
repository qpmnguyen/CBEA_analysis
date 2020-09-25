library(compositions)
library(zCompositions)

#' @title Function to perform simple cilr transformation for a set
#' @param X Matrix of n x p dimensions 
#' @param A Matrix of p x m dimensions
simple_cilr <- function(X, A, abs = FALSE, preprocess = T, pcount = NULL, transform = NULL){
  #TODO: Add zero imputation using zCompositions 
  if(preprocess == T){
    if (missing(transform)){
      warning("Performing default transformations into proportions")
      transform <- "prop"
    } 
    if (missing(pcount)){
      warning("Adding default pseudocount of 1")
      pcount <- 1
    }
    message("Pre-processing...")
    message(glue("Adding pseudocount of {pcount}", pcount = pcount))
    message(glue("Performing transformation of {trans}", trans = transform))
    X <- process(X, pcount = pcount, transform = transform)
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
  if (abs == T){
    R <- abs(R)
  }
  return(R)
}

#' @title Some processing functionalities
#' @param X the data frame
#' @param pcount Pseudocount, defaults to 1
#' @param transform What transformation to do. If NULL, leave everything alone. 
#'   Sample include "prop" for proportion, "clr" for centered log-ratio, "log" for simple 
#'   log transformation 
process <- function(X, pcount = 1, transform = "prop"){
  X[X == 0] <- pcount
  if (transform == "prop"){
    X <- unclass(acomp(X))
  } else if (transform == "clr"){
    X <- unclass(clr(X))
  } else if (transform == "log"){
    X <- log(X)
  }
  return(X)
}

#' @title Generate scores using different models  
generate_alt_scores <- function(X, A, method, preprocess = T, transform = NULL, pcount = NULL){
  if (method %in% c("plage", "zscore") & preprocess == F){
    warning("This requires pre-processing into centered log ratio transformation")
    preprocess <- T
    transform <- "clr"
    pcount <- 1
  }
  # doing some preprocessing
  if(preprocess == T){
    if (missing(transform)){
      warning("Performing default transformations into proportions")
      transform <- "prop"
    } 
    if (missing(pcount)){
      warning("Adding default pseudocount of 1")
      pcount <- 1
    }
    message("Pre-processing...")
    message(glue("Adding pseudocount of {pcount}", pcount = pcount))
    message(glue("Performing transformation of {trans}", trans = transform))
    X <- process(X, pcount = pcount, transform = transform)
  }
  # generate list set for GSVA
  set <- vector(mode = "list", length = ncol(A))
  names(set) <- colnames(A)
  for (i in seq(ncol(A))){ # for each set 
    set[[i]] <- rownames(A)[which(A[,i] == 1)]
  }
  # Select correct kernel for GSVA
  if (method == "gsva"){
    if (is.null(transform)){
      kernel <- "Poisson"
    } else if (transform == "clr" | transform == "log"){
      kernel <- "Gaussian"
    } else {
      stop("GSVA doesn't work with proportion transformation")
    }
  } else {
    kernel <- NULL
  }
  scores <- gsva(expr = t(X), gset.idx.list = set, method = method, kcdf = kernel)
  return(t(scores))
}

#' @title Perform wilcoxon rank sum test per sample 
#' @param ... Pass to the wilcox.test argument 
wc_test <- function(X, A, thresh, alt = "two.sided", preprocess = F, transform=NULL, pcount=NULL, ...){
  #TODO: Make sure transform has an option for "none"
  if(preprocess == T){
    if (missing(transform)){
      warning("Performing default transformations into proportions")
      transform <- "prop"
    } 
    if (missing(pcount)){
      warning("Adding default pseudocount of 1")
      pcount <- 1
    }
    message("Pre-processing...")
    message(glue("Adding pseudocount of {pcount}", pcount = pcount))
    message(glue("Performing transformation of {trans}", trans = transform))
    X <- process(X, pcount = pcount, transform = transform)
  }
  R <- matrix(nrow = nrow(X), ncol = ncol(A))
  for (i in seq(ncol(A))){
    R[,i] <- apply(X, 1, function(x){
      wilcox.test(x = x[which(A[,i] == 1)], y = x[which(A[,i] != 1)], alternative = alt, ...)$p.value
    })
  }
}





