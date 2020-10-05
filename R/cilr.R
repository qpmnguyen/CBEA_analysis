# Script to implement all the modelling and evaluation 
# TODO: Add zero imputation using zCompositions
# TODO: A more elegant way of handling preprocessing options
# TODO: Make loop in cilr faster using Rcpp
# TODO: 
library(compositions)
library(zCompositions)
library(fitdistrplus)

#' @title Function to perform simple cilr transformation for a set
#' @param X Matrix of n x p dimensions 
#' @param A Matrix of p x m dimensions
simple_cilr <- function(X, A, abs = FALSE, preprocess = T, pcount = NULL, transform = NULL, method = "raw"){
  if(preprocess == T){
    if (missing(pcount)){
      message("Adding default pseudocount of 1")
      pcount <- 1
    }
    message("Pre-processing...")
    message(glue("Adding pseudocount of {pcount}", pcount = pcount))
    message(glue("Performing transformation of {trans}", trans = transform))
    X <- process(X, pcount = pcount, transform = transform)
  }
  R <- matrix(0, ncol = ncol(A), nrow = nrow(X))
  p <- ncol(X)
  for (i in seq(ncol(A))){
    size <- length(which(A[,i] == 1))
    scale <- sqrt((size*(p - size))/p)
    if(scale == 0){
      warning("scale is 0 here")
    }
    num <- geometricmeanRow(x = X[,A[,i] == 1])
    if (method == "raw"){
      denom <- geometricmeanRow(x = X[,A[,i] == 0])
    } else if (method == "random"){
      rand_idx <- sample(1:ncol(X), size = sum(data$A[,1] == 0), replace = F)
      denom <- geometricmeanRow(x = X[,rand_idx])
    } else if (method == "sample"){
      denom <- geometricmeanRow(x = X)
    }
    
    R[,i] <- scale * log(num/denom)
  }
  colnames(R) <- colnames(A)
  rownames(R) <- rownames(X)
  if (abs == T){
    R <- abs(R)
  }
  return(R)
}

#' Evaluate cilr based on criteria 
cilr_eval <- function(scores, alt="two.sided", distr = "norm", thresh=0.05, resample = T, X=NULL, A=NULL, 
                      return = "sig"){
  if(resample == T){
    if (missing(X)|missing(A)){
      stop("Using the resampling method to generate p-values and scores requires")
    }
    X_perm <- X[,sample(seq(ncol(X)))]
    cilr_perm <- simple_cilr(X = X_perm, A = A, preprocess = T, pcount = 1, transform = NULL)
    cilr_perm <- as.vector(cilr_perm) # convert matrix to one vector 
    param <- estimate_distr(data = cilr_perm, distr = distr)
  } else { # if there is no resampling 
    distr <- "norm"
    param <- c(mean = 0, sd = 1)
  }
  if (return %in% c('sig', 'p-value')){ # returning p-values or significant 0,1 values 
    p_val <- map_dfc(1:ncol(scores), ~get_p_values(scores = scores[,.x], distr = distr, param = param, alt = alt))
    colnames(p_val) <- colnames(scores)
  } else if (return == "sig"){
    p_val <- ifelse(p_val < thresh, 1,0)
  } 
  return(p_val)
}


estimate_distr <- function(data, distr, init){
  dist <- tryCatch({
    if (missing(init)){
      if(distr == "t"){
        init <- list(df = 13)
      } else {
        init <- list(mean = 0, sd = 1)
      }
      message(glue("Fitting permuted null on the {d} distribution", d = distr))
      message(glue("Using default initialization parameters {init}", init = init))
    }
    fitdistrplus::fitdist(data, distr = distr, method = "mle", start = init)
  }, 
  error = function(cond){
    message("There fitting process cannot identify proper distribution parameters")
    message(cond)
    return(NULL)
  })
  if (is.null(dist)){
    param <- NULL
  } else {
    param <- dist$estimate
  }
  return(param)
}

get_p_values <- function(scores, param, alt, distr){
  if (is.null(param)){
    p_val <- rep(NA, length(scores))
  } else {
    if (alt == "two.sided"){
      message("Using 1-sided test")
      scores <- abs(scores)
      if (distr == "norm"){
        p_val <- 2*(1-pnorm(scores, mean = param['mean'], sd = param['sd']))
      } else if (distr == "t"){
        p_val <- 2*(1 - pt(scores, df = param['df']))
      }
    } else {
      message("Using 1-sided test")
      if (distr == "norm"){
        p_val <- 1 - pnorm(scores, mean = param['mean'], sd = param['sd'])
      }else if (distr == "t"){
        p_val <- 1 - pt(scores, df = param['df'])
      }
    }
  }
  return(p_val)
}
