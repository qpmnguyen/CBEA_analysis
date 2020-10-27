# Script to implement all the modelling and evaluation 
# TODO: Add zero imputation using zCompositions
# TODO: A more elegant way of handling preprocessing options
# TODO: Make loop in cilr faster using Rcpp
# TODO: 
library(compositions)
library(fitdistrplus)
library(sn)
library(mixtools)

#' @title Function to perform simple cilr transformation for a set
#' @param X Matrix of n x p dimensions 
#' @param A Matrix of p x m dimensions
simple_cilr <- function(X, A, resample = T, preprocess = T, pcount = NULL, transform = NULL, 
                        abs = F, method = "cdf", ...){
  if(preprocess == T){
    if (missing(pcount)){
      message("Adding default pseudocount of 1")
      pcount <- 1
    }
    message("Pre-processing...")
    message(glue("Adding pseudocount of {pcount}", pcount = pcount))
    if (missing(transform)){
      message("No transformation!")
    } else {
      message(glue("Performing transformation of {trans}", trans = transform))
    }
    X <- process(X, pcount = pcount, transform = transform)
  }
  M <- X[,sample(1:ncol(X))] # shuffling columns since we use index 
  R <- matrix(0, ncol = ncol(A), nrow = nrow(X))
  p <- ncol(X)
  for (i in seq(ncol(A))){
    size <- length(which(A[,i] == 1))
    scale <- sqrt((size*(p - size))/p)
    if(scale == 0){
      warning("scale is 0 here")
    }
    #TODO: Deal with Singleton Sets 
    num <- geometricmeanRow(x = as.matrix(X[,A[,i] == 1]))
    denom <- geometricmeanRow(x = X[,A[,i] == 0])
    cilr <- scale * log(num/denom)
    if (resample == T){
      if (missing(method)){
        stop("Need to require score adjustment if resample is TRUE")
      }
      cilr_resamp <- scale * log(geometricmeanRow(x = as.matrix(M[,A[,i] == 1]))/geometricmeanRow(x = M[,A[,i] == 0]))
      # Fitting a normal mixture distribution
      fit <- normalmixEM(cilr_resamp, ...)
      parm <- list(mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda)
      if (method == "cdf"){
        R[,i] <- pmnorm(cilr, parm = parm)
      } else if (method == "zscore"){
        mean <- sum(fit$lambda * fit$mu)
        sd <- sqrt(sum((fit$sigma + fit$mu - mean)*fit$lambda))
        R[,i] <- (cilr - mean) * (1/sd)
      }
    } else {
      R[,i] <- cilr
    }
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
  p_val <- map_dfc(1:ncol(scores), ~get_p_values(scores = scores[,.x], distr = distr, param = param, alt = alt))
  colnames(p_val) <- colnames(scores)
  if (return == "sig"){
    return(ifelse(p_val < thresh, 1,0))
  } else {
    return(p_val)
  }
}


estimate_distr <- function(data, distr, init){
  dist <- tryCatch({
    if (missing(init)){
      if(distr == "t"){
        init <- list(df = 1)
      } else if (distr == "norm"){
        init <- list(mean = 0, sd = 1)
      } else if (distr == "st"){
        init <- list(xi = 0, omega = 1, alpha = 1, nu = 1)
      }
    }
    message(glue("Fitting permuted null on the {d} distribution", d = distr))
    if (distr != "mnorm"){
      fitdistrplus::fitdist(data, distr = distr, method = "mle", start = init, control = list(maxit = 1000))
    } else {
      fit <- normalmixEM(x = data)
      list(estimate = list(mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda))
    }
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
      } else if (distr == "mnorm"){
        p_val <- 2*(1 - pmnorm(scores, parm = param))
      } else if (distr == "st"){
        p_val <- 2*(1 - pst(scores, xi = param['xi'], omega = param['omega'], alpha = param['alpha'],
                            nu = param['nu']))
      }
    } else {
      message("Using 1-sided test")
      if (distr == "norm"){
        p_val <- 1 - pnorm(scores, mean = param['mean'], sd = param['sd'])
      }else if (distr == "t"){
        p_val <- 1 - pt(scores, df = param['df'])
      } else if (distr == "mnorm"){
        p_val <- 1- pmnorm(scores, parm = param)
      } else if (distr == "st"){
        p_val <- 1 - pst(scores, xi = param['xi'], omega = param['omega'], alpha = param['alpha'],
                            nu = param['nu'])
      }
    }
  }
  return(p_val)
}


#' Mixture p, d, r, q functions for mixture normals 
#' Quantile
qmnorm <- function(p, parm, log=FALSE){
  p <- as.vector(p)
  if(all(names(parm) == c('mu', 'sigma', 'lambda')) == FALSE){
    stop("Parameters requires mu, sigma and lambda")
  }
  n_components <- length(parm$sigma)
  print(paste(n_components, "components!"))
  comp <- vector(mode = "list", length = n_components)
  for (i in 1:n_components){
    comp[[i]] <- parm$lambda[i] * qnorm(p, parm$mu[i], parm$sigma[i], log.p = log)
  }
  return(Reduce("+", comp))
}

# Density
dmnorm <- function(x, parm, log=FALSE){
  x <- as.vector(x)
  if(all(names(parm) == c('mu', 'sigma', 'lambda')) == FALSE){
    stop("Parameters requires mu, sigma and lambda")
  }
  n_components <- length(parm$sigma)
  print(paste(n_components, "components!"))
  comp <- vector(mode = "list", length = n_components)
  for (i in 1:n_components){
    comp[[i]] <- parm$lambda[i] * dnorm(x, parm$mu[i], parm$sigma[i], log = log)
  }
  return(Reduce("+", comp))
}

# Distribution 
pmnorm <- function(q, parm, log = FALSE){
  q <- as.vector(q)
  if(all(names(parm) == c('mu', 'sigma', 'lambda')) == FALSE){
    stop("Parameters requires mu, sigma and lambda")
  }
  n_components <- length(parm$sigma)
  print(paste(n_components, "components!"))
  comp <- vector(mode = "list", length = n_components)
  for (i in 1:n_components){
    comp[[i]] <- parm$lambda[i] * pnorm(q, parm$mu[i], parm$sigma[i], log.p = log)
  }
  return(Reduce("+", comp))
}


# Random generation 
rmnorm <- function(n, parm){
  if(names(parm) != c('mu', 'sigma', 'lambda')){
    stop("Parameters requires mu, sigma and pmix")
  }
  rnormmix(n = n, lambda = parm$lambda, sigma = parm$sigma, mu = parm$mu)
}

