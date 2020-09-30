library(tidyverse)
library(GSVA)
library(ROCR)
library(MASS)
#' Calculate test statistics  
calculate_statistic <- function(eval, pred, true=NULL){
  if(eval == "fdr"){
    stat <- sum(pred == 1)/length(pred)
  } else if (eval == "pwr"){
    if(missing(true)){
      stop("Need true values")
    }
    stat <- sum(pred == 1)/sum(true == 1)
  } else if (eval == "auc"){
    if(missing(true)){
      stop("Need true values")
    }
    stat <- performance(prediction(pred, true),'auc')@y.values[[1]]
  }
  return(stat)
}


#' @title Generate scores using different models  
generate_alt_scores <- function(X, A, method, preprocess = T, transform = NULL, pcount = NULL){
  if (method %in% c("plage", "zscore") & preprocess == F){
    message("This requires pre-processing into centered log ratio transformation")
    preprocess <- T
    transform <- "clr"
    pcount <- 1
  }
  # doing some preprocessing
  if(preprocess == T){
    if (missing(transform)){
      message("Not performing any transformations and leaving it as raw counts")
    } 
    if (missing(pcount)){
      message("Adding default pseudocount of 1")
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
  if(preprocess == T){
    if (missing(transform)){
      message("Not performing any transformations and leaving it as raw counts")
    } 
    if (missing(pcount)){
      message("Adding default pseudocount of 1")
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
      print(x)
      wilcox.test(x = x[which(A[,i] == 1)], y = x[which(A[,i] != 1)], alternative = alt, ...)$p.value
    })
  }
  R <- ifelse(R < thresh, 1, 0)
  return(R)
}

#' @title Some processing functionalities
#' @param X the data frame
#' @param pcount Pseudocount, defaults to 1
#' @param transform What transformation to do. If NULL, leave everything alone. 
#'   Sample include "prop" for proportion, "clr" for centered log-ratio, "log" for simple 
#'   log transformation 
process <- function(X, pcount = 1, transform = NULL){
  X[X == 0] <- pcount
  if (!is.null(transform)){
    if (transform == "prop"){
      X <- unclass(acomp(X))
    } else if (transform == "clr"){
      X <- unclass(clr(X))
    } else if (transform == "log"){
      X <- log(X)
    } 
  } else {
    X <- as.matrix(X)
  }
  return(X)
}

#' Function to get differential abundance 
#' @param sim The simulated object list 
#' @param method includes "wilcoxon", "welch", "ancom", "deseq2", "voom"
get_diff_ab <- function(sim, method){
  A <- sim$A
  X <- sim$X
  label <- sim$label
  if(method == "wilcox"){
    for (i in seq(ncol(A))){
      invisible()
    }
  } else if (method == "ancom"){
    
  } else if (method == "deseq2"){
    
  } else if (method == "voom"){
    
  }
}


