library(tidyverse)
library(glue)
library(phyloseq)
library(MCMCpack)
library(MASS)
library(VGAM)

#' Simple simulation using the Dirichlet-Multinomial distribution
#' TODO: Add functionality to do block sparsity
#' @param template This is the vector of counts that represent the alpha parameters for the Dirichlet distribution
#' @param n_samp Number of samples
#' @param spar Overall sparsity 
#' @param size_vec Give size vector if do not want to use negative binomial distribution 
#' @return A data frame of data 
dm_simulation <- function(template, n_samp, spar, size_vec = NULL){
  n_tax <- length(template)
  message(glue("Number of taxa generated is {n_tax}", n_tax = length(template)))
  if (is.null(size_vec)){
    message("Defaulting to sampling from a negative binomial distribution of mean 5000 and disp. 25")
    size_vec <- rnegbin(n = n_samp, mu = 5000, theta = 25)
  } 
  # TODO: Add check to see if size_vec is actually a vector of length n_samp
  # generate base data
  data <- matrix(nrow = n_samp, ncol = n_tax)
  for (i in seq(n_samp)){
    prob <- MCMCpack::rdirichlet(1, alpha = template)
    data[i,] <- stats::rmultinom(n = 1, size = size_vec[i], prob = prob)
  }
  sparsity <- matrix(rbinom(n_samp * n_tax, size = 1, prob = 1-spar), nrow = n_samp, ncol = n_tax)
  data <- data * sparsity
  data <- as.data.frame(data)     
  colnames(data) <- glue("tax_{num}", num = 1:ncol(data))
  return(data)
}

#' Simulating according to zero inflated negative binomial distribution 
#' @param n_samp Number of samples
#' @param b_spar Base sparsity
#' @param b_rho Base correlation 
#' @param spar_ratio Ratio of sparsity between the set and baseline 
#' @param rho_ratio Ratio of correlation between the set and baseline 
#' @param n_tax Number of taxa
#' @param n_inflate Number of taxa with inflated counts 
#' @param eff_size Effect size that get's added to the base mu of 6.646 from HMP data
#' TODO: More than one correlation structure
zinb_simulation <- function(n_samp, b_spar, b_rho, eff_size, spar_ratio = 1, 
                            rho_ratio = 1, n_tax = 300, n_inflate = 50, prop_inflate = 1,
                            samp_prop = 0.5){
  # generate the the diagnonal matrix
  sigma <- diag(n_tax)
  sigma[sigma == 0] <- b_rho
  set_size <- seq(n_inflate)
  set_sigma <- sigma[set_size, set_size]
  set_sigma[set_sigma != 1] <- s_rho
  sigma[set_size,set_size] <- set_sigma
  # First create mvnorm variables with correlation set by sigma
  margins <- pnorm(mvrnorm(n = n_samp, mu = rep(0, n_tax), Sigma = sigma))
  # Second, set marginals
  # default marginals for negative binomial is from size of 0.595 and mu of 6.646 from HMP data
  true_size <- round(n_inflate * prop_inflate, 0)
  means <- runif(n_tax, 1,10)
  sizes <- runif(n_tax, 0,1)
  
  # first n_inflate taxa will always be inflated 
  inf_size <- round(n_samp * samp_prop,0)
  print(inf_size)
  inf_samples <- map_dfc(seq(n_tax),.f = function(.x){
    if (.x %in% seq(n_inflate)){
      result <- qzinegbin(p = margins[seq(inf_size),.x], size = sizes[.x], munb = means[.x]*eff_size, pstr0 = spar)
    } else {
      result <- qzinegbin(p = margins[seq(inf_size),.x], size = sizes[.x], munb = means[.x], pstr0 = spar)
    }
    return(result)
  })
  
  notinf_samples <- map_dfc(seq(n_tax), .f = function(.x){
    result <- qzinegbin(p = margins[-seq(inf_size),.x], size = sizes[.x], munb = means[.x], pstr0 = spar)
    return(result)
  })
  
  abundance <- rbind(inf_samples, notinf_samples)

  label <- c(rep(1, inf_size), rep(0, n_samp - inf_size))
  
  colnames(abundance) <- glue("Tax{i}", i = seq(n_tax))
  A <- matrix(0, nrow = n_tax, ncol = 1)
  colnames(A) <- "Set1"
  rownames(A) <- colnames(abundance)
  A[set_size,] <- 1
  output <- list(X = abundance, A = A, label = label)
  return(output)
}


#' TODO: Return phyloseq object 
#' Simple inflation in the method described by McMurdie and Holmes
#' @param data The data fram e
#' @param eff_size Effect size 
#' @param n_inflate Number of inflated taxa 
inflate_simple <- function(data, eff_size, n_inflate, prop = 0.5){
  n_tax <- ncol(data)
  n_samp <- nrow(data)
  tax_inf <- sample(seq(n_tax), size = n_inflate, replace = F)
  samp_inf <- sample(seq(n_samp), size = round(n_samp*prop), replace = F)
  # first n_inflate taxa will always be inflated
  data[samp_inf, seq(n_inflate)] <- data[samp_inf,tax_inf] * eff_size
  # generating tax_table
  A <- matrix(0, ncol = 1, nrow = n_tax)
  A[tax_inf] <- 1 
  colnames(A) <- "Set1"
  rownames(A) <- colnames(data)
  
  obj <- list(X = data, A = A, idx = list(tax = tax_inf, samp = samp_inf))
  return(obj)
}

inflate_weiss <- function(){
  
}


