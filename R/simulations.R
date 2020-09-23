library(tidyverse)
library(glue)
library(phyloseq)
library(MCMCpack)
library(MASS)

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
  attributes(data)$sim_type <- "dm" # dirichlet multinomial distribution
  return(data)
}

#' Simulating according to zero inflated negative binomial distribution 
#' @param n_samp Number of samples
#' @param b_spar Base sparsity
#' @param b_corr Base correlation 
#' @param s_corr Set Correlation
#' @param n_tax Number of taxa
#' @param n_inflate Number of taxa with inflated counts 
#' TODO: More than one correlation structure
zinb_simulation <- function(n_samp, spar, b_cor, s_cor, n_tax = 300, n_inflate = 50){
  cop <- rCopula(n_samp, normalCopula(param ))
}


#' TODO: Return phyloseq object 
#' Simple inflation in the method described by McMurdie and Holmes
#' @param data The data fram e
#' @param eff_size Effect size 
#' @param n_inflate Number of inflated taxa 
inflate_simple <- function(data, eff_size, n_inflate, prop = 0.5){
  n_tax <- ncol(data)
  n_samp <- nrow(data)
  #TODO: Add option to have multiple sets 
  tax_inf <- sample(seq(n_tax), size = n_inflate, replace = F)
  samp_inf <- sample(seq(n_samp), size = round(n_samp*prop), replace = F)
  data[samp_inf, tax_inf] <- data[samp_inf,tax_inf] * eff_size
  # generating tax_table
  tax_tab <- matrix("Set2", ncol = 1, nrow = n_tax) 
  tax_tab[tax_inf,1] <- "Set1"
  colnames(tax_tab) <- "SetLevel"
  rownames(tax_tab) <- colnames(data)
  obj <- list(data = data, tax_tab = tax_tab, idx = list(tax = tax_inf, samp = samp_inf))
  return(obj)
}

