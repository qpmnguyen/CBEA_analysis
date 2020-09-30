library(tidyverse)
library(glue)
library(phyloseq)
library(MCMCpack)
library(MASS)
library(VGAM)
library(furrr)

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
#' @param n_sets Number of sets to simulate 
#' @param prop_set_inflate Number of proportions of sets that are inflated 
#'
#' @param eff_size Effect size that is a multiplier to the base mu of a sample 
#' TODO: More than one correlation structure
#' TODO: Adjust rho_ratio calculation when b_rho = 0 (no background correlation)
#' TODO: Add shuffling to make sure that everything is randomized 
zinb_simulation <- function(n_samp, b_spar, b_rho, eff_size, spar_ratio = 1,
                            rho_ratio = 1, n_tax = 300, n_inflate = 50, n_sets = 1, prop_set_inflate = 1, 
                            prop_inflate = 1,
                            samp_prop = 0.5, parallel = T){
  # generate the the diagnonal matrix
  sigma <- diag(n_tax)
  sigma[sigma == 0] <- b_rho
  set_size <- seq(n_inflate)
  set_sigma <- sigma[set_size, set_size]
  set_sigma[set_sigma != 1] <- b_rho * rho_ratio
  sigma[set_size,set_size] <- set_sigma
  
  # First create mvnorm variables with correlation set by sigma
  margins <- pnorm(mvrnorm(n = n_samp, mu = rep(0, n_tax), Sigma = sigma))
  # Second, set marginals
  # default marginals for negative binomial is from size of 0.595 and mu of 6.646 from HMP data
  true_size <- round(n_inflate * prop_inflate, 0)
  means <- runif(n_tax, 1,10)
  sizes <- runif(n_tax, 0,1)
  
  # first n_samp * samp_prop samples will always be inflated 
  inf_size <- round(n_samp * samp_prop,0)
  
  # Number of inf_taxa 
  inf_tax <- round(n_inflate * n_sets * prop_set_inflate, 0)

  
  # create elevated sample sizes 
  if (parallel == T){
    message("Starting parallel procedure...")
    plan(multiprocess)
  } else {
    plan(sequential)
  }
  suppressMessages(
    inf_samples <- future_map_dfc(seq(n_tax),.f = function(.x){
      if (.x %in% seq(inf_tax)){
        result <- qzinegbin(p = margins[seq(inf_size),.x], size = sizes[.x], 
                              munb = means[.x]*eff_size, pstr0 = b_spar * spar_ratio)
      } else {
        result <- qzinegbin(p = margins[seq(inf_size),.x], size = sizes[.x], 
                            munb = means[.x], pstr0 = b_spar)
      }
      return(result)
    })
  )
  suppressMessages(
    notinf_samples <- future_map_dfc(seq(n_tax), .f = function(.x){
      result <- qzinegbin(p = margins[-seq(inf_size),.x], size = sizes[.x], 
                          munb = means[.x], pstr0 = b_spar)
      return(result)
    })
  )
  plan(sequential)
  message("Completed loop!")
  if (inf_size == n_samp){
    abundance <- inf_samples
  } else {
    abundance <- rbind(inf_samples, notinf_samples)
  }
  label <- c(rep(1, inf_size), rep(0, n_samp - inf_size))
  
  colnames(abundance) <- glue("Tax{i}", i = seq(n_tax))
  A <- diag(n_sets)
  vec <- as.matrix(rep(1, n_inflate))
  A <- kronecker(A,vec)
  print(dim(A))
  colnames(A) <- glue("Set{i}", i = 1:n_sets)
  rownames(A) <- colnames(abundance)
  sets_inf <- rep(0, n_sets)
  sets_inf[seq(round(n_sets * prop_set_inflate,0))] <- 1
  output <- list(X = abundance, A = A, label = label, sets_inf = sets_inf)
  return(output)
}


#' Shorthand to create parameters
create_parameters <- function(params){
  par <- cross_df(params)
  par <- par %>% mutate(id = seq(1:nrow(par))) %>% group_by(id) %>% nest()
  par <- par %>% transmute(param = data)
  return(par)
}

