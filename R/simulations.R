library(methods)
library(VGAM)

setClass("simulation", 
         slots = c(
           "data" = "data.frame",
           "param" = "list",
           "type" = "vector"
         ))


setGeneric()

new("simulation")

#' @title Function to simulate abundances
#' @description Function to get simulated counts with different sparsity levels 
#'   following a zero-inflated negative binomial distribution  
#' @param n_tax Number of taxa, 
#' @param n_samp sample size
#' @param sparsity Sparsity 
#' @param phi Dispersion Parameter
#' @param mu Mean abundance
null_sim <- function(n_tax, n_samp, sparsity, phi = NULL, mu = NULL){
  if(is.null(phi) | is.null(mu)){
    message("Using HMP16S data to simulate non zero entries")
    phi <- 0.5951434
    mu <- 6.4604577
  } 
  spar <- matrix(rbinom(n_tax * n_samp,prob = 1-sparsity, size = 1), nrow = n_samp, ncol = n_tax)
  non_zero <- matrix(rnbinom(n_tax * n_samp, size = phi, mu = mu), nrow = n_samp, ncol = n_tax)
  simulated_null <- spar * non_zero
  return(simulated_null)
}

simulate_dummy_matrix <- function(n_sets, set_size){
  
}

#' @title Simulating taxonomic data set according to different conditions
#' @description A wrapper function across all individual simulation types
#' @param n_tax Integer. Number of taxa
#' @param n_samp Integer. Number of samples 
#' @param type Character. Type of simulation. \code{null} means no elevated counts for type I error 
#'   calibration. 
#' @param sparsity Numeric. Degree of sparsity \code{=1/p} where \code{p} is the probability of 0 
#' @param n_sets Integer. Number of sets to simulate 
#' @param set_size Integer. Size of the taxonomic set. \code{set_size < n_tax}
#' @param phi Numeric. Dispersion parameter of the negative binomial distribution 
#' @param mu Numeric. Mean parameter of the negative binomial distribution
#' @return \code{data} Phyloseq-type object  
simulation <- function(n_tax, n_samp, type, sparsity, set_size, n_sets = 1, phi = NULL, mu = NULL){
  null_data <- null_sim(n_tax = n_tax, n_samp = n_samp, sparsity = sparsity, phi = phi, mu = mu)
  if (type == "null"){
    if (n_sets != 1){
      warning("If null simulation for type I error, number of sets should not be 0")
    }
    otu_table(null_data, taxa_are_rows = F)
  }
  return(data)
}






