library(tidyr)
library(dplyr)
library(rlang)
library(phyloseq)
library(recipes)
library(Rfast)
library(Rcpp) 
library(tibble)
library(glue)

sourceCpp("src/get_s_matrix.cpp")

#TODO: Add phyloseq object support
#TODO Figure out how to do massive multiplication by scale at once and not having to do loops per set
#TODO: Add parallelization for per-sample calculation here 
#TODO: Add a check for the ordering of the rows of A and the columns of X
#TODO: Convert to package format and add all unit tests 

#' @title Generating the a matrix 
#' @param tab Taxa table, preferably of matrix format - however currently support TaxonomyTable object
#' @param taxlevel Tax level to get the dummy variable of 
#' @param drop_unknown Logical. If \code{drop_unknown == TRUE}, all Unknowns are dropped from the table
generate_a_matrix <- function(tab, taxlevel, drop_unknown = FALSE){
  # processing the tab data frame
  if ("taxonomyTable" %in% class(tab)){
    tab <- tab %>% as("matrix")
  }
  
  level <- as.vector(tab[,taxlevel])
  if (drop_unknown == T){
    level <- na.omit(level)
  }
  unq <- unique(level)
  print(unq)
  A <- matrix(0, nrow = length(unq), ncol = nrow(tab))
  for (i in seq(unq)){
    vec <- ifelse(level == unq[i],1,0)
    A[i,] <- vec
  }
  rownames(A) <- unq
  colnames(A) <- rownames(tab)
  return(t(A))
}

#' @title Calculate the ratio geometric mean
#' @param X is a \code{n} by \code{p} matrix of \code{n} samples and \code{p} variables
#' @param A is a \code{p} by \code{m} matrix of \code{p} variables and \code{m} taxonomic sets
#' @return \code{R}, a \code{n} by \code{m} matrix of \code{n} samples and \code{m} taxonomic sets
generate_r_matrix <- function(X, A){
  # First generate matrix S
  if (is.matrix(X) == F){
    X <- as.matrix(X)
  }
  if (is.matrix(A) == F){
    A <- as.matrix(A)
  }
  A[A == 0] <- -1 # set variables not in the set to be -1
  # loop through each A column and get the adjusted scales
  message("Getting the scales ready...")
  S <- generate_s_matrix(A) # cpp function defined prior 
  message("Calculating the R matrix")
  R <- Rfast::Log(X) %*% S # unscaled variables
  set_size <- get_sizes(A)
  R <- R %*% diag(set_size) # cpp function to get the set size 
  return(R)
}

#' @title Taxonomic aggregation using modifled ILR statistic  
#' @param tax_mat Taxonomic counts of matrix type 
#' @param tax_table Taxonomic table of matrix type 
#' @param tax_level Taxonomic level to aggregate too, must be the same name as column in Taxonomoic Table
#' @param resample Whether or not resampled data is used to generate z-scores
#' @param preprocess Optional. If preprocess = T then the standard preprocessing is done. 
#'   This includes simple conversion to compositional form and adding 1 to all counts to avoid zeros
ilr_agg <- function(tax_mat, tax_table, tax_level, resample = F, verbose = T, preprocess = F){
  if(preprocess == T){
    message("Adding standard data preprocessing...")
    tax_mat <- unclass(acomp(tax_mat + 1)) 
  }
  A <- generate_a_matrix(tax_table, taxlevel = tax_level)
  R <- generate_r_matrix(X = tax_mat, A = A)
  rownames(R) <- rownames(tax_mat)
  colnames(R) <- colnames(A)
  if (resample == T){
    X_perm <- tax_mat[,sample(ncol(tax_mat))] # resampling columns
    suppressMessages(R_perm <- generate_r_matrix(X = X_perm, A = A)) # generating column names 
    message("Fitting distributions each column in R_perm...")
    mean_perm <- matrix(0, ncol = ncol(R), nrow = nrow(R))
    sd_perm <- matrix(0, ncol = ncol(R), nrow = nrow(R))
    for (i in 1:nrow(R_perm)){
      if (verbose == T){ # this will print a lot of messages 
        if (i %% 5 == 0){
          message(glue("Currently at {i} over {m}", i = i, m = nrow(R_perm)))  
        }
      }
      
      dist <- fitdistrplus::fitdist(R_perm[i,], method = "mle", distr = "norm")
      mean_perm[i,] <- rep(dist$estimate[1], ncol(R)) 
      sd_perm[i,] <- rep(1/dist$estimate[2], ncol(R))
    }
    R <- (R - mean_perm) * sd_perm # z-score standardizing
  }
  return(R)
}




