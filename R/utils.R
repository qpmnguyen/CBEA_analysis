library(tidyr)
library(dplyr)
library(rlang)
library(phyloseq)
library(recipes)
library(Rfast)
library(Rcpp) 
library(tibble)

sourceCpp("src/get_s_matrix.cpp")



#' @title Generating the a matrix 
#TODO: Add phyloseq object support
#' @param tab Taxa table, preferably of matrix format - however currently support TaxonomyTable object
#' @param taxlevel Tax level to get the dummy variable of 
#' @param drop_unknown Logical. If \code{drop_unknown == TRUE}, all Unknowns are dropped from the table
generate_a_matrix <- function(tab, taxlevel, drop_unknown = FALSE){
  # processing the tab data frame
  if (class(tab) == "taxonomyTable"){
    tab <- tab %>% as("matrix")
  }
  tab <- tab %>% as.data.frame() %>% rownames_to_column(var = "tax_id") %>%
    dplyr::select(c(tax_id,!!quo(taxlevel))) %>%
    mutate_at(vars(-tax_id), function(.x){replace_na(as.character(.x), "Unknown")})
  # make dummy variables
  dummy <- tab %>% recipe(tax_id ~ .) %>% step_dummy(-tax_id) %>% prep(training = tab) %>%
    bake(new_data = tab)
  if (drop_unknown == T){
    dummy <- dummy %>% dplyr::select(-ends_with("Unknown"))
  }
  dummy <- dummy %>% as.data.frame() %>% column_to_rownames(var = "tax_id") %>% data.matrix()
  return(dummy)
}

#' @title Calculate the ratio geometric mean
#' @param X is a \code{n} by \code{p} matrix of \code{n} samples and \code{p} variables
#' @param A is a \code{p} by \code{m} matrix of \code{p} variables and \code{m} taxonomic sets
#' @return \code{R}, a \code{n} by \code{m} matrix of \code{n} samples and \code{m} taxonomic sets
generate_r_matrix <- function(X, A){
  # First generate matrix S
  #TODO Figure out how to do massive multiplication by scale at once and not having to do loops per set
  if (is.matrix(X) == F){
    X <- as.matrix(X)
  }
  if (is.matrix(A) == F){
    A <- as.matrix(A)
  }
  A[A == 0] <- -1 # set variables not in the set to be -1
  # loop through each A column and get the adjusted scales
  message("Getting the scales ready...")
  S <- get_s_matrix(A) # cpp function defined prior 
  message("Calculating the R matrix")
  R <- Rfast::Log(X) %*% S # unscaled variables
  set_size <- get_sizes(A)
  R <- R %*% diag(set_size) # cpp function to get the set size 
  return(R)
}

#' @title Taxonomic aggregation using modifled ILR statistic  
#' @param tax_mat A matrix of taxonomic counts 
#' @param tax_table Taxonomic table
#' TODO: Add support for phyloseq type objects  
ilr_agg <- function(tax_mat, tax_table, tax_level,...,resample = F, 
                    n_resamples = 100, preprocess = F){
  if(preprocess = T){
    message("Adding standard data preprocessing...")
    tax_mat <- unclass(acomp(tax_mat + 1)) 
  }
  A <- generate_a_matrix(tax_table, taxlevel = tax_level,...)
  if (resample = T){
    A_perm <- A
    colnames(A_perm) <- sample(colnames(A), size = 3, replace = F)
  } else {
    R <- generate_r_matrix(X = tax_mat, A = A)
  }
}




