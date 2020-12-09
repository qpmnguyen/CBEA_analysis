library(tidyverse)
library(glue)
library(phyloseq)
library(MASS)
library(VGAM)
library(furrr)


#' Shorthand to create parameters
create_parameters <- function(params){
  par <- cross_df(params)
  par <- par %>% mutate(id = seq(1:nrow(par))) %>% group_by(id) %>% nest()
  par <- par %>% transmute(param = data)
  return(par)
}


#' Simulating according to zero inflated negative binomial distribution 
#' @param n_samp Number of samples
#' @param spar Additive sparsity  
#' @param b_rho Baseline inter-taxa correlation 
#' @param s_rho Correlation within the set 
#' @param rho_ratio Ratio of correlation between the set and baseline   
#' @param n_tax Number of taxa  
#' @param n_inflate Number of taxa with inflated counts   
#' @param n_sets Number of sets to simulate   
#' @param prop_set_inflate Number of proportions of sets that are inflated   
#' @param prop_inflate The number of differentially abundant taxa within a set that is inflated 
#' @param samp_prop The proportion of samples have an inflated taxa 
#' @param eff_size Effect size that is a multiplier to the base mu of a sample 
#' @param method Method of simulation. Can be "compensation", "non-compensation" or "normal"
#' @param parameters Path to rds file containing estimated values 
#' @param vary_params Whether we stochastically draw parameter values for each feature  
#' TODO: More than one correlation structure
#' TODO: Adjust rho_ratio calculation when b_rho = 0 (no background correlation)
#' TODO: Add shuffling to make sure that everything is randomized 
#' TODO: Add a way to sample from empirical distribution of zinb values 
zinb_simulation <- function(n_samp, spar, s_rho, eff_size, 
                            b_rho = 0, n_tax = 300, n_inflate = 50, n_sets = 1, 
                            prop_set_inflate = 1, prop_inflate = 1, samp_prop = 0.5, 
                            method = "compensation", vary_params=TRUE, parameters=NULL){
  # generate the the diagnonal matrix
  sigma <- diag(n_tax)
  sigma[sigma == 0] <- b_rho
  print(dim(sigma))
  set_sigma <- sigma[1:n_inflate, 1:n_inflate]
  set_sigma[set_sigma != 1] <- s_rho
  sigma[1:n_inflate,1:n_inflate] <- set_sigma
  
  # First create mvrnorm variables with correlation set by sigma
  margins <- pnorm(mvrnorm(n = n_samp, mu = rep(0, n_tax), Sigma = sigma))
  # Second, set marginals
  # default marginals for negative binomial is from size of 0.595 and mu of 6.646 from HMP data
  true_size <- round(n_inflate * prop_inflate, 0)
  
  # Create parameters based on situation
  if (vary_params == T){
    if (is.null(parameters)){
      message("Randomly sample means from 1 to 10 and sizes from 1 to 5")
      means <- runif(n_tax, 3,5)
      sizes <- runif(n_tax, 1,3)
    } else {
      estimated <- readRDS(file = parameters)
      means <- sample(estimated$mean, size = n_tax, replace = T)
      sizes <- sample(estimated$size, size = n_tax, replace = T)
    }
  } else {
    message("Setting mean to be constant at 3.06 and size at 1.67 estimated from HMP data...")
    means <- rep(3.05,n_tax)
    sizes <- rep(1.67, n_tax)
  }
  
  
  # first n_samp * samp_prop samples will always be inflated 
  inf_size <- round(n_samp * samp_prop,0)
  
  # Number of inf_taxa 
  inf_tax <- round(n_inflate * n_sets * prop_set_inflate, 0)

  # inflating samples  
  suppressMessages(
  inf_samples <- map_dfc(seq(n_tax), .f = function(.x){
    if(method == "compensation"){
      prop <- 1/(eff_size + 1)
      # randomly select first prop of those selected to be upregulated
      idx <- sample(seq(inf_tax), size = round(prop * inf_tax, 0), replace = F)
      remainder <- seq(inf_tax)[-idx]
      a <- sum(means[idx])
      b <- sum(means[remainder])
      # randomly select second prop of those selected to be downregulated 
      if (.x %in% idx){
        result <- qnbinom(p = margins[seq(inf_size),.x], size = sizes[.x], 
                            mu = means[.x]*eff_size)
      } else if (.x %in% remainder){
        result <- qnbinom(p = margins[seq(inf_size),.x], size = sizes[.x],
                            mu = means[.x]*((a/b)*(1-eff_size) + 1))
      } else {
        result <- qnbinom(p = margins[seq(inf_size), .x], size = sizes[.x], 
                            mu = means[.x])
      }
    } else if (method == "no_compensation"){
      # randomly select first prop of those selected to be upregulated
      idx <- sample(seq(inf_tax), size = inf_tax/2, replace = F)
      remainder <- seq(inf_tax)[-idx]
      # randomly select second prop of those selected to be downregulated 
      if (.x %in% idx){
        result <- qnbinom(p = margins[seq(inf_size),.x], size = sizes[.x], 
                            mu = means[.x]*eff_size)
      } else if (.x %in% remainder){
        result <- qnbinom(p = margins[seq(inf_size),.x], size = sizes[.x],
                            mu = means[.x]/eff_size)
      } else {
        result <- qnbinom(p = margins[seq(inf_size), .x], size = sizes[.x], 
                            mu = means[.x])
      }
      
    } else if (method == "normal"){
      if (.x %in% seq(inf_tax)){
        result <- qnbinom(p = margins[seq(inf_size),.x], size = sizes[.x], 
                            mu = means[.x]*eff_size)
      } else {
        result <- qnbinom(p = margins[seq(inf_size), .x], size = sizes[.x], 
                            mu = means[.x])
      }
    }
    return(result)
  })
  )
  # not inflated samples 
  suppressMessages(
    notinf_samples <- map_dfc(seq(n_tax), .f = function(.x){
      result <- qnbinom(p = margins[-seq(inf_size),.x], size = sizes[.x], 
                          mu = means[.x])
      return(result)
    })
  )
  message("Completed loop!")
  if (inf_size == n_samp){
    abundance <- inf_samples
  } else {
    abundance <- rbind(inf_samples, notinf_samples)
  }
  abundance <- as.matrix(abundance)
  if (!is.null(spar)){
    zeroes <- rbinom(length(abundance), size = 1, prob = 1 - spar)
    abundance <- abundance * zeroes
  }
  label <- c(rep(1, inf_size), rep(0, n_samp - inf_size))
  
  colnames(abundance) <- glue("Tax{i}", i = seq(n_tax))
  rownames(abundance) <- glue("Samp{i}", i = seq(n_samp))
  if (n_sets > 1){
    A <- diag(n_sets)
    vec <- as.matrix(rep(1, n_inflate))
    A <- kronecker(A,vec)
    
  } else {
    #TODO if n_sets * n_inflate != n_tax, then there are issues.  
    message("Only one set!")
    A <- rep(0,n_tax)
    A[1:n_inflate] <- 1
    A <- as.matrix(A)
  }
  colnames(A) <- glue("Set{i}ss", i = 1:n_sets)
  rownames(A) <- colnames(abundance)
  sets_inf <- rep(0, n_sets)
  sets_inf[seq(round(n_sets * prop_set_inflate,0))] <- 1
  output <- list(X = as.data.frame(abundance), A = A, label = label, sets_inf = sets_inf)
  return(output)
}

#'  @param modeltype What type of simulations 
#'  @param snr Signal to noise ratio (or effect size)
#'  @param sat Saturation proportion (proportion of sets associated with outcome)
#'  @param ... parameters to pass to the zinb simulation function
sim_prediction <- function(beta_eff, snr, sat, ...){
    # handle defaults 
    def <- list(n_samp = 300, spar = 0.2, s_rho = 0, eff_size = 1, 
                b_rho = 0, n_tax = 2000, n_inflate = 50, n_sets = 40, 
                prop_set_inflate = 0.5, prop_inflate = 1, samp_prop = 0.5, 
                method = "normal", vary_params=TRUE, parameters=NULL)
    if (missing(...)){
        args <- def
    } else {
        sup <- list(...)
        args <- merge_lists(defaults = def, supplied = sup)
    }
    # generate baseline data based on arguments
    baseline <- do.call(zinb_simulation, args)
    # index of samples that are related to the outcome based on model saturation 
    # TODO: Create exceptions for when sat == 0
    set_names <- colnames(baseline$A) 
    if (sat > 0){
        important_sets <- sample(1:length(set_names), size = length(set_names)*sat)
        index_list <- lapply(important_sets, function(x){
            as.vector(which(baseline$A[,x] == 1))
        })
        col_idx <- do.call(c, index_list)
        beta <- rep(rlnorm(1, meanlog = beta_eff, sdlog = 0.5),length(col_idx))
        y <- as.matrix(baseline$X[,col_idx]) %*% as.matrix(beta)
        noise <- rnorm(nrow(baseline$X))
        k <- sqrt(var(y)/(snr * var(noise)))
        y <- as.vector(k) * noise + y
        y <- log(y)
        results_idx <- rep(0, ncol(baseline$X))
        results_idx[col_idx] <- beta[1]
    } else {
        y <- rnorm(nrow(baseline$X))
        results_idx <- rep(0, ncol(baseline$X))
    }
    output <- list(outcome = y, predictors = baseline, idx = results_idx)
}
