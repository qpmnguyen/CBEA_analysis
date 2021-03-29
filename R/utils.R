# Utility functions  
library(tidyverse)
library(GSVA)
library(pROC)
library(MASS)
library(compositions)
library(mixtools)
library(binom)

#' Calculate test statistics
#' Check if these are vectors
calculate_statistic <- function(eval, pred, true=NULL){
  if (is.vector(pred) == F){
    warning("Coercing pred to vector")
    pred <- as.vector(pred)
  }
  if (is.null(true) == F & is.vector(true) == F){
    warning("Coercing true to vector")
    true <- as.vector(true)
  }
  if(eval %in% c("fdr","pwr")){
    if (missing(true)){
      message("Cannot find true, assuming all values are 1")
      true <- rep(1, length(pred))
    }
    stat <- sum(pred == 1 & true == 1)
    conf <- binom.confint(stat, length(pred), conf.level = 0.95, methods = "ac")
    stat <- conf %>% dplyr::select(c(mean, upper, lower)) %>% as.list() 
    names(stat)[1] <- "estimate"
  } else if (eval == "auc"){
    if(missing(true)){
      stop("Need true values")
    }
    stat <- pROC::ci.auc(response = true, predictor = pred, conf.level = 0.95)
    stat <- list(estimate = stat[2], upper = stat[3], lower = stat[1])
  }
  return(stat)
}


#' @title Generate scores using different models  
generate_alt_scores <- function(X, A, method=c("plage", "zscore", "ssgsea", "gsva", "prop"), 
                                preprocess = T, transform = NULL, pcount = NULL, ...){
  method <- match.arg(method)
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
  if (method %in% c("gsva", "ssgsea", "plage", "z-score")){
    scores <- t(gsva(expr = t(X), gset.idx.list = set, method = method, kcdf = kernel))
  } else if (method == "prop"){
    scores <- aggregate(X = X, A = A)
    colnames(scores) <- colnames(A)
    rownames(scores) <- rownames(X)
  }
  return(as.data.frame(scores))
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
      wilcox.test(x = x[which(A[,i] == 1)], y = x[which(A[,i] != 1)], alternative = alt, ...)$p.value
    })
  }
  R <- ifelse(R <= thresh, 1, 0)
  rownames(R) <- rownames(X)
  colnames(R) <- colnames(A)
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

taxtab_prune <- function(physeq, agg_level){
  # first, grab both otu table and tax tab
  tax <- tax_table(physeq) %>% as("matrix")
  otu <- otu_table(physeq) %>% as("matrix")
  # Generate A matrix as per usual
  A <- taxtab2A(tax, agg_level = agg_level, full = TRUE)
  # Get size of sets  
  set_size <- colsums(A)
  # Separate out index of singletons
  singletons <- which(set_size <= 1)
  
  # Extracting full names from tax table similar to taxtab2A
  id <- which(colnames(tax) == agg_level)
  tax <- tax[,1:id] %>% as.data.frame()
  tax_names <- apply(tax, 1, function(i){
    paste(i, sep = ";_;", collapse = ";_;")
  })
  
  # Match complete full names  
  exclude_ids <- which(tax_names %in% colnames(A)[singletons])
  tax <- tax[-exclude_ids, ]
  
  # Reassign 
  if (taxa_are_rows(physeq)){
    otu <- otu[-exclude_ids,]
  } else {
    otu <- otu[,-exclude_ids]
  }
  
  otu_table(physeq) <- otu_table(otu, taxa_are_rows = taxa_are_rows(physeq))
  tax_table(physeq) <- tax %>% as.matrix()
  return(physeq)
}



#' This function converts a taxonomic table to an A matrix.  
taxtab2A <- function(tax, agg_level, full=TRUE){
  id <- which(colnames(tax) == agg_level)
  tax <- as(tax, "matrix")[,1:id] %>% as.data.frame()
  if (full == TRUE){
    tax_names <- apply(tax,1,function(i){
      paste(i, sep = ";_;", collapse = ';_;')
    })
  } else {
    tax_names <- as.vector(tax[,id])
  }
  
  labels <- unique(tax_names)
  labels <- na.omit(labels)
  labels <- labels[!stringr::str_ends(labels, "NA")] # remove NAs
  A <- matrix(0, ncol = length(labels), nrow = nrow(tax))
  for (i in seq(length(labels))){
    idx <- which(tax_names == labels[i])
    A[idx,i] <- 1
  }
  colnames(A) <- labels
  if (full == T){
    rownames(A) <- names(tax_names)
  } else {
    rownames(A) <- rownames(tax)
  }
  return(A)
}

#' This function aggregates X by simple summation using A matrix 
aggregate <- function(X, A){
  data <- matrix(nrow = nrow(X), ncol = ncol(A))
  for (i in 1:ncol(A)){
    data[,i] <- rowSums(X[,A[,i] == 1])
  }
  colnames(data) <- colnames(A)
  rownames(data) <- rownames(X)
  
  return(data)
}

# Function to convert simulation data to phyloseq type objects to be used in ancom, deseq2 and other related
# packages
sim2phylo <- function(sim){
    tab <- sim$A
    tab <- tab %>% as.data.frame() %>% rownames_to_column(var = "tax") %>% 
            pivot_longer(c(-tax, starts_with("Set")), names_to = "GENUS") %>% 
            filter(value == 1) %>% mutate(DUMMY = paste0("dum", GENUS)) %>% 
            dplyr::select(-value) %>% 
            column_to_rownames(var = "tax") %>%  
            as.matrix() %>% tax_table()
    meta <- sample_data(data.frame(group = sim$label))
    X <- sim$X %>% as.data.frame()
    rownames(X) <- sample_names(meta)
    tax <- otu_table(X, taxa_are_rows = F) 
    physeq <- phyloseq(tax,tab,meta)
  return(physeq)
}


# Function to wrap single sample evaluation into one function  
#' @param n_settings Number of settings, equivalent to number of rows in the parameters data set 
#' @param dir Directory of where the simulation files are generated from 
#' @param method The basic function for the most basic evaluation
#' @param ... Other arguments to functional defined in method 
ss_eval <- function(n_settings, dir, method, ..., ncores=3){
  opt <- furrr_options(seed = TRUE)
  tic()
  plan(multisession, workers = ncores)
  with_progress({
    p <- progressor(steps = n_settings)
    scores <- future_map(1:n_settings, .f = ~{
      p()
      data <- qread(file = glue("{dir}/simulation_{.x}.qs", dir = dir))
      method(X = data$X, A = data$A, ...)
    }, .options = opt)
  })
  plan(sequential)
  toc()
  return(scores)
}


# Fit assuming there is only one set
get_fit <- function(data, adj, distr=c("norm", "mnorm"), init=NULL, ...){
  X <- data$X
  A <- data$A
  print("Estimating distribution")
  # actual scores raw wich is also the distributed 
  scores <- cilr(X = X, A = A, resample = F, preprocess = T, transform = "prop", pcount = 1, adj = F, nperm = 1)
  p <- ncol(X)
  # let's create the permuted and unperumted values 
  perm <- X[,sample(1:p, replace = FALSE)]
  # if adjusted == TRUE generate the unpermuted data with 
  sc_perm <- cilr(X = perm, A = A, resample = F, preprocess = T, transform = "prop", pcount = 1, adj = F, nperm = 1)
  perm_dist <- estimate_distr(as.vector(sc_perm), distr = distr, init = init, ...)
  if (adj == TRUE){ # if correlation adjustment is true 
    print("Adjusting for correlation")
    unperm_dist <- estimate_distr(as.vector(scores), distr = distr, init = init, ...)
    if (distr == "norm"){
      final_distr <- list(mean = perm_dist$mean, sd = unperm_dist$sd)
      dist_name <- "pnorm"
      f <- "dnorm"
      k <- 2
    } else if (distr == "mnorm"){
      final_distr <- get_adj_mnorm(perm = perm_dist, unperm = unperm_dist)
      dist_name <- "pmnorm"
      f <- "dmnorm"
      k <- 5
    }
  } else {
    f <- paste0("d",distr)
    dist_name <- paste0("p", distr)
    final_distr <- perm_dist
    if (distr == "norm"){
      k <- 2
    } else if (distr == "mnorm"){
      k <- 5
    }
  }
  print("Getting output")
  param <- rlist::list.append(final_distr, log = T, x = as.vector(scores))
  loglikelihood <- sum(do.call(f, param))
  ks_param <- rlist::list.append(final_distr, x = as.vector(scores), y = dist_name)
  output <- tibble(ks = do.call(ks.test, ks_param)$statistic, 
                   aic = 2*k - 2*loglikelihood, 
                   bic = k*log(length(scores)) - 2*loglikelihood)
  return(output)
}

#' This function handles the ability to merge supplied and defaults 
merge_lists <- function(defaults, supplied){
    similar_idx <- which(names(defaults) %in% names(supplied))
    if (length(similar_idx) == 0){
      merged <- c(defaults, supplied)
    }
    else {
      merged <- c(defaults[-similar_idx], supplied)
    }
    return(merged)
}

#' Quickly convert phyloseq object to X and A formats 
phylo2cilr <- function(physeq, agg_level){
  X <- otu_table(physeq)
  X <- as(t(X), "matrix")
  A <- taxtab2A(tax = tax_table(physeq), agg_level = agg_level)
  return(list(X = X, A = A))
  
}
