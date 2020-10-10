library(tidyverse)
library(GSVA)
library(ROCR)
library(MASS)
library(limma)


#' Calculate test statistics  
calculate_statistic <- function(eval, pred, true=NULL){
  if(eval %in% c("fdr","pwr")){
    print(length(pred))
    stat <- sum(pred == 1)
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
  if (method %in% c("gsva", "ssgsea", "plage", "z-score")){
    scores <- t(gsva(expr = t(X), gset.idx.list = set, method = method, kcdf = kernel))
  } else if (method == "prop"){
    scores <- aggregate(X = X, A = A)
  }
  return(scores)
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
#' @param X the data set (scores or not scores) aggregated to known sets. This is a data set of n samples and s sets
#' @param labels Sample labels of case and control 
#' @param method includes "wilcoxon", "welch", "ancom", "deseq2", "voom", "corncob", "ancombc"
get_diff_ab <- function(X, A, labels, method, data_type = "simulated"){
  if(method == "wilcox"){
    result <- rep(0, ncol(X))
    for (i in 1:ncol(X)){
      test <- wilcox.test(x = X[labels == 1,i], X[labels == 0,i])
      result[i] <- test$p.value
    }
  } else if (method == "welch"){
    result <- rep(0, ncol(X))
    for (i in 1:ncol(X)){
      test <- t.test(x = X[labels == 1,i], X[labels == 0,i])
      result[i] <- test$p.value
    }
  } else if (method == "ancom"){
    
  } else if (method == "deseq2"){
    
  } else if (method == "voom"){
    
  }
  return(result)
}

#' This function converts a taxonomic table to an A matrix.  
taxtab2A <- function(tax, agg_level){
  id <- which(colnames(tax) == agg_level)
  tax <- as(tax, "matrix")[,1:id] %>% as.data.frame()
  tax_names <- apply(tax,1,function(i){
    paste(i, sep = ";_;", collapse = ';_;')
  })
  labels <- unique(tax_names)
  labels <- labels[!stringr::str_ends(labels, "NA")] # remove NAs
  A <- matrix(0, ncol = length(labels), nrow = nrow(tax))
  for (i in seq(length(labels))){
    idx <- which(tax_names == labels[i])
    A[idx,i] <- 1
  }
  colnames(A) <- labels
  rownames(A) <- names(tax_names)
  return(A)
}

# Function to get data for data type
get_data <- function(type){
  if (type == "16S"){
    library(HMP16SData)
    data <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Stool" & VISITNO == 1) %>% as_phyloseq()
    data <- subset_samples(data,!duplicated(RSID)) %>% 
      filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
    data <- prune_samples(sample_sums(data) >= 1000, data)
  } else if (type == "WGS"){
    # NOT READY
    library(curatedMetagenomicData)
    data <- curatedMetagenomicData(x = "HMP_2012.metaphlan_bugs_list.stool", dryrun = F, bugs.as.phyloseq = T) 
    data <- data[[1]]
    data <- data %>% subset_samples(!duplicated(subjectID)) %>% subset_samples(disease = "healthy") %>%
      filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE) 
  }
  return(data)
}

#' This function aggregates X by simple summation using A matrix 
aggregate <- function(X, A){
  data <- matrix(nrow = nrow(X), ncol = ncol(A))
  for (i in 1:ncol(A)){
    data[,i] <- rowSums(X[,A[,i] == 1])
  }
  return(data)
}

# Function to convert simulation data to phyloseq type objects to be used in ancom, deseq2 and other related
# packages
sim2phylo <- function(sim){
  tab <- sim$A %>% as.data.frame() %>% rownames_to_column(var = "tax") %>% pivot_longer(-tax, "SetLevel") %>% 
    dplyr::select(-value) %>% as.data.frame() %>% column_to_rownames(var = "tax") %>% as.matrix() %>% tax_table()
  meta <- sample_data(data.frame(group = dim$label))
  X <- sim$X %>% as.data.frame()
  rownames(X) <- sample_names(meta)
  tax <- otu_table(X, taxa_are_rows = F) 
  physeq <- phyloseq(tax,tab,meta)
  return(physeq)
}


