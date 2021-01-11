library(phyloseq)
library(here)
library(tidyverse)
here::i_am("functions/enrichment.R")
source(here("functions", "utils.R"))
source(here("functions", "cilr.R"))


enrichment_processing <- function(dat){
  physeq <- dat$physeq
  annotation <- dat$annotation
  label <- physeq@sam_data$HMP_BODY_SUBSITE
  tab <- tax_table(physeq) %>% as('matrix') %>% as.data.frame()
  otu_names <- rownames(tab)
  colnames(annotation) <- c("GENUS", "METAB")
  tab <- dplyr::left_join(tab, annotation, by = "GENUS")
  tab <- tab %>% as.matrix()
  rownames(tab) <- otu_names
  
  # then, convert taxtable to A matrix  
  A <- taxtab2A(tab, "METAB", full = FALSE)
  
  # after than, retrieve X matrix as pivoted OTU table and re-convert to matrix
  X <- as(otu_table(physeq), "matrix") %>% t()
  
  # Shuffling
  idx <- sample(1:nrow(X), size = nrow(X), replace = F)
  X <- X[idx,]
  label <- label[idx]
  
  return(list(X = X, A = A, label = label))
}

enrichment_evaluate <- function(scores, results, metric){
  significance <- ifelse(results == "Supragingival Plaque", 1, 0)
  scores <- calculate_statistic(eval = metric, pred = scores$Aerobic, true = significance)
  return(auc_score)
}

enrichment_analysis <- function(X, A, method, label, type, ...){
  if (method %in% c("ssgsea", "gsva")){
    scores <- generate_alt_scores(X = X, A = A, method = method, preprocess = T, pcount = 1)
  } else {
    scores <- cilr(X = X, A = A, resample = T, ..., maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
  }
  if (type == "auc"){
    output <- enrichment_evaluate(scores = scores, results = label)
  } else if (type == "pwr"){
    output <- enrichment_evaluate(scores = scores, results = label)
  }
  return(output)
}




