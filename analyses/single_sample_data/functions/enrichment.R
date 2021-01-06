library(phyloseq)
library(here)
library(tidyverse)
here::i_am("analyses/single_sample_data/functions/enrichment.R")
source(here("R", "utils.R"))
source(here("R", "cilr.R"))


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

enrichment_auc <- function(scores, results){
  significance <- ifelse(results == "Supragingival Plaque", 1, 0)
  auc_score <- calculate_statistic(eval = "auc", pred = scores$Aerobic, true = significance)
  return(auc_score)
}

enrichment_analysis <- function(X, A, method, label, ...){
  if (method %in% c("ssgsea", "gsva")){
    scores <- generate_alt_scores(X = X, A = A, method = method, preprocess = T, pcount = 1)
  } else {
    scores <- cilr(X = X, A = A, resample = T, ..., maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
  }
  auc <- enrichment_auc(scores = scores, results = label)
  return(auc)
}



