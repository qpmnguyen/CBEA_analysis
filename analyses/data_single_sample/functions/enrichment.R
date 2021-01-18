library(phyloseq)
library(tidyverse)
library(ggsci)
library(showtext)
source("../../R/utils.R")
source("../../R/cilr.R")

font_add_google("Lora")
showtext_auto()

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
  if (metric == "fdr"){
    significance <- ifelse(results != "Supragingival Plaque",1,0)
  } else {
    significance <- ifelse(results == "Supragingival Plaque", 1, 0)  
  }
  scores <- calculate_statistic(eval = metric, pred = scores[,"Aerobic"], true = significance)
  return(scores)
}

enrichment_analysis <- function(X, A, method, label, metric, ...){
  if (method %in% c("ssgsea", "gsva")){
    scores <- generate_alt_scores(X = X, A = A, method = method, preprocess = T, pcount = 1)
  } else if (method %in% c("wilcoxon")){
    scores <- wc_test(X = X, A = A, thresh = 0.05, preprocess = T, pcount = 1)
  } else {
    scores <- cilr(X = X, A = A, resample = T, ..., maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
  }
  is.matrix(scores)
  output <- enrichment_evaluate(scores = scores, results = label, metric = metric)
  return(output)
}

generate_plots <- function(df, metric){
  df <- df %>% mutate(distr = ifelse(is.na(distr), "None", distr), 
                      adj = ifelse(is.na(adj), FALSE, adj))
  showtext_begin()
  if (metric == "fdr"){
    plt <- ggplot(df, aes(x = models, y = fdr, color = adj, shape = distr)) + 
      geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 1)) + 
      ggsci::scale_color_npg(labels = c("Adjusted", "Unadjusted")) + 
      theme_bw() + theme(text = element_text(family = "Roboto", size = 20)) +
      scale_shape_discrete(labels = c("Mixture Normal", "None", "Normal")) + 
      geom_hline(yintercept = 0.05, col = "red") + 
      labs(x = "Models", y = "Type I error", color = "Correlation adjustment", shape = "Distribution")
  } else if (metric == "pwr"){
    plt <- ggplot(df, aes(x = models, y = pwr, color = adj, shape = distr)) + 
      geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 1)) + 
      ggsci::scale_color_npg(labels = c("Adjusted", "Unadjusted")) + 
      theme_bw() + theme(text = element_text(family = "Roboto", size = 20)) +
      scale_shape_discrete(labels = c("Mixture Normal", "None", "Normal")) + 
      geom_hline(yintercept = 0.05, col = "red") + 
      labs(x = "Models", y = "Power", color = "Correlation adjustment", shape = "Distribution")
  } else if (metric == "auc"){
    df <- df %>% group_by(models, fdr, adj, distr) %>% summarise()
  }
  return(plt)
  showtext_en()
}





