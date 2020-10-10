library(tidyverse)
library(ANCOMBC)
library(DESeq2)
library(corncob)
library(phyloseq)
library(GSVA)
library(stringr)
source("R/utils.R")


#' Function to perform different differential abundance tests aggregated to a certain level 
#' @param physeq Phyloseq object 
#' @param method Method to perform evaluation. can be of "cilr_wilcox", "cilr_welch", "gsva_wilcox", "gsva_welch", 
#'   "deseq2", "ancombc", "corncob"
#' @param agg_level Level to aggregate to 
#' @param params Additional parameters (NOT USED)
#' @param thresh threshold for significance after BH correction 
#' @return A sig with names as the taxa level of interest 
physeq_eval <- function(physeq, method, agg_level, params=NULL, thresh=0.05){
  match.arg(method, choices = c("cilr_wilcox", "cilr_welch", "gsva_wilcox", "gsva_welch", 
                                "deseq2", "ancombc", "corncob"))
  # First, extract out matrices X and A 
  if (sum(stringr::str_detect(method, c("cilr","gsva"))) >= 1){
    message("Extracting matrices X and A")
    X <- otu_table(physeq) %>% t() %>% as("matrix") %>% as.data.frame()
    A <- taxtab2A(tax = tax_table(physeq), agg_level = agg_level)
  }
  # Second, generate scores for cilr and gsva if cilr or gsva is detected in method
  if (stringr::str_detect(method, "cilr")){
    message("Fitting cILR scores...")
    scores <- simple_cilr(X = X, A = A, preprocess = T, pcount = 1)
  } else if (stringr::str_detect(method, "gsva")) {
    message("Fitting GSVA scores...")
    gsets <- vector(mode = "list", length = ncol(A))
    for (i in 1:ncol(A)){
      idx <- which(A[,i] == 1)
      gsets[[i]] <- rownames(A)[idx] 
    }
    names(gsets) <- colnames(A)
    # transpose because genes are rows
    scores <- gsva(t(X), gset.idx.list = gsets, kcdf = "Poisson", method = "gsva") %>% t()
  } else {
    # if methods are not gsva or cilr then do normal aggregation 
    message("Performing normal physeq aggregation to set level")
    physeq <- tax_glom(physeq, taxrank = agg_level)
    physeq <- transform_sample_counts(physeq, function(x) ifelse(x == 0, 1, x)) # add pseudocount
  }
  
  # Then run models according to specifications in argument method 
  if (method == "ancombc"){
    message("Running ANCOMBC...")
    # running ancombc model 
    mod <- ANCOMBC::ancombc(physeq, formula = "group", lib_cut = 1000, group = "group", 
                            p_adj_method = "BH", 
                            alpha = thresh) 
    sig <- mod$res$diff_abn[,"group"] * 1 # + 1 to convert from logical to numeric
    names(sig) <- as.vector(tax_table(physeq)[rownames(mod$res$diff_abn), agg_level]) # create names for sig as actual genus names
  } else if (method == "deseq2"){
    message("Running DESeq2...")
    # running convert physeq to desq
    deseq <- phyloseq_to_deseq2(physeq, ~factor(group))
    # TODO: Support multiple deseq2 options
    mod <- DESeq2::DESeq(deseq, test = "LRT", fitType = "parametric", 
                         reduced = ~ 1) # running dseq2
    res <- DESeq2::results(mod) # getting results
    sig <- (res$padj < thresh) * 1
    names(sig) <- as.vector(tax_table(physeq)[rownames(res), agg_level])
  } else if (method == "corncob"){
    # Running corncob 
    # TODO: Support multiple corncob options
    mod <- differentialTest(formula = ~ group, phi.formula = ~ group, formula_null = ~ 1, phi.formula_null = ~ group,
                            test = "LRT", boot = FALSE, data = physeq, fdr_cutoff = thresh)
    sig <- (da_analysis$p_fdr < thresh)*1
    names(sig) <- as.vector(tax_table(physeq)[names(sig), agg_level])
  } else if (stringr::str_detect(method, "wilcox")){
    # Wilcoxon rank-sum test 
    sig <- vector(length = ncol(A))
    labels <- sample_data(physeq)$group
    for (i in 1:ncol(A)){
      test <- wilcox.test(x = scores[labels == 1,i], scores[labels == 0,i], alt = "greater")
      sig[i] <- test$p.value
    }
    sig <- p.adjust(sig, method = "BH")
    sig <- ifelse(sig < thresh, 1, 0)
    names(sig) <- colnames(A)
  } else if (stringr::str_detect(method, "welch")){
    # Welch's t-test
    sig <- vector(length = ncol(A)) # preallocate vector 
    labels <- sample_data(physeq)$group
    for (i in 1:ncol(A)){
      test <- t.test(x = scores[labels == 1,i], scores[labels == 0,i], alternative = "greater")
      sig[i] <- test$p.value
    }
    sig <- p.adjust(sig, method = "BH")
    sig <- ifelse(sig < thresh, 1, 0)
    names(sig) <- colnames(A)
  } 
  return(sig)
}