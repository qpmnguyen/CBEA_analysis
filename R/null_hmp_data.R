library(HMP16SData)
library(phyloseq)
library(tidyverse)
library(ANCOMBC)
library(limma)
library(DESeq2)

# TODO: Use curatedMetagenomicsData
stool_16S <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Stool" & VISITNO == 1) %>% as_phyloseq()
stool_16S <- subset_samples(stool_16S,!duplicated(RSID)) %>% filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
stool_16S <- prune_samples(sample_sums(stool_16S) >= 1000, stool_16S)
for (i in 1:1000){
  
  # first, generate the fake data 
  sample_data(stool_16S)$group <- rbinom(nsamples(stool),1,0.5)
  # Fit models  
  
  
}

#' Function to perform different differential abundance tests aggregated to a certain level 
#' @return A label with names as the taxa level of interest 
physeq_eval <- function(physeq, method, agg_level, params){
  # First, extract out matrices X and A 
  X <- otu_table(physeq) %>% t() %>% as("matrix") %>% as.data.frame()
  A <- taxtab2A(taxtab = tax_tab(physeq), level = agg_level)
  # Second, generate scores for cilr and gsva
  if (stringr::string_detect(method, "cilr")){
    scores <- simple_cilr(X = X, A = A, preprocess = T, pcount = 1)
  } else if (stringr::string_detect(method, "gsva")) {
    gsets <- vector(mode = "list", length = ncol(A))
    for (i in ncol(A)){
      idx <- which(A[,i] == 1)
      gsets[[i]] <- rownames(A)[idx] 
    }
    names(gsets) <- colnames(A)
    scores <- gsva(expr = X, gset.idx.list = gsets, kcdf = "Poisson", method = "gsva")
  # if methods are not gsva or cilr then do normal aggregation 
  } else {
    physeq <- tax_glom(physeq, taxrank = agg_level)
    physeq <- transform_sample_counts(physeq, function(x) ifelse(x == 0, 1, x))
  }
  if (method == "ancombc"){
    mod <- ANCOMBC::ancombc(physeq, formula = "group", lib_cut = 1000, group = "group", p_adj_method = "BH", 
                            alpha = 0.05)
    label <- mod$res$diff_abn[,"group"] * 1 # + 1 to convert from logical to numeric
    names(label) <- as.vector(tax_table(physeq)[rownames(mod$res$diff_abn), agg_level]) # create names for label as actual genus names
  } else if (method == "deseq2"){
    deseq <- phyloseq_to_deseq2(physeq, ~group)
    mod <- DESeq2::DESeq(deseq, test = "Wald", fitType = "parametric")
    res <- DESeq2::result(mod)
    label <- (res$padj < 0.05) * 1
    names(label) <- as.vector(tax_table(physeq)[rownames(res), agg_level])
  } else if (stringr::str_detect(method, "welch")){
    for (i in 1:ncol(A)){
      test <- wilcox.test(x = X[labels == 1,i], X[labels == 0,i])
      result[i] <- test$p.value
    }
  } else if (stringr::str_detect(method, "wilcoxon")){
    
  } else if (stringr::str_detect(method, "voom")){
    
  } 

  return(label)
}


