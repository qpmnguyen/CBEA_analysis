library(HMP16SData)
library(phyloseq)
library(tidyverse)
library(ANCOMBC)
library(limma)
library(DESeq2)
library(corncob)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

# TODO: Use curatedMetagenomicsData
stool_16S <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Stool" & VISITNO == 1) %>% as_phyloseq()
stool_16S <- subset_samples(stool_16S,!duplicated(RSID)) %>% filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
stool_16S <- prune_samples(sample_sums(stool_16S) >= 1000, stool_16S)
sample_data(stool_16S)$group <- rbinom(nsamples(stool_16S),1, 0.5)
stool_16S <- tax_glom(stool_16S, "GENUS")


methods <- c("cilr_wilcox", "cilr_welch", "gsva_wilcox", "gsva_welch", 
             "deseq2", "ancombc")

for (i in 1){
  physeq_eval(stool_16S, method = methods[i], agg_level = "GENUS", param = NULL)
}
 #TODO: The taxtable2A function is a bit too much 

tax <- tax_table(physeq)
id <- which(colnames(tax) == "GENUS")
tax <- as(tax, "matrix")[,1:id]

tax_names <- apply(tax,1, function(i){
  paste(i, sep = ";_;", collapse = ";_;")
})

labels <- unique(tax_names)
A <- matrix(0, ncol = length(labels), nrow = nrow(tax))
for (i in seq(length(labels))){
  idx <- which(tax_names == labels[i])
  A[idx,i] <- 1
}
colnames(A) <- labels
return(A)




for (i in 1:1000){
  
  # first, generate the fake data 
  sample_data(stool_16S)$group <- rbinom(nsamples(stool),1,0.5)
  # Fit models  
  
  
}

#' Function to perform different differential abundance tests aggregated to a certain level 
#' @return A sig with names as the taxa level of interest 
physeq_eval <- function(physeq, method, agg_level, params, thresh=0.05){
  match.arg(method, choices = c("cilr_wilcox", "cilr_welch", "gsva_wilcox", "gsva_welch", 
                                "deseq2", "ancombc"))
  # First, extract out matrices X and A 
  message("Extracting matrices X and A")
  X <- otu_table(physeq) %>% t() %>% as("matrix") %>% as.data.frame()
  A <- taxtab2A(taxtab = tax_table(physeq), level = agg_level)
  
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
    scores <- gsva(t(X), gset.idx.list = gsets, kcdf = "Poisson", method = "gsva", min.sz = 2) %>% t()
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
    sig <- vector(length = ncol(A))
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




