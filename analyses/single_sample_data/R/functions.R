library(phyloseq)
library(here)
here::i_am("analyses/single_sample_data/R/functions.R")
source(here("R", "utils.R"))
source(here("R", "cilr.R"))

# Function to take the compound object, and generate enrichment sets
process_dat <- function(data){
  physeq <- data$physeq
  annotation <- data$annotation
  
  # first, get annotation as a column in A  
  tab <- tax_table(physeq) %>% as('matrix') %>% as.data.frame()
  otu_names <- rownames(tab)
  colnames(annotation) <- c("GENUS", "METAB")
  tab <- left_join(tab, annotation, by = "GENUS")
  tab <- tab %>% as.matrix()
  rownames(tab) <- otu_names

  # then, convert taxtable to A matrix  
  A <- taxtab2A(tab, "METAB", full = FALSE)
  
  # after than, retrieve X matrix as pivoted OTU table and re-convert to matrix
  X <- as(otu_table(physeq), "matrix") %>% t()
  
  # export label for samples  
  label <- physeq@sam_data$HMP_BODY_SUBSITE
  
  # shuffle  
  idx <- sample(seq_len(length(label)), length(label), replace = FALSE)
  X <- X[idx,]
  label <- label[idx]
  return(list(X = X, A = A , label = label))
}

# function to evaluate models 
evalutate <- function(sig, label){
  supra <- ifelse(label == "Supragingival Plaque", 1, 0)
  sub <- ifelse(label == "Subgingival Plaque", 1, 0)
  
  
  
  
}

score <- cilr(X, A, resample = T, distr = "norm", output = "cdf")
score2 <- generate_alt_scores(X, A, method = "gsva", preprocess = T, pcount = 1)

aero_idx <- which(score$Aerobic == 1)
true <- ifelse(label == "Subgingival Plaque", 1, 0)
true
supra <- which(label == "Supragingival Plaque")

calculate_statistic(eval = "auc", pred = score2$Anaerobic, true = true)


