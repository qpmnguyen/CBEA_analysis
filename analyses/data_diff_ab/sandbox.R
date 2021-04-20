library(targets)
library(tidyverse)
library(tarchetypes)
source("functions/diff_ab_functions.R")

set.seed(12345)

agg_level <- "GENUS"
data <- readRDS(file = "../../data/hmp_supergingival_supragingival_16S.rds")
physeq <- data$physeq
annotation <- data$annotation
physeq <- taxtab_prune(physeq, agg_level = agg_level)

# define group 
print("Define Group")
group <- ifelse(sample_data(physeq)$HMP_BODY_SUBSITE == "Supragingival Plaque", 1, 0)
sample_data(physeq)$group <- group %>% as.factor()

# eval_grid <- tibble(
#     methods = "cilr_wilcox", 
#     distr = "mnorm", 
#     adj = TRUE, 
#     output = "cdf"
# )

A <- taxtab2A(tax_table(physeq), agg_level = "GENUS", full = FALSE)
X <- otu_table(physeq) %>% as("matrix")
X <- t(X)


scores <- cilr(X = X, A = A, resample = T, output = "cdf", distr = "mnorm", adj = T, preprocess = T, transform = "prop", 
               pcount = 1, maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)



