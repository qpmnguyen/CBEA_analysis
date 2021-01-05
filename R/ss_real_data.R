library(phyloseq)
source("R/cilr.R")
source("R/utils.R")

import <- readRDS(file = "data/hmp_supergingival_supragingival_16S.rds")
data <- import$physeq
annotations <- import$annotation


df <- phylo2cilr(data, "GENUS")

scores <- generate_alt_scores(X = df$X, A = df$A, method = "ssgsea", preprocess = T, transform = "prop", pcount = 1)
scores <- cilr(X = df$X, A = df$A, output = "sig", distr = "norm", perm = 5, adj = T, resample = T)
scores
