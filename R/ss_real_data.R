library(phyloseq)
source("R/cilr.R")
source("R/utils.R")

import <- readRDS(file = "data/hmp_supergingival_supragingival_16S.rds")
data <- import$physeq
annotations <- import$annotation

tax_table(data)
