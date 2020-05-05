library(HMP16SData)
library(magrittr)
library(phyloseq)
library(fitdistrplus)
library(tidyverse)

stool_16S <- readRDS(file = "data/stool_16S.rds")

raw_tax <- stool_16S %>% otu_table() %>% t() %>% as("matrix")
tax_vector <- as.vector(raw_tax)
dist_fit <- fitdist(tax_vector[tax_vector != 0], method = "mle", distr = "nbinom")
saveRDS(file = "data/HMP16S_stool_dist_fit.rds", dist_fit)


tax <- stool_16S %>% transform_sample_counts(function(x) x + 1) %>% 
  transform_sample_counts(function(x) x/sum(x)) %>%
  otu_table(taxa_are_rows = T) %>% t() %>% as("matrix")
tax_table <- stool_16S %>% tax_table() %>% as("matrix")