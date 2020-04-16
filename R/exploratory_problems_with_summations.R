library(phyloseq)
library(tidyverse)
library(compositions)

data("GlobalPatterns") 
proc_dat <- GlobalPatterns %>% filter_taxa(function(x) sum(x > 0) >= 0.1 * length(x), TRUE) %>% 
  filter_taxa(function(x) mean(x) > 0.005/100, TRUE)  %>% transform_sample_counts(function(x) x + 1)
proc_tax <- otu_table(proc_dat, taxa_are_rows = T)
proc_tab <- tax_table(proc_dat)
proc_samp <- sample_data(proc_dat)

proc_tab <- as.data.frame(proc_tab) %>% rownames_to_column(var = "taxa_index") %>% as_tibble()
unq_fam <- as.vector(unique(pull(proc_tab, "Family")))
idx <- proc_tab %>% filter(Family == unq_fam[2]) %>% pull(taxa_index)

proc_tax <- unclass(acomp(proc_tax))
proc_tax <- as.data.frame(proc_tax) %>% rownames_to_column(var = "taxa_index") %>% as_tibble()
proc_tax %>% filter(taxa_index %in% idx)[,1 ] 
