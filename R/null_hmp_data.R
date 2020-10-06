library(HMP16SData)
library(phyloseq)
library(tidyverse)
library(ANCOMBC)
source("R/")
# TODO: Use curatedMetagenomicsData



stool_16S <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Stool" & VISITNO == 1) %>% as_phyloseq()
stool_16S <- subset_samples(stool_16S,!duplicated(RSID)) %>% filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
stool_16S <- prune_samples(sample_sums(stool_16S) >= 1000, stool_16S)







sample_data(stool)$group <- rbinom(nsamples(stool), 1, 0.5)

stool <- prune_samples(sample_sums(stool)>= 1000, stool)
mod <- ANCOMBC::ancombc(stool, formula = "group", lib_cut = 1000, group = "group")
which(mod$res$diff_abn$group == T)
