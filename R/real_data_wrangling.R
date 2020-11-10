library(tidyverse)
library(HMP16SData)
library(curatedMetagenomicData)
library(phyloseq)

# 16S data stool
data <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Stool" & VISITNO == 1) %>% as_phyloseq()
data <- subset_samples(data,!duplicated(RSID)) %>% 
  filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
data <- prune_samples(sample_sums(data) >= 1000, data)

# WGS data stool
data <- curatedMetagenomicData(x = "HMP_2012.metaphlan_bugs_list.stool", dryrun = F, bugs.as.phyloseq = T) 
data <- data[[1]]
data <- data %>% subset_samples(!duplicated(subjectID)) %>% subset_samples(disease = "healthy") %>%
  filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE) 

# 16S data supragingival and subgingival
annotation <- read.table(file = "../sc2meta/data/genera_methabolism.tsv", sep = "\t", header = TRUE)
data <- V35() %>% subset(select = HMP_BODY_SUBSITE %in% c("Supragingival Plaque", "Subgingival Plaque") & VISITNO == 1) %>% as_phyloseq()
data <- subset_samples(data, !duplicated(RSID)) %>% filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
data <- prune_samples(sample_sums(data) >= 1000, data)

which(data@sam_data@.Data[[6]] == "Supragingival Plaque")
