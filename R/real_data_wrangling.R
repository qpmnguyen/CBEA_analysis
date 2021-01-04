library(tidyverse)
library(HMP16SData)
library(curatedMetagenomicData)
library(phyloseq)

# 16S data stool
data <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Stool" & VISITNO == 1) %>% as_phyloseq()
data <- subset_samples(data,!duplicated(RSID)) %>% 
  filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
data <- prune_samples(sample_sums(data) >= 1000, data)

saveRDS(data, "data/hmp_stool_16S.rds")

# WGS data stool
data <- curatedMetagenomicData(x = "HMP_2012.metaphlan_bugs_list.stool", dryrun = F, bugs.as.phyloseq = T) 
data <- data[[1]]
# subset out duplicated subject ids, keep only those who are healthy, filter for those who are not missing for more than 
# 90%, transform percentage back to proportions 
data <- data %>% subset_samples(!duplicated(subjectID)) %>% subset_samples(disease = "healthy") %>%
  filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE) %>% transform_sample_counts(function(x) x/100)
# Then use number_reads to extrapolate back the absolute counts per bug 
otu_tab <- otu_table(data)
samp_names <- colnames(otu_tab) # keep sample names 
otu_tab <- otu_tab %>% as("matrix") %>% t() %>% as.data.frame() %>% 
  mutate(num_reads = data@sam_data$number_reads) %>% 
  summarise(across(-num_reads, ~ round(.x * num_reads, digits = 0))) %>% t()
colnames(otu_tab) <- samp_names
otu_table(data) <- otu_table(otu_tab, taxa_are_rows = TRUE)
saveRDS(data, "data/hmp_stool_wgs.rds")


# 16S data supragingival and subgingival
annotation <- read.table(file = "../sc2meta/data/genera_methabolism.tsv", sep = "\t", header = TRUE)
data <- V35() %>% subset(select = HMP_BODY_SUBSITE %in% c("Supragingival Plaque", "Subgingival Plaque") & VISITNO == 1) %>% as_phyloseq()
data <- subset_samples(data, !duplicated(RSID)) %>% filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
data <- prune_samples(sample_sums(data) >= 1000, data)

which(data@sam_data@.Data[[6]] == "Supragingival Plaque")
