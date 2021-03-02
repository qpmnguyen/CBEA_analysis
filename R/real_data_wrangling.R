library(tidyverse)
library(HMP16SData)
library(curatedMetagenomicData)
library(phyloseq)
library(biomformat)
library(stringr)

#' @title Function to process 16S data 
process_16S <- function(physeq){
  physeq <- physeq %>% filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
  physeq <- prune_samples(sample_sums(physeq) >= 1000, physeq)
  return(physeq)
}

#' @title Function to process WGS data 
process_WGS <- function(physeq){
  physeq <- physeq %>% filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE) %>% 
    transform_sample_counts(function(x) x/100)
  # Then use number_reads to extrapolate back the absolute counts per bug 
  otu_tab <- otu_table(physeq)
  samp_names <- colnames(otu_tab) # keep sample names 
  otu_tab <- otu_tab %>% as("matrix") %>% t() %>% as.data.frame() %>% 
    mutate(num_reads = physeq@sam_data$number_reads) %>% 
    summarise(across(-num_reads, ~ round(.x * num_reads, digits = 0))) %>% t()
  colnames(otu_tab) <- samp_names
  otu_table(physeq) <- otu_table(otu_tab, taxa_are_rows = TRUE)
  physeq <- prune_samples(sample_sums(physeq) >= 1000, physeq)
  return(physeq)
}



# 16S data stool
data <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Stool" & VISITNO == 1) %>% as_phyloseq()
data <- subset_samples(data,!duplicated(RSID)) %>% process_16S()
saveRDS(data, "data/hmp_stool_16S.rds")

# 16S data stool and tongue 
stool <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Stool" & VISITNO == 1) %>% as_phyloseq()
tongue <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Tongue Dorsum" & VISITNO == 1) %>% as_phyloseq()
data <- merge_phyloseq(stool, tongue)
data <- data %>% #subset_samples(!duplicated(RSID)) %>% 
  process_16S()
saveRDS(data, "data/hmp_stool_tongue_16S.rds")


# 16S data supragingival and subgingival
annotation <- read.table(file = "../sc2meta/data/genera_methabolism.tsv", sep = "\t", header = TRUE)

supra <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Supragingival Plaque" & VISITNO == 1) %>% 
  as_phyloseq() 

sub <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Subgingival Plaque" & VISITNO == 1) %>% 
  as_phyloseq() 

merged <- merge_phyloseq(supra, sub)
merged <- merged %>% subset_samples(!duplicated(RSID)) %>% process_16S()

output <- list(
  physeq = merged, 
  annotation = annotation
)
saveRDS(output, "data/hmp_supergingival_supragingival_16S.rds")


# WGS data stool
data <- curatedMetagenomicData(x = "HMP_2012.metaphlan_bugs_list.stool", dryrun = F, bugs.as.phyloseq = T) 
data <- data[[1]]
data <- data %>% subset_samples(!duplicated(subjectID)) %>% subset_samples(disease = "healthy") %>% process_WGS()
saveRDS(data, "data/hmp_stool_wgs.rds")

# wgs tongue dorsum  
oral <- curatedMetagenomicData(x = "HMP_2012.metaphlan_bugs_list.oralcavity", dryrun = F, bugs.as.phyloseq = T) 
stool <- curatedMetagenomicData(x = "HMP_2012.metaphlan_bugs_list.stool", dryrun = F, bugs.as.phyloseq = T) 
oral <- oral[[1]]
stool <- stool[[1]]
oral <- oral %>% subset_samples(body_subsite == "tongue_dorsum")
data <- merge_phyloseq(oral, stool)
data <- data %>% #subset_samples(!duplicated(subjectID)) %>% 
  subset_samples(disease = "healthy") %>% process_WGS()
saveRDS(data, "data/hmp_stool_tongue_wgs.rds")



# Prediction data sets  
# Easy prediction using HMP 16S tongue dorsum and stool  
# Disease prediction using WGS IBD
data <- curatedMetagenomicData(x = "NielsenHB_2014.metaphlan_bugs_list.stool", dryrun = F, bugs.as.phyloseq = T)
data <- data[[1]]
data <- data %>% subset_samples(!duplicated(subjectID)) %>% process_WGS()
saveRDS(data, "data/nielsen_ibd_wgs.rds")


# qiita data sets  
biom <- import_biom("../../datasets/study_1939_020121-133947/processed_data/441_otu_table.biom")
samp_dat <- read.table(file = "../../datasets/study_1939_020121-133947/1939_20180418-110402.txt", sep = "\t", 
                       header = T)
rownames(samp_dat) <- samp_dat$sample_name
sample_data(biom) <- samp_dat

biom <- biom %>% process_16S()
colnames(tax_table(biom)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
map_df(tax_table(biom), ~str_remove(.x, "[a-z]__"))
otu_id <- rownames(tax_table(biom))
taxtab <- apply(tax_table(biom),2,str_remove, "[a-z]__")
taxtab <- ifelse(taxtab == "", NA, taxtab)
tax_table(biom) <- taxtab

saveRDS(biom, "data/ackerman_ibd_16S.rds")
