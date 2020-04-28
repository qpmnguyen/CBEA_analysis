library(phyloseq)
library(compositions)
library(tidyverse)
library(fitdistrplus)
library(randomForest)
library(patchwork)
source("R/utils.R")

data(GlobalPatterns)

dummy <- generate_a_matrix(tab = tax_table(GlobalPatterns), taxlevel = "Genus")
tax <- GlobalPatterns %>% transform_sample_counts(function(x) x + 1) %>% 
  transform_sample_counts(function(x) x/sum(x)) %>%
  otu_table(taxa_are_rows = T) %>% t() %>% as("matrix")

glom <- GlobalPatterns %>% transform_sample_counts(function(x) x+1) %>% 
  tax_glom(taxrank = "Genus") %>%  
  transform_sample_counts(function(x) x/sum(x)) %>% 
  otu_table(taxa_are_rows = T) %>% t() %>% as("matrix")

output <- GlobalPatterns %>% sample_data()

R <- gen_r_matrix(X = tax, A = dummy)

human <- get_variable(GlobalPatterns, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")
sample_data(GlobalPatterns)$human <- factor(human)


rownames(R) <- rownames(glom)


glom_phyloseq <- phyloseq(otu_table(unclass(clr(glom)), taxa_are_rows = F), 
                          sample_data(sample_data(GlobalPatterns)))

R_phyloseq <- phyloseq(otu_table(R, taxa_are_rows = F), 
                       sample_data(sample_data(GlobalPatterns)))

glom_dist <- ordinate(glom_phyloseq, method = "MDS", distance = "euclidean") 
R_dist <- ordinate(R_phyloseq, method = "MDS", distance = "euclidean")

p1 <- plot_ordination(glom_phyloseq, glom_dist, color = "SampleType", shape = "human") + 
  theme(legend.position = "None") + labs(title = "Summation") + geom_point(size = 4)
p2 <- plot_ordination(R_phyloseq, R_dist, color = "SampleType", shape = "human") +
  labs(title = "Modified ILR") + geom_point(size = 4)



ggsave(p1+p2, filename = "prelim_results.png", dpi = 300, 
       width = 8, height = 9)


