library(tidyverse)
library(compositions)
library(phyloseq)
library(patchwork)
library(vegan)

# loading default data 
data("GlobalPatterns") 

# normal processing and then transform to composition
proc_dat <- GlobalPatterns %>% filter_taxa(function(x) sum(x > 0) >= 0.1 * length(x), TRUE) %>% 
  filter_taxa(function(x) mean(x) > 0.005/100, TRUE)  %>% transform_sample_counts(function(x) x + 1)
glom <- proc_dat %>% tax_glom(taxrank =  "Family") # aggregate to family level 

# after aggregation - do centered log ratio transform to get aichison distances 
glom_tax <- otu_table(glom) %>% as("matrix") %>% t() %>% clr() %>% unclass()
proc_tax <- otu_table(proc_dat) %>% as("matrix") %>% t() %>% clr() %>% unclass()

# re-add back to phyloseq object
otu_table(proc_dat) <- otu_table(proc_tax, taxa_are_rows = F)
otu_table(glom) <- otu_table(glom_tax, taxa_are_rows = F)

# calculate distances (euclidean under clr is equivalent to aichison)
not_glom_dist <- distance(proc_dat, method = "euclidean")
glom_dist <- distance(glom, method = "euclidean")

# ordinations (MDS)
ord_glom <- ordinate(glom, "MDS", glom_dist)
ord_not_glom <- ordinate(proc_dat, "MDS", not_glom_dist)

# plotting individual plots
plt1 <- plot_ordination(glom, ord_glom, color = "SampleType") + theme_bw() + geom_point(size = 4) +
  theme(legend.position = "None") + labs(title = "Aggregated Composition")
plt2 <- plot_ordination(proc_dat, ord_not_glom, color = "SampleType") + theme_bw() + geom_point(size = 4) +
  labs(title = "Original Composition")
plot(plt1 + plt2)

# mantel test 
mantel(glom_dist, not_glom_dist)

# procrustes 
procrust <- protest(ord_not_glom$vectors[,c(1,2)], ord_glom$vectors[,c(1,2)])
procrust
plot(procrust)

# procrustes superimposition plots 
plot1 <- ggplot() + geom_point(aes(x = V1, y = V2, color = "blue"), size = 4, data = as.data.frame(procrust$Yrot)) + 
  geom_point(aes(x = Axis.1, y = Axis.2, color = "red"), size = 4, data = as.data.frame(procrust$X)) + 
  labs(x = "Axis.1", y = "Axis.2", title = "Superimposed Procrustes Plot", subtitle = "p = 0.001*") + 
  theme_bw() + 
  scale_color_discrete(name = "Data type", labels = c("red" = "Original Composition", 
                                                       "blue" = "Aggregated Composition")) 

# scaled euclidean distances 
plot2 <- ggplot() + geom_histogram(aes(x = scale(as.vector(not_glom_dist)), fill = "red")) + 
  geom_histogram(aes(x = scale(as.vector(glom_dist)), fill = "blue")) + 
  scale_fill_discrete(name = "Data type", labels = c("red" = "Original Composition",
                                                     "blue" = "Aggregated Composition")) + 
  labs(x = "Scaled euclidean distances", y = "Counts") + theme_bw() + theme(legend.position = "None")

# raw euclidean distances 
plot3 <- ggplot() + geom_histogram(aes(x = as.vector(not_glom_dist), fill = "red")) + 
  geom_histogram(aes(x = as.vector(glom_dist), fill = "blue")) + 
  scale_fill_discrete(name = "Data type", labels = c("red" = "Original Composition",
                                                     "blue" = "Aggregated Composition")) + 
  labs(x = "Raw euclidean distances", y = "Counts") + theme_bw() + theme(legend.position = "None")

aggregation_distance_plot <- (plot2/plot3) | (plot1/(plt1 + plt2)) + plot_annotation(tag_levels = "A")
ggsave(aggregation_distance_plot, filename = "docs/aggregation_distance_plot.png", dpi = 300, width = 15, height = 10)
