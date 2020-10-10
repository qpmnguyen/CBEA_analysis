library(HMP16SData)
library(phyloseq)
library(tidyverse)
library(optparse)
library(furrr)
library(glue)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")
source("R/diff_abundance.R")


data <- qread(file = "data/hmp_stool_16S.qs")
labels <- matrix(nrow = nsamples(data), ncol = 1000)
for (i in 1:ncol(labels)){
  labels[,i] <- rbinom(nsamples(data),1,0.5)
}

methods <- c("cilr_wilcox", "cilr_welch", "gsva_wilcox", "gsva_welch", 
                                 "deseq2", "corncob")

for (i in 1:length(methods)){
  print(methods[i])
  plan(multiprocess, workers = availableCores()/2)
  results <- future_map_dfc(1:ncol(labels), .f = ~{
    temp <- data
    sample_data(temp)$group <- labels[,.x]
    physeq_eval(temp, method = methods[i], agg_level = "GENUS")
  }, .progress = T)
  name <- glue("objects/type_i_16S_{method}.qs", method = methods[i])
  qsave(results, file = name)
  plan(sequential)
}



