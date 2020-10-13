library(qs)
library(phyloseq)
library(tidyverse)
library(optparse)
library(furrr)
library(glue)
library(optparse)
source("../R/cilr.R")
source("../R/simulations.R")
source("../R/utils.R")
source("../R/diff_abundance.R")

option_list <- list(
  make_option("--ncores", type="integer", help="Number of cores to run parallel jobs on")
)

opt <- parse_args(OptionParser(option_list=option_list))

data <- qread(file = "../data/hmp_stool_16S.qs")
labels <- matrix(nrow = nsamples(data), ncol = 1000)
for (i in 1:ncol(labels)){
  labels[,i] <- rbinom(nsamples(data),1,0.5)
}

methods <- c("cilr_wilcox", "cilr_welch", "gsva_wilcox", "gsva_welch", 
                                 "deseq2", "corncob")

for (i in 1:length(methods)){
  print(methods[i])
  if (file.exists(glue("./type_i_16S_{method}.qs", method = methods[i]))){
    print("Moving to the next iteration")
    next
  } else {
    plan(multiprocess, workers = opt$ncores)
    results <- future_map_dfc(1:ncol(labels), .f = ~{
      temp <- data
      sample_data(temp)$group <- labels[,.x]
      physeq_eval(temp, method = methods[i], agg_level = "GENUS")
    }, .progress = T)
    name <- glue("./type_i_16S_{method}.qs", method = methods[i])
    qsave(results, file = name)
    plan(sequential)
  }
}



