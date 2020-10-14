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
  make_option("--data_path", type = "character", help = "Path to the data set to be loaded"),
  make_option("--ncores", type="integer", help="Number of cores to run parallel jobs on"),
  make_option("--method", type ="character", help="What methods to evaluate"),
  make_option("--label_path", type = "character", help = "Path to label file")
)

opt <- parse_args(OptionParser(option_list=option_list))

path <- opt$data_path
data_name <- strsplit(path, "/")[[1]]
dir_name <- strsplit(data_name[length(data_name)],".", fixed = T)[[1]][1]

if (!dir.exists(dir_name)){
  dir.create(dir_name)
}

method <- opt$method
data <- qread(file = path)
labels <- qread(file = opt$label_path)

plan(multiprocess, workers = opt$ncores)
results <- future_map_dfc(1:ncol(labels), .f = ~{
  temp <- data
  sample_data(temp)$group <- labels[,.x]
  physeq_eval(temp, method = method, agg_level = "GENUS")
}, .progress = T)
output_name <- glue("{dir}/type_i_{method}.qs", method = method, dir = dir_name)
qsave(results, file = output_name)
plan(sequential)



