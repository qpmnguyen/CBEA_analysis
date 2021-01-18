library(targets)
library(tarchetypes)
library(tidyverse)
library(optparse)
# setwd(file.path("D:/research/teailr/analyses/data_single_sample/")) 


option_list <- list( 
    make_option(c("-c", "--ncores"), type="integer", default=5, 
                help="Number of random normals to generate [default %default]",
                metavar="number")
)

opt <- parse_args(OptionParser(option_list=option_list))

# targets::tar_make() # debug locally
# targets::tar_visnetwork() # visualize network to debug locally
targets::tar_make_future(workers = opt$ncores) # run using futures on a cluster 
