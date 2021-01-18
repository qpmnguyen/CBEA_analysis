library(targets)
library(tarchetypes)
library(tidyverse)
library(optparse)

# directory for local debug 
# setwd(file.path("D:/research/teailr/analyses/data_single_sample/")) 


option_list <- list( 
    make_option(c("-c", "--ncores"), type="integer", default=5, 
                help="Number of cores",
                metavar="NUMBER"),
    make_option(c("-d", "--dir"), type = "character",
                help = "Directory to move to",
                metavar="DIR")
)

opt <- parse_args(OptionParser(option_list=option_list))

dir <- opt$dir
print(glue("Currently in directory {dir}", dir = dir))

setwd(dir)


# targets::tar_make() # debug locally
# targets::tar_visnetwork() # visualize network to debug locally
targets::tar_make_future(workers = opt$ncores) # run using futures on a cluster 
