# this script is the new run script that takes advantage of the
# project based approaches introduced in the latest targets release 
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
    make_option(c("-a", "--analysis"), type = "character",
                help = "What type of analysis to perform. Coincide with project name",
                metavar="ANALYSIS"),
    make_option(c("-r", "--remove"), type = "logical", default=TRUE,
                help = "Restart pipeline from scratch"),
    make_option(c("-p", "--parallel"), type = "logical", default=TRUE,
                help = "Running job in parallel")
)

opt <- parse_args(OptionParser(option_list=option_list))


proj_name <- opt$analysis

match.arg(proj_name, c("data_ss"))

Sys.setenv(TAR_PROJECT = proj_name)

if (opt$remove == TRUE){
    tar_destroy()
}

if (opt$parallel == TRUE){
    tar_make_future(workers = opt$ncores)
} else {
    tar_make()
}
# targets::tar_make() # debug locally
# targets::tar_visnetwork() # visualize network to debug locally
# targets::tar_make_future(workers = opt$ncores) # run using futures on a cluster 
