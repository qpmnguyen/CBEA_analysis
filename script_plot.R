options(tidyverse.quiet = TRUE)
library(targets)
library(tarchetypes)
library(tidyverse)


list(
    tar_rds(data_fdr_ss, ~{
        test <- readRDS("analyses/data_single_sample/output/fdr_randset_comparison.rds")
    })
    
    
)
