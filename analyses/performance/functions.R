library(teaR)
library(tidyverse)
library(bench)
source("../../R/utils.R")
source("../../R/simulations.R")

# convert A matrix to list 
A2list <- function(A){
    sets <- colnames(A)
    set_list <- map(sets, ~{
        rownames(A)[which(A[,.x] == 1)] %>% as.vector()
    })
    names(set_list) <- sets
    return(set_list)
}

# Adding a pseudocount 
add_pseudocount <- function(data){
    data <- data + 1 
    return(data)
}

