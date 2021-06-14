library(teaR)
library(tidyverse)
library(phyloseq)
library(bench)
source("R/utils.R")
source("R/simulations.R")

grid <- tibble(
    n_samp = c(1000, 2000, 3000, 1000, 1000, 1000),
    n_tax = c(500, 500, 500, 500, 2500, 5000),
    eval = c("samp", "samp", "samp", "tax", "tax", "tax")
)

A2list <- function(A){
    sets <- colnames(A)
    set_list <- map(sets, ~{
        rownames(A)[which(A[,.x] == 1)] %>% as.vector()
    })
    names(set_list) <- sets
    return(set_list)
}

add_pseudocount <- function(data){
    data$X <- data$X + 1 
    return(data)
}

grid <- grid %>% rowwise() %>% 
    mutate(data = list(zinb_simulation(n_samp = n_samp, n_tax = n_tax, n_sets = n_tax/50, 
                                  spar = 0.3, s_rho = 0.2, eff_size = 3, b_rho = 0)))
grid <- crossing(grid, distr = c("norm", "mnorm"), adj = c(T,F)) 

grid <- grid %>% mutate(data = map(data, add_pseudocount))

grid <- grid %>% rowwise() %>% 
    mutate(eval = list(mark(cilr(ab_tab = data$X, set_list = A2list(data$A), distr = distr, 
                                       adj = adj)))) 
    
saveRDS(file = "data/performance_grid.rds", grid)

grid <- grid %>% unnest(time)
ggplot(grid %>% filter(n_tax == 500), aes(x = n_samp, y = time, col = distr, shape = adj)) + geom_point() + geom_line()
