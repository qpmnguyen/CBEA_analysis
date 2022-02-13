library(targets)
library(tarchetypes)
library(CBEA)
library(tidyverse)
library(bench)
source("R/simulations.R")


values <- cross_df(list(
    n_perm = c(50, 100), 
    adj = c(TRUE, FALSE), 
    distr = c("mnorm", "norm")
))

data <- tar_target(df, {
    zinb_simulation(n_samp = 500, spar = 1, s_rho = 0, eff_size = 1, 
                    n_tax = 800, n_inflate = 20, n_sets = 40)
})
mapping <- tar_map(values = values, unlist = FALSE, 
    tar_target(benchmarking, {
        df <- bench::mark( 
            cbea(df$obj, df$set, abund_values = "Counts", 
                 distr = distr, adj = adj, output = "sig", 
                 n_perm = n_perm, control = list(fix_comp = "large")),
            memory = FALSE, iterations = 1) 
        df %>% dplyr::select()
    })
) 
combine <- tar_combine(combine_runtime, mapping[[1]], command = dplyr::bind_rows(!!!.x))
save <- tar_rds(save_runtime, saveRDS(combine_runtime, file = "output/runtime.rds"))


list(data, mapping, combine, save)
