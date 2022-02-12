library(CBEA)
library(tidyverse)
library(bench)
library(furrr)
library(future)
source("R/simulations.R")
source("R/plot_utils.R")

# benchmarking the same data set 
df <- zinb_simulation(n_samp = 500, spar = 0, s_rho = 0.2, eff_size = 3, 
                      b_rho = 0.1, n_tax = 800,
                      n_inflate = 20, n_sets = 40, prop_set_inflate = 0.5)
results <- bench::press(adj = c(TRUE, FALSE), distr = c("norm", "mnorm"), n_perm = c(50,100,200), {
    bench::mark(
        cbea(df$obj, df$set, abund_values = "Counts", distr = distr, adj = adj, output = "sig", 
             n_perm = n_perm)
    )    
})
res_plot <- results %>% dplyr::select(adj, distr, n_perm, median)
saveRDS(res_plot, "output/runtime_standard.rds")
