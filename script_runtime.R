library(targets)
library(tarchetypes)
library(CBEA)
library(tidyverse)
library(bench)
library(mia)
source("R/simulations.R")


values <- cross_df(list(
    n_perm = c(50, 100), 
    adj = c(TRUE, FALSE), 
    distr = c("mnorm", "norm")
))

data <- tar_target(df, {
    sim_dat <- zinb_simulation(n_samp = 500, spar = 1, s_rho = 0, eff_size = 1, 
                    n_tax = 800, n_inflate = 20, n_sets = 40, vary_params = TRUE)
    assay(sim_dat$obj, "Counts") <- assay(sim_dat$obj, "Counts") + 1
    sim_dat$obj <- mia::transformSamples(sim_dat$obj, abund_values = "Counts", method = "relabundance")
    sim_dat
})
mapping <- tar_map(values = values, unlist = FALSE, 
    tar_target(benchmarking, {
	print(df$obj)
        results <- bench::mark( 
            cbea(df$obj, df$set, abund_values = "relabundance", 
                 distr = distr, adj = adj, output = "sig", 
                 n_perm = n_perm, control = list(fix_comp = "large")),
            memory = FALSE, iterations = 1) 
        results %>% dplyr::select(adj, distr, n_perm, median)
    })
) 
combine <- tar_combine(combine_runtime, mapping[[1]], command = dplyr::bind_rows(!!!.x))
save <- tar_rds(save_runtime, saveRDS(combine_runtime, file = "output/runtime.rds"))


list(data, mapping, combine, save)
