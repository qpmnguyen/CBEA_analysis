library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
library(glue)
library(MASS)

source("../data_diff_ab/functions/diff_ab_functions.R")
source("../../R/simulations.R")
tar_option_set(error = "workspace", memory = "transient", garbage_collection = TRUE)
set.seed(1020)

plan(multicore)

# first, define simulation grid 
sim_grid <- cross_df(list(
    rep = seq(1,10),
    n_samp = 2000, 
    spar = c(0.2,0.4,0.6), 
    s_rho = c(0,0.2,0.5), 
    eff_size = c(1, 1.5, 2, 3),
    b_rho = 0, 
    n_tax = 5000, 
    n_inflate = 100, 
    n_sets = 50,
    prop_set_inflate = 0.5, 
    prop_inflate = 1, 
    samp_prop = 0.5,
    method = "compensation", 
    vary_params=FALSE
))

sim_grid$id <- seq(1, nrow(sim_grid))
saveRDS(sim_grid, file = "output/sim_diff_ab_grid.rds")

eval_grid <- tar_target(sim_eval_grid, {
  eval_settings <- cross_df(list(
    model = c("cilr_wilcox", "cilr_welch"),
    distr = c("mnorm", "norm"),
    adj = c(TRUE, FALSE),
    output = c("zscore", "cdf")
  ))
  other <- tibble(model = c("deseq2", "corncob"))
  eval_settings <- dplyr::bind_rows(other, eval_settings)
  eval_settings
})



# define an evaluation grid for cilr 
# define function that performs simulation across the defined grid  
# define function that performs the differential abundance testing across simulation grid
# define a function that performs evaluation across the differential abundance results  
analysis <- tar_map(values = sim_grid, unlist = TRUE, names = c("id"), 
        tar_target(simulation_dat, {
            zinb_simulation(n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                            n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                            method = method, samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                            prop_inflate = prop_inflate)
        }),
        tar_target(transform_dat, {
            sim2phylo(simulation_dat)
        }),
        tar_target(analysis_res,{
            print(sim_eval_grid)
            diff_ab(transform_dat, method = sim_eval_grid$model, 
                    agg_level = "GENUS", data_type = "16S", prune = FALSE, 
                    adj = sim_eval_grid$adj, distr = sim_eval_grid$distr, 
                    output = sim_eval_grid$output)
        }, pattern = map(sim_eval_grid)),
        tar_target(evaluation_res,{
            res <- eval_function(analysis_res)
            dplyr::bind_cols(id = id, sim_eval_grid, est = res)
        })
)

combined <- tar_combine(combined_results, analysis[[4]])
file <- tarchetypes::tar_rds(save_file, saveRDS(combined_results, file = "output/sim_diff_ab.rds"))

list(eval_grid, analysis, combined, file)
