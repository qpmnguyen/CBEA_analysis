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

plan(multisession)

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

# testing sim_grid 
# sim_grid <- sim_grid[1,]

saveRDS(sim_grid, file = "output/sim_diff_ab_grid.rds")

# eval_grid_mini <- tar_target(sim_eval_grid, {
#   tibble(model = c( "cilr_welch", "deseq2", "corncob"), distr = c("norm", NA, NA), 
#          adj = c(TRUE, NA, NA))
# })





# define an evaluation grid for cilr 
# define function that performs simulation across the defined grid  
# define function that performs the differential abundance testing across simulation grid
# define a function that performs evaluation across the differential abundance results  
analysis <- tar_map(values = sim_grid, unlist = FALSE, names = c("id"), 
        tar_target(simulation_dat, {
            zinb_simulation(n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                            n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                            method = method, samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                            prop_inflate = prop_inflate)
        }),
        tar_target(transform_dat, {
            sim2phylo(simulation_dat)
        }),
        tar_target(eval_grid, {
          eval_settings <- cross_df(list(
            model = c("cilr_wilcox", "cilr_welch"),
            distr = c("mnorm", "norm"),
            adj = c(TRUE, FALSE),
            output = c("zscore", "cdf")
          ))
          other <- tibble(model = c("deseq2", "corncob"))
          eval_settings <- dplyr::bind_rows(other, eval_settings)
          eval_settings
        }),
        tar_target(analysis_res,{
            print(eval_grid)
            diff_ab(transform_dat, method = eval_grid$model, 
                    agg_level = "GENUS", data_type = "16S", prune = FALSE, 
                    adj = eval_grid$adj, distr = eval_grid$distr, 
                    output = eval_grid$output, return = "sig")
        }, pattern = map(eval_grid)),
        tar_target(evaluation_res,{
            res <- eval_function(analysis_res, ci = FALSE)
            tibble(id = id, res = res, eval_grid)
        }, pattern = map(analysis_res, eval_grid))
)

combined <- tar_combine(combined_results, analysis[[5]], 
                        command = dplyr::bind_rows(!!!.x, .id = "id"))

file <- tarchetypes::tar_rds(save_file, saveRDS(combined_results, file = "output/sim_diff_ab.rds"))

list(analysis, combined, file)
