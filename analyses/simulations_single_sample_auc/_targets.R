# Targets file for auc sim
library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
plan(multisession)
source("../simulations_single_sample_functions/sim_test.R")

set.seed(1020)  

tar_option_set(error = "workspace")

sim_auc <- generate_grid(eval = "auc")
saveRDS(sim_auc, file = "output/simulation_grid_auc.rds")

auc_test_grid <- tar_target(auc_test_grid, {
    eval_settings <- cross_df(list(
        model = "cilr",
        distr = c("mnorm", "norm"),
        adj = c(TRUE, FALSE),
        output = c("zscore", "cdf")
    ))
    other <- tibble(model = c("ssgsea", "gsva"))
    eval_settings <- full_join(eval_settings, other, by = "model")
    eval_settings
})

auc_jobs <- tar_map(unlist = FALSE, values = sim_auc, names = c("id"),
                    tar_target(simulation_auc, {
                        zinb_simulation(n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                                        samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                                        n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                                        method = "normal")}),
                    tar_target(result_auc, {
                        print("Currently mapping")
                        res <- analysis(sim = simulation_auc, 
                                        model = auc_test_grid$model, 
                                        distr = auc_test_grid$distr, 
                                        eval = "auc", 
                                        adj = auc_test_grid$adj, 
                                        output = auc_test_grid$output)
                        dplyr::bind_cols(id = id, auc_test_grid, res)
                    }, pattern = map(auc_test_grid))
)

combine_auc <- tar_combine(combine_auc, auc_jobs[[2]], command = dplyr::bind_rows(!!!.x))
save_auc <- tarchetypes::tar_rds(save_auc, saveRDS(combine_auc, file = "output/sim_ss_auc.rds"))


list(auc_test_grid, auc_jobs, combine_auc, save_auc)

