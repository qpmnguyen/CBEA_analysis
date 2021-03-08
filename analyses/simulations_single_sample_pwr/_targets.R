# targets for pwr analyses  
library(targets)
library(tarchetypes)
library(tidyverse)
source("../simulations_single_sample_functions/sim_test.R")

set.seed(1020)  

tar_option_set(error = "workspace")

sim_pwr <- generate_grid(eval = "pwr")  
saveRDS(sim_pwr, file = "output/simulation_grid_pwr.rds")
hypo_test_grid <- tar_target(hypo_test_grid, {
    eval_settings <- cross_df(list(
        model = "cilr",
        distr = c("mnorm", "norm"),
        adj = c(TRUE, FALSE)
    ))
    wlcx <- tibble(model = "wilcox")
    eval_settings <- full_join(eval_settings, wlcx, by = c("model"))
    eval_settings
})

pwr_jobs <- tar_map(unlist = FALSE, values = sim_pwr, names = c("id"),
                    tar_target(simulation_pwr, {
                        zinb_simulation(n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                                        samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                                        n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                                        method = "normal")}),
                    tar_target(result_pwr, {
                        print("Currently mapping")
                        res <- analysis(sim = simulation_pwr, 
                                        model = hypo_test_grid$model, 
                                        distr = hypo_test_grid$distr, 
                                        eval = "pwr", 
                                        adj = hypo_test_grid$adj)
                        dplyr::bind_cols(id = id, hypo_test_grid, res)
                    }, pattern = map(hypo_test_grid))
)

combine_pwr <- tar_combine(combine_pwr, pwr_jobs[[2]], command = dplyr::bind_rows(!!!.x))
save_pwr <- tarchetypes::tar_rds(save_pwr, saveRDS(combine_pwr, file = "output/sim_ss_pwr.rds"))

list(hypo_test_grid, pwr_jobs, combine_pwr, save_pwr)