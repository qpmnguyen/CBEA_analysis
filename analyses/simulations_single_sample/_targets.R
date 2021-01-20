library(targets)
library(tarchetypes)
library(tidyverse)
source("functions/sim_test.R")
# Setting up simulation settings  
set.seed(2105)
# targets, similar to drake, will use this seed in order to generate the pool of seeds per target 
# ensuring reproducibility 
tar_option_set(error = "workspace")
grid <- generate_grid(eval = "fdr")

sim <- grid$sim

# evaluation grid for hypothesis testing (pwr/fdr)
hypo_test_grid <- tar_target(fdr_grid, {
    eval_settings <- cross_df(list(
        model = "cilr",
        distr = c("mnorm", "norm"),
        adj = c(TRUE, FALSE)
    ))
    wlcx <- tibble(model = "wilcox")
    eval_settings <- full_join(eval_settings, wlcx, by = c("model"))
    eval_settings
})

# evaluation grid for auc 
auc_grid <- tar_target(auc_grid,{
    eval_settings <- cross_df(list(
        model = "cilr",
        distr = c("mnorm", "norm"),
        adj = c(TRUE, FALSE),
        output = c("zscore", "cdf")
    ))
    other <- tibble(model = c("gsva", "ssgsea"))
    eval_settings <- full_join(eval_settings, other, by = "model")
    eval_settings
})


fdr_jobs <- tar_map(unlist = FALSE, values = sim, 
        tar_target(simulation, {
            zinb_simulation(n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                            samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                            n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                            method = "normal")}),
        tar_target(result, {
            fdr_grid$eval <- "fdr"
            fdr_grid$sim <- list(simulation)
            print("Currently mapping")
            res <- analysis(sim = fdr_grid$sim[[1]], model = fdr_grid$model, distr = fdr_grid$distr, 
                     eval = "fdr", adj = fdr_grid$adj)
            dplyr::bind_cols(sim, fdr_grid %>% dplyr::select(-c(eval,sim)), res %>% dplyr::select(-eval))
        }, pattern = map(fdr_grid))
)

fdr_combine <- tar_combine(combined_analysis, fdr_jobs[[2]], command = dplyr::bind_rows(!!!.x))

list(fdr_grid, fdr_jobs, fdr_combine)



