library(targets)
library(tarchetypes)
library(tidyverse)
library(qs)
source("functions/sim_test.R")

set.seed(2105)
# targets, similar to drake, will use this seed in order to generate the pool of seeds per target 
# ensuring reproducibility 
tar_option_set(error = "workspace")

sim_fdr <- generate_grid(eval = "fdr")
sim_pwr <- generate_grid(eval = "pwr")
sim_auc <- generate_grid(eval = "auc")


# evaluation grid for hypothesis testing (pwr/fdr)
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

# Grid for auc models 
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


fdr_jobs <- tar_map(unlist = FALSE, values = sim_fdr, names = c("id"),
        tar_target(simulation_fdr, {
            zinb_simulation(n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                            samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                            n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                            method = "normal")}),
        tar_target(result_fdr, {
            hypo_test_grid$eval <- "fdr"
            hypo_test_grid$sim <- list(simulation_fdr)
            print("Currently mapping")
            res <- analysis(sim = hypo_test_grid$sim[[1]], model = hypo_test_grid$model, distr = hypo_test_grid$distr, 
                     eval = "fdr", adj = hypo_test_grid$adj)
            dplyr::bind_cols(sim_fdr, hypo_test_grid %>% dplyr::select(-c(eval,sim)), res %>% dplyr::select(-eval))
        }, pattern = map(hypo_test_grid))
)


pwr_jobs <- tar_map(unlist = FALSE, values = sim_pwr, names = c("id"),
                    tar_target(simulation_pwr, {
                        zinb_simulation(n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                                        samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                                        n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                                        method = "normal")}),
                    tar_target(result_pwr, {
                        hypo_test_grid$eval <- "pwr"
                        hypo_test_grid$sim <- list(simulation_pwr)
                        print("Currently mapping")
                        res <- analysis(sim = hypo_test_grid$sim[[1]], model = hypo_test_grid$model, 
                                        distr = hypo_test_grid$distr, 
                                        eval = "pwr", adj = hypo_test_grid$adj)
                        dplyr::bind_cols(sim_pwr, hypo_test_grid %>% dplyr::select(-c(eval,sim)), res %>% dplyr::select(-eval))
                    }, pattern = map(hypo_test_grid))
)

auc_jobs <- tar_map(unlist = FALSE, values = sim_auc, names = c("id"),
                    tar_target(simulation_auc, {
                        zinb_simulation(n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                                        samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                                        n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                                        method = "normal")}),
                    tar_target(result_auc, {
                        auc_test_grid$eval <- "auc"
                        auc_test_grid$sim <- list(simulation_auc)
                        print("Currently mapping")
                        res <- analysis(sim = auc_test_grid$sim[[1]], 
                                        model = auc_test_grid$model, distr = auc_test_grid$distr, 
                                        eval = "auc", adj = auc_test_grid$adj, output = auc_test_grid$output)
                        dplyr::bind_cols(sim_auc, auc_test_grid %>% dplyr::select(-c(eval,sim)), res %>% dplyr::select(-eval))
                    }, pattern = map(auc_test_grid))
)

# combine all the results into one file 
combine_fdr <- tar_combine(combine_fdr, fdr_jobs[[2]], command = dplyr::bind_rows(!!!.x))
combine_auc <- tar_combine(combine_auc, auc_jobs[[2]], command = dplyr::bind_rows(!!!.x))
combine_pwr <- tar_combine(combine_pwr, pwr_jobs[[2]], command = dplyr::bind_rows(!!!.x))

save_fdr <- tarchetypes::tar_qs(save_fdr, qsave(combine_fdr, file = "output/sim_ss_fdr.qs"))
save_auc <- tarchetypes::tar_qs(save_auc, qsave(combine_auc, file = "output/sim_ss_auc.qs"))
save_pwr <- tarchetypes::tar_qs(save_pwr, qsave(combine_pwr, file = "output/sim_ss_pwr.qs"))


list(hypo_test_grid, auc_test_grid, fdr_jobs, auc_jobs, pwr_jobs, 
     combine_fdr, combine_auc, combine_pwr, save_fdr, save_auc, save_pwr)


