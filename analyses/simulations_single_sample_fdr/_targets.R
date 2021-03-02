library(targets)
library(tarchetypes)  
source("../simulations_single_sample_functions/sim_test.R")

set.seed(1020)  

tar_option_set(error = "workspace")
sim_fdr <- generate_grid(eval = "fdr")

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

combine_fdr <- tar_combine(combine_fdr, fdr_jobs[[2]], command = dplyr::bind_rows(!!!.x))

save_fdr <- tarchetypes::tar_rds(save_fdr, saveRDS(combine_fdr, file = "output/sim_ss_fdr.rds"))


list(hypo_test_grid, fdr_jobs,combine_fdr, save_fdr)