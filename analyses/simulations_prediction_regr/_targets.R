library(targets)
library(tarchetypes)
library(future)

plan(multisession)
source("../simulations_prediction_functions/pred_functions.R")

set.seed(1020)

tar_option_set(error = "workspace")


sim_regr <- cross_df(list(
    type = c("regr"),
    snr = c(1.5, 2, 3),
    sat = c(0.1, 0.5),
    n_samp = 1000, 
    spar = c(0.2,0.4,0.6), 
    s_rho = c(0,0.2,0.5), 
    eff_size = 1,
    b_rho = 0, 
    n_tax = 5000, 
    n_inflate = 100, 
    n_sets = 50,
    prop_set_inflate = 0.5, 
    prop_inflate = 1, 
    samp_prop = 0.5,
    method = "normal", 
    vary_params=FALSE
))
sim_regr$id <- seq(1, nrow(sim_regr))
saveRDS(sim_regr, file = "output/simulation_grid_regr.rds")


sim_eval_grid <- tar_target(sim_eval_grid, {
    eval_settings <- cross_df(list(
        model = "cilr",
        distr = c("mnorm", "norm"),
        adj = c(TRUE, FALSE),
        output = c("zscore", "cdf")
    ))
    other <- tibble(model = c("ssgsea", "gsva"))
    eval_settings <- dplyr::bind_rows(other, eval_settings)
    eval_settings
})



regr_jobs <- tar_map(unlist = FALSE, values = sim_regr, names = c("id"),
                    tar_target(simulation_regr, {
                        sim_prediction(type = type, snr = snr, sat = sat, 
                                       n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                                        samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                                        n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                                        method = "normal")}),
                    tar_target(agg_regr, {
                        print("Currently running")
                        generate_aggregation(proc_object = simulation_regr, 
                                                method = sim_eval_grid$model, 
                                                distr = sim_eval_grid$distr, 
                                                adj = sim_eval_grid$adj, 
                                                output = sim_eval_grid$output)
                    }, pattern = map(sim_eval_grid)),
                    tar_target(eval_regr, {
                        print("Currently evaluating")
                        res <- fit_and_eval(agg_regr, nfolds = 10, task = "regression")
                        dplyr::bind_cols(id = id, sim_eval_grid, res)
                    }, pattern = map(agg_regr, sim_eval_grid))
)

combine_regr <- tar_combine(combine_regr, regr_jobs[[3]], command = dplyr::bind_rows(!!!.x))

save_regr <- tarchetypes::tar_rds(save_regr, saveRDS(combine_regr, file = "output/sim_pred_regr.rds"))

#list(grid, sim, hypo_grid, eval)
list(sim_eval_grid, regr_jobs, combine_regr, save_regr)


