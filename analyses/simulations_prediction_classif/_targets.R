# Target files to perform prediction analyses 

library(targets)
library(tarchetypes)
library(future)

plan(multisession)

source("../simulations_prediction_functions/pred_functions.R")
set.seed(1020)

tar_option_set(error = "workspace", memory = "transient", garbage_collection = TRUE)


sim_classif <- cross_df(list(
    type = c("classif"),
    snr = c(1.5, 2, 3),
    sat = c(0.1, 0.5),
    n_samp = 2000, 
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
sim_classif$id <- seq(1, nrow(sim_classif))
saveRDS(sim_classif, file = "output/simulation_grid_classif.rds")

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



classif_jobs <- tar_map(unlist = FALSE, values = sim_classif, names = c("id"),
                     tar_target(simulation_classif, {
                         sim_prediction(type = type, snr = snr, sat = sat, 
                                        n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                                        samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                                        n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                                        method = "normal")}),
                     tar_target(agg_classif, {
                         print("Currently running")
                         generate_aggregation(proc_object = simulation_classif, 
                                              method = sim_eval_grid$model, 
                                              distr = sim_eval_grid$distr, 
                                              adj = sim_eval_grid$adj, 
                                              output = sim_eval_grid$output)
                     }, pattern = map(sim_eval_grid)),
                     tar_target(eval_classif, {
                         print("Currently evaluating")
                         res <- fit_and_eval(agg_classif, nfolds = 10, task = "classification")
                         dplyr::bind_cols(id = id, sim_eval_grid, res)
                     }, pattern = map(agg_classif, sim_eval_grid))
)

combine_classif <- tar_combine(combine_classif, classif_jobs[[3]], command = dplyr::bind_rows(!!!.x))

save_classif <- tarchetypes::tar_rds(save_classif, saveRDS(combine_classif, file = "output/sim_pred_classif.rds"))

#list(grid, sim, hypo_grid, eval)
list(sim_eval_grid, classif_jobs, combine_classif, save_classif)


