#TODO: TEST and ADAPT this code first 
library(targets)
library(tarchetypes)
library(future)
library(mia)
library(phyloseq)
library(BiocSet)

plan(multisession)
source("R/functions_data_pred.R")
source("R/simulations.R")
set.seed(1020)
tar_option_set(workspace_on_error = TRUE)

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

classif_jobs <- tar_map(unlist = FALSE, values = sim_classif, names = c("id"),
                     tar_target(classif_dat, {
                         sim_prediction(type = type, snr = snr, sat = sat, 
                                        n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                                        samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                                        n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                                        method = "normal")}),
                     tar_target(generate_grid, {
                         settings <- cross_df(list(
                             models = c("cbea"),
                             distr = c("norm", "morm"),
                             adj = c(TRUE, FALSE),
                             output = c("zscore", "cdf")
                         ))
                         addition <- cross_df(list(
                             models = c("ssgsea", "gsva", "clr")
                         ))
                         addition_2 <- cross_df(list(
                             models = "cbea", distr = "norm", adj = FALSE, output = "raw"))
                         settings <- full_join(settings, addition)
                         settings <- full_join(settings, addition_2)
                         settings
                     }),
                     tar_target(proc_classif, proc_sim(classif_dat)),
                     tar_target(agg_classif, {
                         print("Currently running")
                         generate_aggregation(proc_classif$physeq, proc_classif$set, 
                                              method = generate_grid$models, 
                                              distr = generate_grid$distr, 
                                              adj = generate_grid$adj, 
                                              output = generate_grid$output)
                     }, pattern = map(generate_grid, index = c(3,7))),
                     tar_target(fit_classif, {
                         print("Currently evaluating")
                         fit_and_eval(agg_classif$aug_df, nfolds = 10, task = "classification")
                     }, pattern = map(agg_classif)),
                     tar_target(eval_classif, {
                         # evaluation_code_here return tibble to bind rows
                         tibble(
                             estimate = fit_classif %>% dplyr::filter(.metric == "rmse") %>% pull(mean),
                             stderr = fit_classif %>% dplyr::filter(.metric == "rmse") %>% pull(std_err),
                             distr = distr, 
                             output = output, 
                             adj = adj, 
                             models = models
                         )
                     }, pattern = map(fit_classif, generate_grid))
)

summary_classif <- tar_combine(combine_classif, classif_jobs[[6]], command = dplyr::bind_rows(!!!.x))

rds_classif <- tarchetypes::tar_rds(save_classif, saveRDS(combine_classif, file = "output/sim_pred_classif.rds"))

#list(grid, sim, hypo_grid, eval)
list(classif_jobs, summary_classif, rds_classif)


