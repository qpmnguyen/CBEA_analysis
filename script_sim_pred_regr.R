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

regr_jobs <- tar_map(unlist = FALSE, values = sim_regr, names = c("id"),
                        tar_target(regr_dat, {
                            sim_prediction(type = type, snr = snr, sat = sat, 
                                           n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                                           samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                                           n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                                           method = "normal")}),
                        tar_target(generate_grid, {
                            settings <- cross_df(list(
                                models = c("cbea"),
                                distr = c("norm", "mnorm"),
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
                        tar_target(proc_regr, proc_sim(regr_dat)),
                        tar_target(agg_regr, {
                            print("Currently running")
                            print(proc_regr)
                            generate_aggregation(proc_regr$physeq, proc_regr$set, 
                                                 method = generate_grid$models, 
                                                 distr = generate_grid$distr, 
                                                 adj = generate_grid$adj, 
                                                 output = generate_grid$output)
                        }, pattern = map(generate_grid)),
                        tar_target(fit_regr, {
                            print("Currently evaluating")
                            fit_and_eval(agg_regr$aug_df, nfolds = 10, task = "regression", unbal_class = FALSE)
                        }, pattern = map(agg_regr)),
                        tar_target(eval_regr, {
                            # evaluation_code_here return tibble to bind rows
                            tibble(
                                id = id,
                                estimate = fit_regr %>% dplyr::filter(.metric == "rmse") %>% pull(mean),
                                stderr = fit_regr %>% dplyr::filter(.metric == "rmse") %>% pull(std_err),
                                distr = generate_grid$distr, 
                                output = generate_grid$output, 
                                adj = generate_grid$adj, 
                                models = generate_grid$models
                            )
                        }, pattern = map(fit_regr, generate_grid))
)

summary_regr <- tar_combine(combine_regr, regr_jobs[[6]], command = dplyr::bind_rows(!!!.x))

rds_regr <- tarchetypes::tar_rds(save_regr, saveRDS(combine_regr, file = "output/sim_pred_regr.rds"))

#list(grid, sim, hypo_grid, eval)
list(regr_jobs, summary_regr, rds_regr)
