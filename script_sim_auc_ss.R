# Targets file for auc sim
library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
library(pROC)
source("R/functions_simulations_ss.R")
source("R/functions_data_ss.R")

tar_option_set(workspace_on_error = TRUE)


set.seed(1020)

plan(multisession)

sim_auc <- generate_grid(eval = "auc")
saveRDS(sim_auc, file = "output/simulation_grid_auc.rds")

auc_jobs <- tar_map(unlist = FALSE, values = sim_auc, names = c("id"),
                     tar_target(simulation_auc, {
                         zinb_simulation(n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                                         samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                                         n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                                         method = "normal")
                     }),
                     tar_target(generate_grid, {
                         eval_settings <- cross_df(list(
                             models = "cbea", 
                             distr = c("mnorm", "norm"),
                             adj = c(TRUE, FALSE),
                             output = c("zscore", "cdf"),
                             fix_comp = "large"
                         ))
                         add <- tibble(models = c("ssgsea", "gsva", "wilcoxon"))
                         eval_settings <- full_join(eval_settings, add, by = "models")
                         eval_settings
                     }),
                     tar_target(result_auc, {
                         print("Currently mapping")
                         res <- enrichment_analysis(physeq = simulation_auc$obj, 
                                                    set = simulation_auc$set, 
                                                    method = generate_grid$models, 
                                                    metric = "auc", 
                                                    abund_values = "Counts", 
                                                    distr = generate_grid$distr, 
                                                    adj = generate_grid$adj, 
                                                    output = generate_grid$output, 
                                                    parametric = TRUE, n_perm = 100, 
                                                    control = list(fix_comp = generate_grid$fix_comp))
                         auc_conf <- pROC::ci.auc(predictor = res %>% pull(2), 
                                                  response = simulation_auc$label, conf.level = 0.95, 
                                                  method = "delong")
                         
                         df <- dplyr::bind_cols(id = id, generate_grid, 
                                                estimate = auc_conf[2], 
                                                lower = auc_conf[1], 
                                                upper = auc_conf[3])
                         df
                     }, pattern = map(generate_grid))
)

auc_summary <- tar_combine(combine_auc, auc_jobs[[3]], command = dplyr::bind_rows(!!!.x))

auc_save <- tarchetypes::tar_rds(save_auc, saveRDS(combine_auc, file = "output/sim_ss_auc.rds"))

#list(grid, sim, hypo_grid, eval)
list(auc_jobs, auc_summary, auc_save)
