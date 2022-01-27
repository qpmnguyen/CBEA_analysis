# targets file for pwr
library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
source("R/functions_simulations_ss.R")
source("R/functions_data_ss.R")

set.seed(1020)

plan(multisession)

tar_option_set(workspace_on_error = TRUE)
sim_pwr <- generate_grid(eval = "pwr")
saveRDS(sim_pwr, file = "output/simulation_grid_pwr.rds")


pwr_jobs <- tar_map(unlist = FALSE, values = sim_pwr, names = c("id"),
                    tar_target(simulation_pwr, {
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
                            fix_comp = "large"
                        ))
                        wlcx <- tibble(models = "wilcoxon")
                        eval_settings <- full_join(eval_settings, wlcx, by = "models")
                        eval_settings
                    }),
                    tar_target(result_pwr, {
                        print("Currently mapping")
                        res <- enrichment_analysis(physeq = simulation_pwr$obj, 
                                                   set = simulation_pwr$set, 
                                                   method = generate_grid$models, 
                                                   metric = "fdr", abund_values = "Counts", 
                                                   distr = generate_grid$distr, 
                                                   adj = generate_grid$adj, output = "sig", 
                                                   parametric = TRUE, n_perm = 100, 
                                                   control = list(fix_comp = generate_grid$fix_comp))
                        
                        b_test <- binom.confint(x = sum(res[,2]), n = nrow(res), methods = "ac")
                        df <- dplyr::bind_cols(id = id, generate_grid, estimate = b_test$mean, 
                                               lower = b_test$lower, 
                                               upper = b_test$upper)
                        df
                    }, pattern = map(generate_grid))
)

pwr_summary <- tar_combine(combine_pwr, pwr_jobs[[3]], command = dplyr::bind_rows(!!!.x))

pwr_save <- tarchetypes::tar_rds(save_pwr, saveRDS(combine_pwr, file = "output/sim_ss_pwr.rds"))

#list(grid, sim, hypo_grid, eval)
list(pwr_jobs, pwr_summary, pwr_save)