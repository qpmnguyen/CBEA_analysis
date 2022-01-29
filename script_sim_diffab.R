library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
library(glue)
library(MASS)
library(BiocSet)
library(mia)
library(future.callr)
library(future.batchtools)

source("R/functions_data_diffab.R")
source("R/simulations.R")
source("R/functions_data_ss.R")
tar_option_set(workspace_on_error = TRUE)
tar_option_set(memory = "transient", garbage_collection=TRUE)

set.seed(1020)
# SET PLAN ####
# plan(multisession)
# plan(batchtools_slurm, template = "batchtools.slurm.tmpl")
plan(callr)


# SIMULATION GRID ####
# first, define simulation grid 
sim_grid <- cross_df(list(
    rep = seq(1,10),
    n_samp = 2000, 
    spar = c(0.2,0.4,0.6), 
    s_rho = c(0,0.2,0.5), 
    eff_size = c(1, 1.5, 2, 3),
    b_rho = 0, 
    n_tax = 5000, 
    n_inflate = 50, 
    n_sets = 100,
    prop_set_inflate = 0.5, 
    prop_inflate = 1, 
    samp_prop = 0.5,
    method = "compensation", 
    vary_params=FALSE
))

sim_grid$id <- seq(1, nrow(sim_grid))
saveRDS(sim_grid, file = "output/sim_diff_ab_grid.rds")


# define an evaluation grid for cilr 
# define function that performs simulation across the defined grid  
# define function that performs the differential abundance testing across simulation grid
# define a function that performs evaluation across the differential abundance results  
analysis <- tar_map(values = sim_grid, unlist = FALSE, names = c("id"), 
                    tar_target(simulation_dat, {
                        zinb_simulation(n_samp = n_samp, spar = spar, s_rho = s_rho, eff_size = eff_size, 
                                        n_inflate = n_inflate, n_sets = n_sets, prop_set_inflate = prop_set_inflate, 
                                        method = method, samp_prop = samp_prop, vary_params = vary_params, n_tax = n_tax, 
                                        prop_inflate = prop_inflate)
                    }),
                    # processing the simulated data 
                    tar_target(proc_sim , ~{
                        rdata <- S4Vectors::DataFrame(row.names = colnames(assay(obj$obj, "Counts")),
                                                      group = obj$label)
                        SummarizedExperiment::colData(obj$obj) <- rdata
                        
                    }),
                    tar_target(eval_grid, {
                        eval_settings <- cross_df(list(
                            models = c("cbea"),
                            distr = c("mnorm", "norm"),
                            adj = c(TRUE, FALSE),
                            output = c("zscore", "cdf")
                        ))
                        other <- tibble(models = c("deseq2", "corncob"))
                        eval_settings <- dplyr::bind_rows(other, eval_settings)
                        eval_settings
                    }),
                    tar_target(analysis_res,{
                        print(eval_grid)
                        # TODO: Fix this Quang
                        diff_ab(obj = simulation_dat$obj, 
                                sets = simulation_dat$set, 
                                abund_values = "Counts", method = "models", 
                                eval = "fdr")
                    }, pattern = map(eval_grid)),
                    tar_target(evaluation_res,{
                        res <- eval_function(analysis_res, ci = FALSE)
                        tibble(id = id, res = res, eval_grid)
                    }, pattern = map(analysis_res, eval_grid))
)

combined <- tar_combine(combined_results, analysis[[5]], 
                        command = dplyr::bind_rows(!!!.x, .id = "id"))

file <- tarchetypes::tar_rds(save_file, saveRDS(combined_results, file = "output/sim_diff_ab.rds"))

list(analysis, combined, file)
