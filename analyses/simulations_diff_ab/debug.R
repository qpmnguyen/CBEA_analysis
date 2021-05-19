library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
library(glue)
library(MASS)


source("../data_diff_ab/functions/diff_ab_functions.R")
source("../../R/simulations.R")

sim_grid <- cross_df(list(
    rep = seq(1,10),
    n_samp = 2000, 
    spar = c(0.2,0.4,0.6), 
    s_rho = c(0,0.2,0.5), 
    eff_size = c(1, 1.5, 2, 3),
    b_rho = 0, 
    n_tax = 5000, 
    n_inflate = 100, 
    n_sets = 50,
    prop_set_inflate = 0.5, 
    prop_inflate = 1, 
    samp_prop = 0.5,
    method = "compensation", 
    vary_params=FALSE
))

sim_grid$id <- seq(1, nrow(sim_grid))

eval_settings <- cross_df(list(
    model = c("cilr_wilcox", "cilr_welch"),
    distr = c("mnorm", "norm"),
    adj = c(TRUE, FALSE),
    output = c("zscore", "cdf")
))
other <- tibble(model = c("deseq2", "corncob"))
eval_settings <- dplyr::bind_rows(other, eval_settings)



test <- sim_grid[174,]
sim_eval_grid <- eval_settings %>% filter(model == "deseq2")

dat <- zinb_simulation(n_samp = test$n_samp, spar = test$spar, s_rho = test$s_rho, eff_size = test$eff_size, 
                n_inflate = test$n_inflate, n_sets = test$n_sets, prop_set_inflate = test$prop_set_inflate, 
                method = test$method, samp_prop = test$samp_prop, vary_params = test$vary_params, n_tax = test$n_tax, 
                prop_inflate = test$prop_inflate)

transform_dat <- sim2phylo(dat)
analysis_res <- diff_ab(transform_dat, method = sim_eval_grid$model, 
    agg_level = "GENUS", data_type = "16S", prune = FALSE, 
    adj = sim_eval_grid$adj, distr = sim_eval_grid$distr, 
    output = sim_eval_grid$output, return = "sig")


res <- eval_function(analysis_res, ci = TRUE)
res
    