# Analyzing the impact of sample size to power and type I error  
library(tidyverse)
library(furrr)
library(qs)
library(progressr)
library(optparse)
source("../R/cilr.R")
source("../R/simulations.R")
source("../R/utils.R")


option_list <- list(
  make_option("--ncores", type = "integer", help="Number of workers to use for parallelization")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Generating data with different sample sizes 
sim <- list(
  rep = 1,
  n_samp = c(200, 2000, 20000),
  n_tax = 2000,
  spar = c(0.2, 0.4, 0.6, 0.8), 
  s_rho = 0,
  n_inflate = 100,
  eff_size = c(1.5,2,2.5,3),
  n_sets = 1,
  vary_params = FALSE,
  samp_prop = 1,
  prop_set_inflate = 1,
  method = "normal"
)
sim <- create_parameters(sim)
dir <- "samp_size"
# if these files do not exist (else rerun again)
if (!file.exists(glue("{dir}/simulation_1.qs", dir = dir))){
  plan(multisession, workers = opt$ncores)
  # plan(sequential)
  sim$sim <- future_map(1:nrow(sim), .f = ~{
    print(.x)
    param <- sim$param[[.x]]
    data <- zinb_simulation(n_samp = param$n_samp, spar = param$spar, s_rho = param$s_rho, 
                            eff_size = param$eff_size, n_inflate = param$n_inflate, n_tax = param$n_tax, b_rho = 0,  
                            method = param$method, samp_prop = param$samp_prop, prop_set_inflate = param$prop_set_inflate,
                            n_sets = param$n_sets, parameters = param_file, vary_params = param$vary_params)
    name <- paste0(dir,"/simulation_",.x,".qs")
    print(name)
    print("Saving file...")
    qsave(data, file = name)
    return("Done!")
  },.options = furrr_options(seed = TRUE), .progress = TRUE)
  plan(sequential)
}

eval_settings <- cross_df(list(
  distr = c("mnorm", "norm"),
  adj = c(TRUE, FALSE),
  id = sim$id
))

sim <- left_join(sim, eval_settings, by = "id")
plan(multisession, workers = 3)
sim$pwr <- future_map(1:nrow(sim), .f = ~{
  dat <- qread(file = glue("{dir}/simulation_{i}.qs", dir = dir, i = sim$id[.x]))
  score <- cilr(X = dat$X, A = dat$A, resample = T, output = "sig", nperm = 5, distr = sim$distr[.x], 
                adj = sim$adj[.x], maxrestarts=1000, epsilon = 1e-6, maxit= 1e5) 
  return(calculate_statistic(eval = "pwr", pred = score))
}, .options = furrr_options(seed = TRUE), .progress = TRUE)
plan(sequential)

qsave(sim, file = glue("{dir}/samp_eval.qs", dir = dir))
