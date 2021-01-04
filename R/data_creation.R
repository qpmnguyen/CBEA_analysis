library(tidyverse)
library(furrr)
library(progressr)
library(tictoc)
library(glue)
library(qs)
library(optparse)
source("../R/simulations.R")

option_list <- list(
  make_option("--ncores", type = "integer", help="Number of workers to use for parallelization"),
  make_option("--setting", type = "character", help = "Simulation Settings to use"),
  make_option("--parallel", type = "logical", help = "Whether to perform this in parallel or not"),
  make_option("--paramfile", type = "character", help = "Directory of the parameter file", default = NULL)
)

opt <- parse_args(OptionParser(option_list=option_list))

if (opt$setting == "fdr"){
  sim <- list(
    rep = 1,
    n_samp = 20000,
    n_tax = 2000,
    spar = c(0.2, 0.5, 0.8),
    s_rho = c(0, 0.2, 0.5),
    n_inflate = c(50,100,150),
    eff_size = 1,
    samp_prop = 1,
    n_sets = 1, 
    vary_params = FALSE,
    prop_set_inflate = 1,
    method = "normal"
  )
  dir <- "fdr_sim"
} else if (opt$setting == "pwr"){
  sim <- list(
    rep = 1,
    n_samp = 20000,
    n_tax = 2000,
    spar = c(0.2, 0.5, 0.8),
    s_rho = c(0, 0.2, 0.5),
    n_inflate = 100,
    eff_size = c(1.5,2,3),
    samp_prop = 1,
    n_sets = 1,
    vary_params = FALSE,
    prop_set_inflate = 1,
    method = "normal"
  )
  dir <- "pwr_sim"
} else if (opt$setting == "auc"){
  sim <- list(
    rep = seq(1,50),
    n_samp = 2000, 
    n_tax = 2000, 
    spar = c(0.2, 0.5, 0.8),
    s_rho = c(0, 0.2, 0.5),
    n_inflate = 100,
    eff_size = c(1.5,2,3),
    samp_prop = 0.5,
    n_sets = 1,
    vary_params = TRUE,
    prop_set_inflate = 1,
    method = "normal"
  )
  dir <- "auc_sim"
} else if (opt$setting == "diff_ab"){
  sim <- list(
    rep = seq(1,10),
    n_samp = 2000,
    n_tax = 5000,
    spar = c(0.2, 0.5, 0.8), 
    s_rho = c(0, 0.2, 0.5),
    n_inflate = 100,
    eff_size = c(1.5,2,3),
    n_sets = 50,
    vary_params = TRUE,
    samp_prop = 0.5, 
    prop_set_inflate = 0.5,
    method = "compensation"
  )
  dir <- "diff_ab_sim"
}

param_file <- opt$paramfile
cores <- opt$ncores
print("Creating parameter list")
sim <- create_parameters(sim)

print("Creating folder...")
if (dir.exists(dir) == F){
    dir.create(dir)
}

saveRDS(sim %>% unnest(param), file = glue("{dir}/parameters.rds", dir = dir))

print("Getting furrr going")
tic()
if (opt$parallel == TRUE){
  plan(multicore, workers = cores)
} else {
  plan(sequential)
}
sim$sim <- future_map(1:nrow(sim), .f = ~{
  print(.x)
  param <- sim$param[[.x]]
  data <- zinb_simulation(n_samp = param$n_samp, spar = param$spar, s_rho = param$s_rho, 
                  eff_size = param$eff_size, n_inflate = param$n_inflate, n_tax = param$n_tax, b_rho = 0,  
                  method = param$method, samp_prop = param$samp_prop, prop_set_inflate = param$prop_set_inflate,
                  n_sets = param$n_sets, parameters = param_file, vary_params = param$vary_params)
  name <- paste0(dir,"/simulation_",.x,".rds")
  print(name)
  print("Saving file...")
  saveRDS(data, file = name)
  return("Done!")
},.options = furrr_options(seed = TRUE), .progress = TRUE)
plan(sequential)
toc()
print("Done with furrr...")


