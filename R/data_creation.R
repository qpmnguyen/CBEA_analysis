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
  make_option("--paramfile", type = "character", help = "Directory of the parameter file", default = NULL)
)

opt <- parse_args(OptionParser(option_list=option_list))

if (opt$setting == "fdr"){
  sim <- list(
    rep = 1,
    n_samp = 20000,
    n_tax = 2000,
    spar = c(0.2, 0.4, 0.6, 0.8),
    b_rho = c(0.1, 0.2, 0.5),
    n_inflate = c(50,100,150,200),
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
    spar = c(0.2, 0.4, 0.6, 0.8),
    b_rho = c(0.1, 0.2, 0.5),
    n_inflate = 50,
    eff_size = c(1.5,2,2.5,3),
    samp_prop = 1,
    n_sets = 1,
    vary_params = FALSE,
    prop_set_inflate = 1,
    method = "normal"
  )
  dir <- "pwr_sim"
} else if (opt$setting == "auc"){
  sim <- list(
    rep = seq(1,100),
    n_samp = 2000, 
    ntax = 2000, 
    spar = c(0.2, 0.4, 0.6, 0.8),
    b_rho = c(0.1, 0.2, 0.5),
    n_inflate = 100,
    eff_size = c(1.5,2,2.5,3),
    samp_prop = 0.5,
    n_sets = 1,
    vary_params = TRUE,
    prop_set_inflate = 1,
    method = "normal"
  )
  dir <- "auc_sim"
}

param_file <- opt$paramfile
cores <- opt$ncores
print("Creating parameter list")
sim <- create_parameters(sim)

print("Creating folder...")
if (dir.exists(dir) == F){
    dir.create(dir)
}

qsave(sim %>% unnest(param), file = glue("{dir}/parameters.qs", dir = dir))

print("Getting furrr going")
tic()
plan(multisession, workers = cores)
opt <- furrr_options(seed = T)

sim$sim <- furrr::future_map(1:nrow(sim), .f = ~{
  param <- sim$param[[.x]]
  data <- zinb_simulation(n_samp = param$n_samp, spar = param$spar, b_rho = param$b_rho, 
                  eff_size = param$eff_size, n_inflate = param$n_inflate, n_tax = param$n_tax, 
                  method = param$method, samp_prop = param$samp_prop, prop_set_inflate = param$prop_set_inflate,
                  n_sets = param$n_sets, parameters = param_file, vary_params = param$vary_params)
  return(data)
}, .options = opt, .progress = T)

plan(sequential)
toc()
print("Done with furrr. Splitting data into multiple files...")

# Save data 
sim_vec <- sim$sim 
for (i in 1:length(sim_vec)){
   name <- paste0(dir,"/simulation_",i,".qs")
   print(name)
   qsave(sim_vec[[i]], file = name)
}
print("Done!")


