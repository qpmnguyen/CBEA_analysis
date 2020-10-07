library(tidyverse)
library(furrr)
library(progressr)
library(optparse)
library(tictoc)
library(glue)
source("simulations.R")

sim <- list(
  rep = seq(1,100,1),
  b_spar = c(0.2, 0.4, 0.8),
  b_rho = c(0.1, 0.2, 0.5),
  eff_size = c(2,4,6)
)

print("Creating parameter list")
sim <- create_parameters(sim)

print("Creating folder...")
dir.create("auc_sim")

saveRDS(sim %>% unnest(param), file = "auc_sim/parameters.rds")

print("Getting furrr going")
tic()
plan(multicore, workers = 5)
opt <- furrr_options(seed = T)
with_progress({
  p <- progressor(steps = nrow(sim))
  sim$sim <- furrr::future_map(1:nrow(sim), .f = ~{
    p()
    param <- sim$param[[.x]]
    sim <- zinb_simulation(n_samp = 1000, b_spar = param$b_spar, b_rho = param$b_rho, 
                    eff_size = param$eff_size, n_inflate = 100, n_tax = 1000, method = "normal", 
                    samp_prop = 0.5)
    saveRDS(sim, file = glue("auc_sim/simulation_{i}", i = .x))
    return(sim)
  }, .options = opt)
})
plan(sequential)
toc()
print("Done with furrr. Saving data now...")

# Save data 
#saveRDS(object = sim, file = "./parameters_fdr_sim.rds")