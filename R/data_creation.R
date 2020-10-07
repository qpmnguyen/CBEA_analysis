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
  n_inflate = c(50,100,150)
)

print("Creating parameter list")
sim <- create_parameters(sim)

print("Getting furrr going")
tic()
plan(multicore, workers = 5)
opt <- furrr_options(seed = T)
with_progress({
  p <- progressor(steps = nrow(sim))
  sim$sim <- furrr::future_map(sim$param, ~{
    p()
    zinb_simulation(n_samp = 1000, b_spar = .x$b_spar, b_rho = .x$b_rho, 
                    eff_size = 1, n_inflate = .x$n_inflate, n_tax = 1000, method = "normal", 
                    samp_prop = 0.5)
  }, .options = opt)
})
plan(sequential)
toc()
print("Done with furrr. Saving data now...")

# Save data 
saveRDS(object = sim, file = "./parameters_fdr_sim.rds")