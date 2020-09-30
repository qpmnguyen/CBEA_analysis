library(tidyverse)
library(furrr)
library(progressr)
source("R/simulations.R")

setwd()

pwr_sim <- create_parameters(list(
  rep = seq(1,100),
  b_spar = c(0.2, 0.5, 0.8),
  b_rho = c(0.1, 0.2, 0.5),
  prop_inflate = c(0.5, 0.8),
  eff_size = c(2,4,6)
))


plan(sequential)
tic()
opt <- furrr_options(seed = T)
with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$sim <- furrr::future_map(pwr_sim$param, ~{
    p()
    zinb_simulation(n_samp = 300, b_spar = .x$b_spar, b_rho = .x$b_rho, 
                    eff_size = .x$eff_size, n_inflate = 50, rho_ratio = 1, n_tax = 5000, 
                    n_sets = 100, prop_set_inflate = 0.5, prop_inflate = .x$prop_inflate, parallel = T)
  }, .options = opt)
})
toc()
plan(sequential)
