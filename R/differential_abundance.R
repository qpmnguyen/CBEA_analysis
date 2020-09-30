library(tidyverse)
library(furrr)
library(tictoc)
library(ggsci)
library(progressr)
library(patchwork)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

pwr_sim <- create_parameters(list(
  rep = seq(1,100),
  b_spar = c(0.2, 0.4, 0.6, 0.8),
  b_rho = c(0.1, 0.2, 0.5),
  eff_size = c(2,3,4,5)
))

plan(multiprocess, workers = round(availableCores()/2,0))
tic()
opt <- furrr_options(seed = T)
with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$sim <- furrr::future_map(pwr_sim$param, ~{
    p()
    zinb_simulation(n_samp = 300, b_spar = .x$b_spar, b_rho = .x$b_rho, 
                    eff_size = .x$eff_size, n_inflate = 50, rho_ratio = 1)
  }, .options = opt)
})
toc()
plan(sequential)
# Generate scores 
tic()
plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$scores_cilr <- future_map(pwr_sim$sim, .f = ~{
    p()
    simple_cilr(X = .x$X, A = .x$A, preprocess = T, pcount = 1, transform = "prop", abs = F)
  })
})
toc()
plan(sequential)


