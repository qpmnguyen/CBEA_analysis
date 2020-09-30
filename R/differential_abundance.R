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
  b_spar = c(0.2, 0.5, 0.8),
  b_rho = c(0.1, 0.2, 0.5),
  prop_inflate = c(0.5, 0.8),
  eff_size = c(2,4,6)
))


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


