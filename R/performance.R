library(tidyverse)
library(bench)
library(profvis)
library(furrr)
source("R/cilr.R")
source("R/utils.R")
source("R/simulations.R")

set.seed(1020)

parameters <- create_parameters(list(
    n_samp = c(1e2,1e3,1e4), 
    n_tax = c(1000, 2000, 3000)
))

# Generating data 
parameters$data <- map(parameters$param, .f = ~{
    zinb_simulation(n_samp = .x$n_samp, spar = 0.1, s_rho = 0.1, b_rho = 0, eff_size = 1, 
                    vary_params = FALSE, n_tax = .x$n_tax, 
                    n_inflate = 100, n_sets = .x$n_tax/100, samp_prop = 1, method = "normal")
})

parameters <- crossing(parameters, distr = c("norm", "mnorm"))

saveRDS(parameters, file = "data/perf_simulated_dsets.rds")

parameters$scores <- map2(parameters$data, parameters$distr, .f = ~{
    bench::mark(cilr(X = .x$X, A = .x$A, distr = .y, output = "cdf", 
         resample = T, pcount = 1, transform = "prop", preprocess = T, maxrestarts=1000, epsilon = 1e-06, maxit= 1e5 ))
})


parameters <- readRDS(file = "data/perf_simulated_dsets.rds")
data <- parameters$data[[1]]

profvis(cilr(X = data$X, A = data$A, distr = "norm", output = "cdf", pcount = 1, resample = T,
                 transform = "prop", 
                 preprocess = T, maxrestarts=1000, epsilon = 1e-06, maxit= 1e5 ))

profvis(cilr(X = data$X, A = data$A, distr = "mnorm", output = "cdf", pcount = 1, resample = T,
             transform = "prop", 
             preprocess = T, maxrestarts=1000, epsilon = 1e-06, maxit= 1e5 ))





         