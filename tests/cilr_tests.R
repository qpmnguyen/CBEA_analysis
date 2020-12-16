# Title     : cILR unit tests
# Objective : Testing for cilr functions
# Created by: Quang
# Created on: 12/16/2020

library(testthat)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")
dat <- zinb_simulation(n_samp = 300, n_tax = 200, s_rho = 0, eff_size = 1)

# reference implementation for cilr
ref_implement <- function(X, idx){
    num <- mean(log(X[,idx == 1]))
    denom <- mean(log(X[,idx == 0]))
    p <- ncol(X)
	size <- sum(idx == 1)
	scale <- sqrt(size * (p - size)/p)
    return(scale*log(num/denom))
}


# testing input
