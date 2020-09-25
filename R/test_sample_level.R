# Sample level inference comparing to standard wilcoxon rank-sum test 
library(tidyverse)
library(GSVA)
library(ROCR)
library(furrr)
library(tictoc)
source("R/cilr.R")
source("R/simulations.R")

params <- cross_df(list(
  rep = seq(1,10),
  n_samp = 200, 
  b_spar = c(0.2, 0.4, 0.6, 0.8), 
  b_rho = 0.2, 
  rho_ratio = c(1,2,3), 
  eff_size = c(1, 1.2, 1.4, 1.6, 1.8, 2), 
  n_inflate = c(50,100,150,200)
))
params <- params[,-1]

params$sim <- pmap(params, ~zinb_simulation(n_samp = ..1, b_spar = ..2, b_rho = ..3,  
                              eff_size = ..4, rho_ratio = ..5, 
                              n_inflate = ..6)) 
sim_cache <- saveRDS(params, file = "objects/parameters.rds")

# False Discovery Rate
fdr_sim <- params %>% filter(eff_size == 1)

tic()
cilr_scores <- map(fdr_sim$sim, .f = ~simple_cilr(X = .x$X, A = .x$A, preprocess = T, pcount = 0.5, transform = "prop", 
                                   abs = F))
toc()

fdr_sim$cilr_scores <- cilr_scores
fdr_sim$wc_test <- 






