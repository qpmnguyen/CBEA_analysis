library(tidyverse)
library(furrr)
library(tictoc)
library(ggsci)
library(progressr)
library(patchwork)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")




data <- zinb_simulation(n_samp = 300, b_spar = 0.2, b_rho = 0.1, eff_size = 3, n_tax = 5000, n_sets = 100, 
                        prop_set_inflate = 0.1, parallel = T, n_inflate = 50, prop_inflate = 0.5, samp_prop = 0.5)

gsva_scores <- generate_alt_scores(X = data$X, A = data$A, method = "ssgsea", preprocess = T, pcount = 1)



shuffle <- data$X[,sample(1:ncol(data$X),replace = F)]
A <- as.matrix(data$A[,51])

shuffled_dat <- simple_cilr(as.matrix(shuffle), A = A, preprocess = T, pcount = 1)

test <- ecdf(shuffled_dat)

test(cilr[151,51])
test(cilr[1,51])




cilr <- simple_cilr(X = data$X, A = data$A, preprocess = T, pcount = 1, method = "random")
heatmap(gsva_scores)
heatmap(cilr)
mod <- get_diff_ab(cilr, labels = data$label, method = "wilcox")
which(mod < 0.05)

eval <- cilr_eval(cilr, resample = F)
