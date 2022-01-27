library(phyloseq)
library(tidyverse)
library(CBEA)
# loading functions from different scripts

source("R/simulations.R")
source("R/utils.R")

# define simulation grid
generate_grid <- function(eval) {
    if (eval == "fdr") {
        grid <- list(
            n_samp = 10000,
            n_tax = 2000,
            spar = c(0.2, 0.4, 0.6),
            s_rho = c(0, 0.2, 0.5),
            n_inflate = c(50, 100, 150),
            eff_size = 1,
            samp_prop = 1,
            n_sets = 1,
            vary_params = FALSE,
            prop_set_inflate = 1,
            method = "normal"
        )
    } else if (eval == "pwr") {
        grid <- list(
            n_samp = 10000,
            n_tax = 2000,
            spar = c(0.2, 0.4, 0.6),
            s_rho = c(0, 0.2, 0.5),
            n_inflate = 100,
            eff_size = c(1.5, 2, 3),
            samp_prop = 1,
            n_sets = 1,
            vary_params = FALSE,
            prop_set_inflate = 1,
            method = "normal"
        )
    } else if (eval == "auc") {
        grid <- list(
            n_samp = 10000,
            n_tax = 2000,
            spar = c(0.2, 0.4, 0.6),
            s_rho = c(0, 0.2, 0.5),
            n_inflate = 100,
            eff_size = c(1.5, 2, 3),
            samp_prop = 0.5,
            n_sets = 1,
            vary_params = TRUE,
            prop_set_inflate = 1,
            method = "normal"
        )
    } else {
        stop("Unrecognized eval")
    }
    sim_settings <- cross_df(grid)
    sim_settings$id <- seq(1, nrow(sim_settings))
    return(sim_settings)
}
