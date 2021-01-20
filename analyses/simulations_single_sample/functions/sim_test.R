library(phyloseq)
library(tidyverse)

# loading functions from different scripts  
source("../../R/cilr.R")
source("../../R/simulations.R")
source("../../R/utils.R")

# define simulation grid
generate_grid <- function(eval){
    if (eval == "fdr"){
        grid <- list(
            n_samp = 1000,
            n_tax = 100,
            spar = c(0.2, 0.4, 0.6),
            s_rho = c(0, 0.2, 0.5),
            n_inflate = c(50,100,150),
            eff_size = 1,
            samp_prop = 1,
            n_sets = 1, 
            vary_params = FALSE,
            prop_set_inflate = 1,
            method = "normal"
        )
    } else if (eval == "pwr"){
        grid <- list(
            n_samp = 10000,
            n_tax = 2000,
            spar = c(0.2, 0.4, 0.6),
            s_rho = c(0, 0.2, 0.5),
            n_inflate = 100,
            eff_size = c(1.5,2,3),
            samp_prop = 1,
            n_sets = 1,
            vary_params = FALSE,
            prop_set_inflate = 1,
            method = "normal"
        )
    } else if (eval == "auc"){
        grid <- list(
            n_samp = 10000, 
            n_tax = 2000, 
            spar = c(0.2, 0.4, 0.6),
            s_rho = c(0, 0.2, 0.5),
            n_inflate = 100,
            eff_size = c(1.5,2,3),
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
    if (eval %in% c("fdr", "pwr")){
        eval_settings <- cross_df(list(
            model = "cilr",
            distr = c("mnorm", "norm"),
            adj = c(TRUE, FALSE),
            id = sim_settings$id
        ))
        wlcx <- tibble(id = 1:100, model = "wilcox")
        eval_settings <- full_join(eval_settings, wlcx, by = c("id", "model"))
    } else if (eval == "auc"){
        eval_settings <- cross_df(list(
            method = "cilr",
            distr = c("mnorm", "norm"),
            adj = c(TRUE, FALSE),
            output = c("zscore", "cdf"),
            id = sim_settings$id
        ))
    }
    #sim_settings <- left_join(sim_settings, eval_settings, by = "id")
    return(list(sim = sim_settings, eval = eval_settings))
}


# This function performs the analysis and evaluation 
analysis <- function(sim, model=c("wilcox", "cilr"), distr, adj, eval, output = NULL, ...){
    match.arg(model)
    if (eval == "auc" & is.null(output)){
        stop("Require output argument if evaluation is null")
    }
    X <- sim$X 
    A <- sim$A
    if (model == "wilcox"){
        scores <- wc_test(X = X, A = A, thresh = 0.05, preprocess = T, pcount = 1)
    } else {
        if (eval != "auc"){
            output <- "sig"
        }
        print(distr)
        scores <- cilr(X = X, A = A, resample = T, output = output, adj = adj, distr = distr, ..., 
                       maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
    }
    # simulation objects usually have only one set of interest 
    scores <- scores[, 1] %>% as.vector()
    evaluate <- calculate_statistic(eval = eval, pred = scores, true = sim$label)
    return(tibble(eval = eval, est = evaluate$estimate, lower = evaluate$lower, upper = evaluate$upper))
}
