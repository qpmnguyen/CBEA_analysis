library(tidyverse)
library(furrr)
library(tictoc)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

sim <- list(
    rep = 1,
    n_samp = 2000,
    n_tax = 2000,
    spar = c(0.2, 0.4, 0.6),
    s_rho = c(0, 0.2, 0.5),
    n_inflate = c(100),
    eff_size = c(2, 3, 4),
    samp_prop = 1,
    n_sets = 1, 
    vary_params = FALSE,
    prop_set_inflate = 1,
    method = "normal"
)
sim <- create_parameters(sim)
dir <- "objects/pwr_sim"
saveRDS(sim %>% unnest(param), glue("{dir}/parameters.rds", dir = dir))

if(file.exists(glue("{dir}/simulation_1.rds", dir = dir)) == FALSE){
    plan(multisession, workers = 3)
    # plan(sequential)
    sim$sim <- future_map(1:nrow(sim), .f = ~{
    print(.x)
    param <- sim$param[[.x]]
    data <- zinb_simulation(n_samp = param$n_samp, spar = param$spar, s_rho = param$s_rho, 
                    eff_size = param$eff_size, n_inflate = param$n_inflate, n_tax = param$n_tax, b_rho = 0,  
                    method = param$method, samp_prop = param$samp_prop, prop_set_inflate = param$prop_set_inflate,
                    n_sets = param$n_sets, vary_params = param$vary_params)
    name <- paste0(dir,"/simulation_",.x,".rds")
    print(name)
    print("Saving file...")
    saveRDS(data, file = name)
    return("Done!")
    },.options = furrr_options(seed = TRUE), .progress = TRUE)
    plan(sequential)
}

eval_settings <- cross_df(list(
    distr = c("mnorm", "norm", "Wilcoxon"),
    adj = c(TRUE, FALSE),
    id = sim$id
))

eval_settings <- eval_settings %>% slice(-which(eval_settings$distr == "Wilcoxon" & eval_settings$adj == TRUE))

sim <- left_join(sim, eval_settings, by = "id")

plan(multisession, workers = 3)
tic()
sim$eval <- future_map(1:nrow(sim), .f = ~{
    source("R/cilr.R")
    data <- readRDS(file = glue("{dir}/simulation_{i}.rds", 
                            dir = dir, i = sim$id[.x]))
    if (sim$distr[.x] == "Wilcoxon"){
        score <- wc_test(X = data$X, A = data$A, thresh = 0.05, preprocess = T, pcount = 1, transform = "prop", alt = "greater")
    } else {
        score <- cilr(X = data$X, A = data$A, resample = T, 
                output = "sig", nperm = 5, distr = sim$distr[.x], 
                adj = sim$adj[.x], maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
    }
    return(calculate_statistic(eval = "pwr", pred = as.vector(score)))
}, .options = furrr_options(seed = TRUE), .progress = TRUE)
toc()
plan(sequential)

saveRDS(sim, file = glue("{dir}/pwr_eval.rds", dir = dir))
