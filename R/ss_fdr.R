library(tidyverse)
library(furrr)
library(tictoc)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

sim <- list(
    rep = 1,
    n_samp = 2000,
    n_tax = 1000,
    spar = c(0.2, 0.4, 0.6, 0.8),
    s_rho = c(0, 0.2, 0.5),
    n_inflate = c(50,100,150),
    eff_size = 1,
    samp_prop = 1,
    n_sets = 1, 
    vary_params = TRUE,
    prop_set_inflate = 1,
    method = "normal"
)

sim <- create_parameters(sim)
dir <- "objects/fdr_sim"
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
      score <- wc_test(X = data$X, A = data$A, thresh = 0.05, preprocess = T, pcount = 1, alt = "greater")
    } else {
      score <- cilr(X = data$X, A = data$A, resample = T, output = "sig", nperm = 5, distr = sim$distr[.x], 
                adj = sim$adj[.x], maxrestarts=1000, epsilon = 1e-6, maxit= 1e5)
    }
    return(calculate_statistic(eval = "fdr", pred = as.vector(score)))
}, .options = furrr_options(seed = TRUE), .progress = TRUE)
toc()
plan(sequential)

saveRDS(sim, file = "objects/fdr_sim/fdr_eval.rds")


sim <- readRDS(file = "objects/fdr_sim/parameters.rds")
data <- readRDS(file = "objects/fdr_sim/simulation_11.rds")

parm <- cross_df(list(
    distr = c('norm', 'mnorm'),
    adj = c(TRUE, FALSE)
))

plot_list <- map(1:nrow(parm),  ~ {
    scores <- cilr(data$X, data$A, resample = T, output = "pval", nperm = 5, 
                distr = parm$distr[.x], adj = parm$adj[.x], maxrestarts = 1000, epsilon=1e-6, 
                maxit=1e5)
    if (parm$adj[.x] == TRUE){
        title <- "Adjusted"
    } else {
        title <- "Unadjusted"
    }
    qplot(scores, geom = "histogram", color = I("black"), fill = I("steelblue")) + 
        labs(x = "p-values", y = "Frequency", title = title, subtitle = parm$distr[.x]) + 
        theme_bw()
})
plot_list
wrap_plots(plot_list)
