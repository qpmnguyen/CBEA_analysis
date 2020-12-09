library(tidyverse)
library(furrr)
library(tictoc)
library(optparse)
source("../R/cilr.R")
source("../R/simulations.R")
source("../R/utils.R")


option_list <- list(
    make_option("--ncores", type="integer", help="Number of cores going to be used"),
    make_option("--eval", type="character", help="What is the evaluation criteria")
)

ncores <- opt$ncores
dir <- paste0(opt$eval, "_sim")
type <- opt$eval

sim <- readRDS(file = glue("{dir}/parameters.rds", dir = dir))
eval_settings <- cross_df(list(
    distr = c("mnorm", "norm", "Wilcoxon"),
    adj = c(TRUE, FALSE),
    id = sim$id
))
eval_settings <- eval_settings %>% slice(-which(eval_settings$distr == "Wilcoxon" & eval_settings$adj == TRUE))

sim <- left_join(sim, eval_settings, by = "id")

plan(multisession, workers = ncores)
tic()
sim$eval <- future_map(1:nrow(sim), .f = ~{
    source("../R/cilr.R")
    data <- readRDS(file = glue("{dir}/simulation_{i}.rds", 
                            dir = dir, i = sim$id[.x]))
    if (sim$distr[.x] == "Wilcoxon"){
        score <- wc_test(X = data$X, A = data$A, thresh = 0.05, preprocess = T, pcount = 1, transform = "prop", alt = "greater")
    } else {
        score <- cilr(X = data$X, A = data$A, resample = T, 
                output = "sig", nperm = 5, distr = sim$distr[.x], 
                adj = sim$adj[.x], maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
    }
    return(calculate_statistic(eval = type, pred = as.vector(score)))
}, .options = furrr_options(seed = TRUE), .progress = TRUE)
toc()
plan(sequential)

saveRDS(sim, file = glue("{dir}/{type}_eval.rds", dir = dir, type = type))