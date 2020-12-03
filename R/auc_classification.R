library(tidyverse)
library(furrr)
library(tictoc)
library(optparse)
source("../R/cilr.R")
source("../R/simulations.R")
source("../R/utils.R")

option_list <- list(
  make_option("--ncores", type="integer", help="Number of cores going to be used")
)

opt <- parse_args(OptionParser(option_list=option_list))

ncores <- opt$ncores
dir <- "auc_sim"

sim <- readRDS(file = glue("{dir}/parameters.rds", dir = dir))

eval_settings <- cross_df(list(
  distr = c("mnorm", "norm", "gsva", "ssgsea", "none"),
  output = c("cdf", "zscore", "scores"),
  adj = c(TRUE, FALSE),
  id = sim$id
))

eval_settings <- eval_settings %>% 
    filter(!(distr %in% c("gsva", "ssgsea", "none") & adj == TRUE)) %>% 
    filter(!(distr %in% c("gsva", "ssgsea", "none") & output %in% c("cdf", "zscore"))) %>% 
    filter(!(distr %in% c("norm", "mnorm") & output == "scores"))

sim <- left_join(sim, eval_settings, by = "id")

plan(multisession, workers = ncores)
tic()
sim$eval <- future_map(1:nrow(sim), .f = ~{
    source("../R/cilr.R")
    data <- readRDS(file = glue("{dir}/simulation_{i}.rds", 
                            dir = dir, i = sim$id[.x]))
    if (distr %in% c("gsva", "ssgsea")){
        scores <- generate_alt_scores(X = data$X, A = data$A, method = sim$distr[.x], 
                                    preprocess = T, transform=NULL, pcount=1)
    } else if (distr == "none"){
        scores <- cilr(X = data$X, A = data$A, resample = F, output = "cdf")        
    } else {
        scores <- cilr(X = data$X, A = data$A, resample = T, 
                output = sim$output[.x], nperm = 5, distr = sim$distr[.x], 
                adj = sim$adj[.x], maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
    }
    return(calculate_statistic(eval = "auc", pred = as.vector(scores), 
                                true = as.vector(data$label)))
}, .options = furrr_options(seed = TRUE), .progress = TRUE)
toc()
plan(sequential)

saveRDS(sim, file = glue("{dir}/auc_eval.rds", dir = dir))
