# Script to execute functions in pred_functions

library(tidyverse)
library(optparse)
library(furrr)
source("../R/cilr.R")
source("../R/utils.R")
source("../R/pred_functions.R")

option_list <- list(
  make_option("--ncores", type="integer", help="Number of cores going to be used"),
  make_option("--dir", type="character", help="Directory to add numbers"),
  make_option("--type", type= "character", help="Type of prediction task")
)

ncores <- opt$ncores
dir <- opt$dir
type <- opt$type
sim <- readRDS(file = glue("{dir}/parameters.rds", dir = dir))

eval_settings <- cross_df(list(
    distr = c("norm", "mnorm"),
    output = c("cdf", "zscore"),
    adj = c(TRUE, FALSE),
    id = sim$id
  )
)
sim <- left_join(sim, eval_settings)

plan(multisession, workers = ncores)
tic()
sim$eval <- future_map(1:nrow(sim), .f = ~{
    source("../R/cilr.R")
    data <- readRDS(file = glue("{dir}/simulation_{i}.rds",
                                dir = dir, i = sim$id[.x]))
    scores <- fit_prediction(data, type = type, distr = sim$distr[.x], adj = sim$adj[.x], output = sim$output[.x])
    return(scores)
}, .options = furrr_options(seed = TRUE), .progress = TRUE)
toc()
plan(sequential)

saveRDS(sim, file = glue("{dir}/pred_{type}_eval.rds", dir = dir, type = type))