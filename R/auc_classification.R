library(tidyverse)
library(furrr)
library(tictoc)
library(MASS)
library(progressr)
library(qs)
library(optparse)
source("../R/cilr.R")
source("../R/simulations.R")
source("../R/utils.R")


option_list <- list(
  make_option("--ncores", type="integer", help="Number of cores going to be used")
)

opt <- parse_args(OptionParser(option_list=option_list))

print("Loading parameters....")
parameters <- qread(file = "auc_sim/parameters.qs")
cores <- opt$ncores
furopt <- furrr::furrr_options(seed = TRUE)

print("Estimating GSVA Pois scores")
plan(multicore, workers = cores)
gsva_pois <- future_map(1:nrow(parameters), .f = ~{
data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
scores <- generate_alt_scores(X = data$X, A = data$A, method = "gsva", 
                              preprocess = T, transform=NULL, pcount=1)
calculate_statistic(eval = "auc", pred = scores, true = data$label)
  }, .options = furopt,  .progress = TRUE)

plan(sequential)

print("Estimating GSVA Gauss scores")
plan(multicore, workers = cores)
with_progress({
  p <- progressor(steps = nrow(parameters))
  gsva_gauss <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
    scores <- generate_alt_scores(X = data$X, A = data$A, method = "gsva", preprocess = T, 
                                  transform="clr", pcount=1)
    calculate_statistic(eval = "auc", pred = scores, true = data$label)
  }, .options = furopt, .progress = TRUE)
})
plan(sequential)

print("Estimating GSEA scores")
plan(multicore, workers = cores)
with_progress({
  p <- progressor(steps = nrow(parameters))
  ssgsea <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
    scores <- generate_alt_scores(X = data$X, A = data$A, method = "ssgsea", 
                                  preprocess = T, transform=NULL, pcount=1)
    calculate_statistic(eval = "auc", pred = scores, true = data$label)
  }, .options = furopt, .progress = TRUE)
})
plan(sequential)

print("Estimating proportional counts as scores")
plan(multicore, workers = cores)
with_progress({
  p <- progressor(steps = nrow(parameters))
  prop <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
    scores <- generate_alt_scores(X = data$X, A = data$A, method = "prop", preprocess = T, 
                                  transform="prop", pcount=1)
    calculate_statistic(eval = "auc", pred = scores, true = data$label)
  }, .options = furopt, .progress = TRUE)
})
plan(sequential)

print("Estimating cilr scores")
plan(multicore, workers = cores)

cilr_raw <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
  scores <- cilr(X = data$X, A = data$A, resample = F, nperm = 5, preprocess = T, pcount = 1, 
                 transform = "prop")
  calculate_statistic(eval = "auc", pred = scores, true = data$label)
}, .options = furopt, .progress = TRUE)

plan(sequential)

# Estimating scores 
eval_settings <- cross_df(list(
  distr = c("mnorm", "norm"),
  adj = c(TRUE, FALSE),
  output = c("cdf", "zscore"),
  rep = unique(parameters$rep)
))

parameters <- left_join(parameters, eval_settings, by = "rep")

print("Estimating cilr scores with resampling and adjustment")
plan(multicore, workers = cores)
parameters$auc <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
  scores <- cilr(X = data$X, A = data$A, resample = T, nperm = 5, preprocess = T, pcount = 1, 
                 transform = "prop", output = parameters$output[.x], 
                 distr = parameters$distr[.x], adj = parameters$adj[.x], thresh = 0.05)
  calculate_statistic(eval = "auc", pred = scores, true = data$label)
}, .options = furopt, .progress = TRUE)
plan(sequential)

output <- list(
  param = parameters,
  cilr_raw = cilr_raw,
  gsva_pois = gsva_pois,
  gsva_gauss = gsva_gauss,
  ssgsea = ssgsea,
  prop = prop
)

print("Saving files")
qsave(output, "auc_sim/auc_rank_evaluation.qs")






