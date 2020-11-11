library(tidyverse)
library(tictoc)
library(MASS)
library(progressr)
library(qs)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

# Loading object files ####
parameters <- qread(file = "objects/fdr_sim/parameters.qs")

opt <- furrr_options(seed = T)
# Generating scores ####
tic()
plan(multisession, workers = 3)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_cilr <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("objects/fdr_sim/simulation_{.x}.qs"))
    simple_cilr(X = data$X, A = data$A, preprocess = T, pcount = 1, resample = F, method = "raw",
                transform = "prop")
  }, .options = opt)
})
plan(sequential)
toc()

# Getting wilcoxon rank sum test statistic  
tic()
plan(multisession, workers = 3)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$label_wc <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("objects/fdr_sim/simulation_{.x}.qs"))
    wc_test(X = data$X, A = data$A, thresh = 0.05, preprocess = T, pcount = 1, alt = "greater")
  })
})
plan(sequential)
toc()


plan(multisession, workers = 3)
parameters$label_cilr_norm <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/fdr_sim/simulation_{i}.qs",i = .x))
  cilr_eval(scores = parameters$scores_cilr[[.x]], resample = T,
            distr = "norm", alt = "greater", thresh = 0.05,
            X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)

plan(multisession, workers = 3)
parameters$label_cilr_norm_adj <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/fdr_sim/simulation_{i}.qs",i = .x))
  cilr_adj_eval(scores = parameters$scores_cilr[[.x]], 
                distr = "norm", alt = "greater", thresh = 0.05,
                X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)


plan(multisession, workers = 3)
parameters$label_cilr_mnorm <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/fdr_sim/simulation_{i}.qs", i = .x))
  cilr_eval(scores = parameters$scores_cilr[[.x]], resample = T, 
            distr = "mnorm", alt = "greater", thresh = 0.05, 
            X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)

plan(multisession, workers = 3)
parameters$label_cilr_mnorm_adj <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/fdr_sim/simulation_{i}.qs", i = .x))
  cilr_adj_eval(scores = parameters$scores_cilr[[.x]], 
                distr = "mnorm", alt = "greater", thresh = 0.05, 
                X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)


# generating fdr comparable statistics ####
parameters$fdr_cilr_norm <- map(parameters$label_cilr_norm, 
                                .f = ~calculate_statistic(eval = "fdr", pred = .x))
parameters$fdr_cilr_norm_adj <- map(parameters$label_cilr_norm_adj, 
                                .f = ~calculate_statistic(eval = "fdr", pred = .x))
parameters$fdr_wc <- map(parameters$label_wc, .f = ~calculate_statistic(eval = "fdr", pred = .x))
parameters$fdr_cilr_mnorm <- map(parameters$label_cilr_mnorm, 
                                 .f = ~calculate_statistic(eval = "fdr", pred = .x))
parameters$fdr_cilr_mnorm_adj <- map(parameters$label_cilr_mnorm_adj, 
                                 .f = ~calculate_statistic(eval = "fdr", pred = .x))

qs::qsave(parameters, "objects/fdr_ss_eval.qs")
