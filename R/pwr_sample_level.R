# Script to get power at the single sample level. 
library(tidyverse)
library(furrr)
library(tictoc)
library(MASS)
library(qs)
library(progressr)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")


# object file is large and not included. 
parameters <- qread(file = "objects/pwr_sim/parameters.qs")

opt <- furrr_options(seed = TRUE)
tic()
plan(multiprocess, workers = 3)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_cilr <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("objects/pwr_sim/simulation_{.x}.qs"))
    simple_cilr(X = data$X, A = data$A, preprocess = T, pcount = 1, resample = F)
  }, .options = opt)
})
plan(sequential)
toc()

tic()
plan(multiprocess, workers = 3)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$label_wc <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("objects/pwr_sim/simulation_{.x}.qs"))
    wc_test(X = data$X, A = data$A, thresh = 0.05, preprocess = T, pcount = 1, alt = "greater")
  })
})
plan(sequential)
toc()



opt <- furrr_options(seed = T)
plan(multisession, workers = 3)
parameters$label_cilr_norm <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/pwr_sim/simulation_{i}.qs",i = .x))
  cilr_eval(scores = parameters$scores_cilr[[.x]], resample = T, 
            distr = "norm", alt = "greater", thresh = 0.05, 
            X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)

opt <- furrr_options(seed = T)
plan(multisession, workers = 3)
parameters$label_cilr_norm_adj <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/pwr_sim/simulation_{i}.qs",i = .x))
  cilr_adj_eval(scores = parameters$scores_cilr[[.x]], 
                distr = "norm", alt = "greater", thresh = 0.05, 
                X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)


plan(multisession, workers = 3)
parameters$label_cilr_mnorm <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/pwr_sim/simulation_{i}.qs", i = .x))
  cilr_eval(scores = parameters$scores_cilr[[.x]], resample = T, 
            distr = "mnorm", alt = "greater", thresh = 0.05,
            X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)

plan(multisession, workers = 3)
parameters$label_cilr_mnorm_adj <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/pwr_sim/simulation_{i}.qs", i = .x))
  cilr_adj_eval(scores = parameters$scores_cilr[[.x]], 
                distr = "mnorm", alt = "greater", thresh = 0.05,
                X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)


# generating fdr comparable statistics 

parameters$pwr_cilr_norm <- map(parameters$label_cilr_norm, .f = ~calculate_statistic(eval = "pwr", pred = .x, 
                                                                                      true = rep(1, length(.x))))
parameters$pwr_cilr_norm_adj <- map(parameters$label_cilr_norm_adj, .f = ~calculate_statistic(eval = "pwr", pred = .x, 
                                                                                      true = rep(1, length(.x))))
parameters$pwr_wc <- map(parameters$label_wc, .f = ~calculate_statistic(eval = "pwr", pred = .x, 
                                                                        true = rep(1, length(.x))))
parameters$pwr_cilr_mnorm <- map(parameters$label_cilr_mnorm, .f = ~calculate_statistic(eval = "pwr", pred = .x, 
                                                                                        true = rep(1, length(.x))))
parameters$pwr_cilr_mnorm_adj <- map(parameters$label_cilr_mnorm_adj, .f = ~calculate_statistic(eval = "pwr", pred = .x, 
                                                                                        true = rep(1, length(.x))))

qsave(parameters, file = "objects/pwr_ss_eval.qs")


                                                                