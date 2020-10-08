# Script to get power at the single sample level. 
library(tidyverse)
library(furrr)
library(tictoc)
library(MASS)
library(fst)
library(progressr)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")


# object file is large and not included. 
parameters <- readRDS(file = "objects/pwr_sim/parameters.rds")

tic()
plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_cilr <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- readRDS(file = glue("objects/pwr_sim/simulation_{.x}.rds"))
    simple_cilr(X = data$X, A = data$A, preprocess = T, pcount = 1)
  })
})
plan(sequential)
toc()

tic()
plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$label_wc <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- readRDS(file = glue("objects/pwr_sim/simulation_{.x}.rds"))
    wc_test(X = data$X, A = data$A, thresh = 0.05, preprocess = T, pcount = 1, alt = "greater")
  })
})
plan(sequential)
toc()



parameters$label_cilr_raw <- map(parameters$scores_cilr, .f = ~cilr_eval(scores = .x, alt = "greater", 
                                                                         thresh = 0.05, resample = F, 
                                                                         return = "sig"))
opt <- furrr_options(seed = T)
plan(multisession, workers = round(availableCores()/2,0))
parameters$label_cilr_norm <- future_map(1:nrow(parameters), .f = ~{
  data <- readRDS(file = glue("objects/pwr_sim/simulation_{i}.rds",i = .x))
  cilr_eval(scores = parameters$scores_cilr[[.x]], 
            distr = "norm", alt = "greater", thresh = 0.05, resample = T, 
            X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)

plan(multisession, workers = round(availableCores()/2,0))
parameters$label_cilr_t <- future_map(1:nrow(parameters), .f = ~{
  data <- readRDS(file = glue("objects/pwr_sim/simulation_{i}.rds", i = .x))
  cilr_eval(scores = parameters$scores_cilr[[.x]], 
            distr = "t", alt = "greater", thresh = 0.05, resample = T, 
            X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)



# generating fdr comparable statistics 
parameters$fdr_cilr_raw <- map(parameters$label_cilr_raw, .f = ~calculate_statistic(eval = "pwr", pred = .x, 
                                                                                    true = rep(1, length(.x))))
parameters$fdr_cilr_norm <- map(parameters$label_cilr_norm, .f = ~calculate_statistic(eval = "pwr", pred = .x, 
                                                                                      true = rep(1, length(.x))))
parameters$fdr_wc <- map(parameters$label_wc, .f = ~calculate_statistic(eval = "pwr", pred = .x, 
                                                                        true = rep(1, length(.x))))
parameters$fdr_cilr_t <- map(parameters$label_cilr_t, .f = ~calculate_statistic(eval = "pwr", pred = .x,
                                                                                true = rep(1, length(.x))))

write_fst(parameters, path = "objects/pwr_ss_eval.fst")


                                                                