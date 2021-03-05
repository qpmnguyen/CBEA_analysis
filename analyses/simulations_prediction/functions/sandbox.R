library(tidyverse)
library(tidymodels)
source("../../R/cilr.R")
source("../../R/utils.R")
source("../../R/simulations.R")

sim <- sim_prediction(type = "classif")

sim$outcome %>% plogis()
sim$predictors



dat <- generate_aggregation(sim, method = "cilr")
result <- fit_and_eval(dat, task = "classification")
