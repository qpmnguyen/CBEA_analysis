library(tidyverse)
library(furrr)
library(tictoc)
library(MASS)
library(progressr)
library(qs)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

parameters <- readRDS(file = "auc_sim/parameters.rds")

tic()
plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_cilr <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- readRDS(file = glue("auc_sim/simulation_{.x}"))
    simple_cilr(X = data$X, A = data$A, preprocess = T, pcount = 1)
  })
})
plan(sequential)
toc()


plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_gsva_pois <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- readRDS(file = glue("auc_sim/simulation_{.x}"))
    generate_alt_scores(X = data$X, A = data$A, method = "gsva", preprocess = T, transform=NULL, pcount=1)
  })
})
plan(sequential)

plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_gsva_gauss <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- readRDS(file = glue("auc_sim/simulation_{.x}"))
    generate_alt_scores(X = data$X, A = data$A, method = "gsva", preprocess = T, transform="clr", pcount=1)
  })
})
plan(sequential)

plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_ssgsea <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- readRDS(file = glue("objects/auc_sim/simulation_{.x}"))
    generate_alt_scores(X = data$X, A = data$A, method = "ssgsea", preprocess = T, transform=NULL, pcount=1)
  })
})
plan(sequential)


plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_prop <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- readRDS(file = glue("objects/auc_sim/simulation_{.x}"))
    generate_alt_scores(X = data$X, A = data$A, method = "prop", preprocess = T, transform="prop", pcount=1)
  })
})
plan(sequential)


scores <- parameters %>% dplyr::select(starts_with("scores"))
auc <- vector(mode = "list", ncol(scores)-1)
for (i in 2:ncol(scores)){
  new_name <- glue("auc_{name}", name = colnames(scores)[i])
  values <- map_dbl(1:nrow(scores), .f = function(.x){
    data <- readRDS(file = glue("objects/auc_sim/simulation_{idx}",idx =.x))
    stat <- calculate_statistic(eval = "auc", pred = scores[.x,i][[1]],
                                true = data$label)
    return(stat)
  })
  auc[[new_name]] <- values
}
auc <- do.call(cbind, auc)
auc <- as_tibble(auc)

parameters <- cbind(parameters, auc)

qsave(parameters, "auc_evaluation.fst")






