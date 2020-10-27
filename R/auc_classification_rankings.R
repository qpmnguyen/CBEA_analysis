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

print("Estimating cilr scores")
tic()
plan(multisession, workers = cores)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_cilr <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
    simple_cilr(X = data$X, A = data$A, preprocess = T, pcount = 1, resample = F)
  }, .progress = TRUE)
})
plan(sequential)
toc()

print("Estimating cilr scores with zscore implementation")
tic()
plan(multisession, workers = cores)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_cilr <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
    simple_cilr(X = data$X, A = data$A, preprocess = T, pcount = 1, resample = T, method = "zscore", on_error = "NA",
                maxrestarts = 50, maxit=2000)
  }, .progress = TRUE)
})
plan(sequential)
toc()

print("Estimating cilr scores with cdf implementation")
tic()
plan(multisession, workers = cores)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_cilr <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
    simple_cilr(X = data$X, A = data$A, preprocess = T, pcount = 1, resample = T, method = "cdf", on_error = "NA",
                maxrestarts = 50, maxit=2000)
  }, .progress = TRUE)
})
plan(sequential)
toc()

print("Estimating GSVA Pois scores")
plan(multisession, workers = cores)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_gsva_pois <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
    generate_alt_scores(X = data$X, A = data$A, method = "gsva", preprocess = T, transform=NULL, pcount=1)
  }, .progress = TRUE)
})
plan(sequential)

print("Estimating GSVA Gauss scores")
plan(multisession, workers = cores)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_gsva_gauss <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
    generate_alt_scores(X = data$X, A = data$A, method = "gsva", preprocess = T, transform="clr", pcount=1)
  },.progress = TRUE)
})
plan(sequential)

print("Estimating GSEA scores")
plan(multisession, workers = cores)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_ssgsea <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
    generate_alt_scores(X = data$X, A = data$A, method = "ssgsea", preprocess = T, transform=NULL, pcount=1)
  }, .progress = TRUE)
})
plan(sequential)

print("Estimating proportional counts as scores")
plan(multisession, workers = cores)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_prop <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("auc_sim/simulation_{.x}.qs"))
    generate_alt_scores(X = data$X, A = data$A, method = "prop", preprocess = T, transform="prop", pcount=1)
  })
})
plan(sequential)


print("Gathering all results into one file...")
scores <- parameters %>% dplyr::select(starts_with("scores"))
auc <- vector(mode = "list", ncol(scores)-1)
for (i in 2:ncol(scores)){
  new_name <- glue("auc_{name}", name = colnames(scores)[i])
  values <- map_dbl(1:nrow(scores), .f = function(.x){
    data <- qread(file = glue("auc_sim/simulation_{idx}.qs",idx =.x))
    stat <- calculate_statistic(eval = "auc", pred = scores[.x,i][[1]],
                                true = data$label)
    return(stat)
  })
  auc[[new_name]] <- values
}
auc <- do.call(cbind, auc)
auc <- as_tibble(auc)

parameters <- dplyr::select(-starts_with("scores"))
parameters <- bind_cols(parameters, auc)

print("Saving files")
qsave(parameters, "auc_rank_evaluation.qs")






