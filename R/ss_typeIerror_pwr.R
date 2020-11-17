library(tidyverse)
library(progressr)
library(qs)
source("../R/cilr.R")
source("../R/simulations.R")
source("R/utils.R")

option_list <- list(
  make_option("--ncores", type = "integer", help="Number of workers to use for parallelization"),
  make_option("--eval", type = "character", help="What is the evaluation method")
)

opt <- parse_args(OptionParser(option_list=option_list))
eval <- opt$eval
workers <- opt$workers
dir <- paste0(eval,"_sim")

# Loading object files ####
sim <- qread(file = glue("{dir}/parameters.qs", dir = dir))

# wilcoxon rank per id 
plan(multisession, workers = workers)
with_progress({
  p <- progressor(steps = nrow(sim))
  wc_eval <- future_map(1:nrow(sim), .f = ~{
    p()
    data <- qread(file = glue("{dir}/simulation_{.x}.qs", dir = dir, .x = .x))
    label <- wc_test(X = data$X, A = data$A, thresh = 0.05, preprocess = T, 
                     pcount = 1, alt = "greater")
    return(calculate_statistic(eval = eval, pred = label))
  })
})
plan(sequential)

eval_settings <- cross_df(list(
  distr = c("mnorm", "norm"),
  adj = c(TRUE, FALSE),
  id = sim$id
))

sim <- left_join(sim, eval_settings, by = "id")
plan(multisession, workers = workers)
sim$eval <- future_map(1:nrow(sim), .f = ~{
  dat <- qread(file = glue("{dir}/simulation_{i}.qs", dir = dir, i = sim$id[.x]))
  score <- cilr(X = dat$X, A = dat$A, resample = T, output = "sig", nperm = 5, distr = sim$distr[.x], 
                adj = sim$adj[.x], maxrestarts=1000, epsilon = 1e-6, maxit= 1e5) 
  return(calculate_statistic(eval = eval, pred = score))
}, .options = furrr_options(seed = TRUE), .progress = TRUE)

plan(sequential)




obj <- list(wc = wc_eval, sim = sim)

qsave(sim, file = glue("{dir}/samp_eval.qs", dir = dir))






















