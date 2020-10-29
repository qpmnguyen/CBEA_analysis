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


# Generating scores ####
tic()
plan(multisession, workers = 3)
with_progress({
  p <- progressor(steps = nrow(parameters))
  parameters$scores_cilr <- future_map(1:nrow(parameters), .f = ~{
    p()
    data <- qread(file = glue("objects/fdr_sim/simulation_{.x}.qs"))
    simple_cilr(X = data$X, A = data$A, preprocess = T, pcount = 1, resample = F)
  })
})
plan(sequential)
toc()

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

# Evaluate scores based on different conditions ####
parameters$label_cilr_raw <- map(parameters$scores_cilr, .f = ~cilr_eval(scores = .x, alt = "greater", 
                                                                         thresh = 0.05, resample = F, 
                                                                         return = "sig"))
opt <- furrr_options(seed = T)
plan(multisession, workers = 3)
parameters$label_cilr_norm <- future_map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/fdr_sim/simulation_{i}.qs",i = .x))
  cilr_eval(scores = parameters$scores_cilr[[.x]], 
            distr = "norm", alt = "greater", thresh = 0.05, resample = T, 
            X = data$X, A = data$A, return = "sig")
}, .options = opt, .progress = T)
plan(sequential)


parameters$label_cilr_st <- map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/fdr_sim/simulation_{i}.qs", i = .x))
  cilr_eval(scores = parameters$scores_cilr[[.x]], 
            distr = "st", alt = "greater", thresh = 0.05, resample = T, 
            X = data$X, A = data$A, return = "sig")
})



parameters$label_cilr_mnorm <- map(1:nrow(parameters), .f = ~{
  data <- qread(file = glue("objects/fdr_sim/simulation_{i}.qs", i = .x))
  cilr_eval(scores = parameters$scores_cilr[[.x]], 
            distr = "mnorm", alt = "greater", thresh = 0.05, resample = T, 
            X = data$X, A = data$A, return = "sig")
})

# theoretical_se <- function(label, ngroups = 200){
#   label <- sample(label, size = length(label), replace = F)
#   splits <- split(label, rep(1:ngroups, length.out = length(label), each = ceiling(length(label)/ngroups)))
#   prop <- vector(length = length(splits))
#   se <- vector(length = length(splits))
#   for (i in 1:length(splits)){
#     labs <- splits[[i]]
#     counts <- sum(labs == 1)
#     prop[i] <- p <- sum(labs == 1)/length(labs) # assigning prop to both p and prop[i]
#     se[i] <- sqrt((p * (1 - p))/length(labs))
#   }
#   return(mean(se))
# }

# empirical_se <- function(label, ngroups = 200){
#   label <- sample(label, size = length(label), replace = F)
#   splits <- split(label, rep(1:ngroups, length.out = length(label), each = ceiling(length(label)/ngroups)))
#   prop <- vector(length = length(splits))
#   for (i in 1:length(splits)){
#     labs <- splits[[i]]
#     counts <- sum(labs == 1)
#     prop[i] <- counts/length(labs)
#   }
#   return(sd(prop))
# }



# generating fdr comparable statistics ####
parameters$fdr_cilr_raw <- map(parameters$label_cilr_raw, .f = ~calculate_statistic(eval = "fdr", pred = .x))
parameters$fdr_cilr_norm <- map(parameters$label_cilr_norm, .f = ~calculate_statistic(eval = "fdr", pred = .x))
parameters$fdr_wc <- map(parameters$label_wc, .f = ~calculate_statistic(eval = "fdr", pred = .x))
parameters$fdr_cilr_st <- map(parameters$label_cilr_st, .f = ~calculate_statistic(eval = "fdr", pred = .x))
parameters$fdr_cilr_mnorm <- map(parameters$label_cilr_mnorm, .f = ~calculate_statistic(eval = "fdr", pred = .x))
# parameters$fdr_se_theo_cilr_raw <- do.call(rbind, map(parameters$label_cilr_raw, .f = ~theoretical_se(.x)))
# parameters$fdr_se_emp_cilr_raw <- do.call(rbind, map(parameters$label_cilr_raw, .f = ~empirical_se(.x)))
# 
# parameters$fdr_se_theo_cilr_norm <- do.call(rbind, map(parameters$label_cilr_norm, .f = ~theoretical_se(.x)))
# parameters$fdr_se_emp_cilr_norm <- do.call(rbind, map(parameters$label_cilr_norm, .f = ~empirical_se(.x)))
# 
# parameters$fdr_se_theo_wc <- do.call(rbind, map(parameters$label_wc, .f = ~theoretical_se(.x)))
# parameters$fdr_se_emp_wc <- do.call(rbind, map(parameters$label_wc, .f = ~empirical_se(.x)))
# 
# parameters$fdr_se_theo_cilr_t <- do.call(rbind, map(parameters$label_cilr_t, .f = ~theoretical_se(.x)))
# parameters$fdr_se_emp_cilr_t <- do.call(rbind, map(parameters$label_cilr_t, .f = ~empirical_se(.x)))



qs::qsave(parameters, "objects/fdr_ss_eval.qs")
