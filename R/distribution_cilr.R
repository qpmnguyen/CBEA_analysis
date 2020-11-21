# Distribution of the cilr test statistic under the null 
# Quang Nguyen
# Last updated 11/21

library(tidyverse)
library(fitdistrplus)
library(glue)
library(goftest)
library(furrr)
library(parameters)
library(mixtools)
library(qs)
library(goftest)

source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

parameters <- create_parameters(list(
  rep = seq(1,10),
  spar = c(0.6),
  n_inflate = c(100),
  s_rho = c(0.0, 0.2, 0.4, 0.6)
))

# Generating data 
plan(multisession, workers = 3)
opt <- furrr_options(seed = TRUE)
parameters$data <- future_map(parameters$param, .f = ~{
  zinb_simulation(n_samp = 1e3, spar = .x$spar, s_rho = .x$s_rho, b_rho = 0, eff_size = 1, 
                  vary_params = FALSE, n_tax = 2000, 
                  n_inflate = .x$n_inflate, n_sets = 1, samp_prop = 1, method = "normal")
},.options = opt, .progress = TRUE)
plan(sequential)

parameters <- parameters %>% unnest(param)
eval_df <- cross_df(list(
  distr = c('norm', 'mnorm'),
  adj = c(TRUE, FALSE),
  rep = unique(parameters$rep)
))

parameters <- left_join(parameters, eval_df)

# perform evaluations 
plan(multisession, workers = 3)
parameters$evaluation <- future_map(1:nrow(parameters), .f = ~{
  source("R/cilr.R")
  get_fit(data = parameters$data[[.x]], adj = parameters$adj[.x], distr = parameters$distr[.x],
          maxit=1e5, maxrestarts=1e3, epsilon=1e-6)
}, .options = furrr_options(seed = TRUE), .progress = TRUE)
plan(sequential)

saveRDS(parameters, file = "cache/parameter_fit.rds")

parameters <- parameters %>% unnest(evaluation)

parameters <- parameters %>% 
    pivot_longer(c(aic, bic, ad), names_to = "eval") %>%
    mutate(value = na_if(value, Inf)) %>% 
    group_by(s_rho, distr, adj, eval) %>% 
    summarise(mean = mean(value, na.rm = TRUE), 
        upper = mean(value, na.rm = TRUE) + sd(value, na.rm = TRUE), 
        lower = mean(value, na.rm = TRUE) - sd(value, na.rm = TRUE)) 

ggplot(parameters, aes(x = s_rho, y = mean, col = distr)) + 
    geom_line(aes(linetype = adj)) + geom_point(aes(shape = adj)) + 
    facet_grid(eval~., scales = "free_y")
