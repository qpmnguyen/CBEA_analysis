library(tidyverse)
library(binom)
source("R/utils.R")
source("R/simulations.R")
source("R/cilr.R")

# Creating simple simulations ####
parameters <- create_parameters(params = list(
  rep = seq(1,50),
  spar = c(0.1,0.5,0.8),
  n_inflate = c(50,100,150)
))

wc <- vector(mode = "list", length = nrow(parameters))
cilr <- vector(mode = "list", length = nrow(parameters))

data <- vector(mode = "list", length = nrow(parameters))
for (i in 1:nrow(parameters)){
  params <- parameters$param[[i]]
  sim <- quick_sim(n_samp = 2000, spar = params$spar, n_tax = 1000, eff_size = 1, 
                   n_inflate = params$n_inflate, samp_prop = 1)
  gc()
  data[[i]] <- sim
}


for (i in 1:length(data)){
  wc[[i]] <- wc_test(X = data[[i]]$X, A = data[[i]]$A, thresh = 0.05, alt = "greater", preprocess = T, pcount = 1)
  gc()
}


for (i in 1:length(data)){
  cilr[[i]] <- simple_cilr(X = data[[i]]$X, A = data[[i]]$A, preprocess = T, pcount = 1)
  gc()
}

eval <- map(cilr, .f = ~cilr_eval(.x, alt = "greater", resample = F, thresh = 0.05))
fdr_cilr <- map_int(eval, .f = ~calculate_statistic(eval = "fdr", pred = .x))
fdr_wc <- map_int(wc, .f = ~calculate_statistic(eval = "fdr", pred = .x))

parameters <- parameters %>% unnest(param)
parameters$fdr_wc <- fdr_wc
parameters$fdr_cilr <- fdr_cilr
parameters <- parameters %>% pivot_longer(starts_with('fdr'))

parameters_2 <- parameters %>% group_by(spar, n_inflate, name) %>% mutate(value = value/2e3) %>% summarise(mean = mean(value),
                                                                                           upper = mean(value) + sd(value),
                                                                                           lower = mean(value) - sd(value)) 

parameters <- parameters %>% group_by(spar, n_inflate, name) %>% summarise(mean = binom.confint(value, 2e4, methods = "ac")$mean, 
                                                                     upper = binom.confint(value, 2e4, methods = "ac")$upper,
                                                                     lower = binom.confint(value, 2e4, methods = "ac")$lower)
plot <- ggplot(parameters_2, aes(x = as.factor(spar), y = mean, col = name, group = name)) + geom_point() + geom_line() + 
  geom_pointrange(aes(ymax = upper, ymin = lower)) + facet_grid(~n_inflate) + ylim(0.02,0.075)
plot
#ggsave(plot, filename = "docs/manuscript/figures/independent_samples_se_scales.png", dpi = 300, width = 5, height = 5)
