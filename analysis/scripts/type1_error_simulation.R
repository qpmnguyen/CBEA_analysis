library(compositions)
library(tidyverse)
library(patchwork)
library(GSVA)
library(ROCR)
source("R/simulations.R")
#source("R/utils.R")


simple_ilr <- function(tax, idx){
  tax <- unclass(acomp(tax + 1))
  n_tax <- ncol(tax)
  inset <- length(idx)
  outset <- n_tax - inset
  coef <- sqrt((inset*outset)/(inset + outset))
  num <- geometricmeanRow(tax[,idx])
  denom <- geometricmeanRow(tax[,-idx])
  result <- coef * log(num/denom)
  return(result)
}

calculate_fdr <- function(result, cutoff = 0.05){
  zthresh <- qnorm(1 - cutoff)
  label <- rep(0, length(result))
  label[result >= zthresh] <- 1
  #label[result >= zthresh | result <= -zthresh] <- 1
  fdr <- length(label[label == 1])/length(label)
  return(fdr)
}

calculate_power <-  function(result, truth_idx, cutoff = 0.05){
  zthresh <- qnorm(1 - cutoff)
  label <- rep(0, length(result))
  #label[result >= zthresh | result <= -zthresh] <- 1
  label[result >= zthresh] <- 1
  truth <- rep(0, length(result))
  truth[truth_idx] <- 1
  
  power <- length(which(truth == 1 & label == 1))/length(truth[truth == 1])
  return(power)
}

calculate_auc <- function(result, truth_idx){
  result <- abs(result) # convert all to positive
  truth <- rep(0, length(result))
  truth[truth_idx] <- 1
  mod <- performance(prediction(result, truth), "auc")
  return(mod@y.values[[1]])
}
calibration <- readRDS(file = "data/hmp_stool_calibration.rds")
top_300 <- calibration[1:300]

param <- cross_df(list(rep = seq(10), spar = seq(0.1,0.8,0.1)))
param <- param %>% mutate(data = map(spar,~dm_simulation(template = top_300, n_samp = 200, spar = .x)))
param <- param %>% mutate(inf = map(data, ~inflate_simple(data = .x, eff_size = 1, n_inflate = 50)))
param <- param %>% mutate(result = map(inf, ~simple_ilr(tax = .x$data, idx = .x$idx$tax))) 
param <- param %>% mutate(fdr = map(result, ~calculate_fdr(.x)))
param <- param %>% unnest(fdr) %>% group_by(spar) %>% 
  summarise(lower  = min(fdr), median = median(fdr), upper = max(fdr))
fdr_plot <- param %>% ggplot(aes(x = spar, y = median)) + geom_line(color = "steelblue") + 
  geom_pointrange(aes(ymin = lower, ymax =  upper), color = "steelblue") + geom_hline(yintercept = 0.05, col = "red") + 
  theme_bw() + labs(x = "Sparsity", y = "False Discovery Rate")

param <- cross_df(list(rep = seq(10), eff_size = c(2,3,4,5,6), spar = c(0.2, 0.4, 0.6, 0.8)))
param <- param %>% mutate(data = map(spar,~dm_simulation(template = top_300, n_samp = 200, spar = .x)))
param <- param %>% mutate(inf = map2(data, eff_size, ~inflate_simple(data = .x, eff_size = .y, n_inflate = 50)))
param <- param %>% mutate(result = map(inf, ~simple_ilr(tax = .x$data, idx = .x$idx$tax))) 
param <- param %>% mutate(auc = map2(result, inf, ~calculate_auc(result = .x, truth_idx = .y$idx$samp)))
param <- param %>% mutate(power = map2(result, inf, ~calculate_power(result = .x, truth_idx = .y$idx$samp)))
param <- param %>% unnest(c(power, auc)) %>% group_by(eff_size, spar) %>% 
  summarise(low_pow  = min(power), med_pow = median(power), up_pow = max(power), 
            low_auc = min(auc), med_auc = median(auc), up_auc = max(auc))
power_plot <- param %>% ggplot(aes(x = eff_size, y = med_pow)) + geom_line(col="steelblue") + facet_wrap(~spar) + 
  geom_pointrange(aes(ymin = low_pow, ymax =  up_pow), color = "steelblue") + geom_hline(yintercept = 0.8, col = "red") +
  theme_bw() + labs(x = "Effect Size", y = "Power")

auc_plot <- param %>% ggplot(aes(x = eff_size, y = med_auc)) + geom_line(col="steelblue") + facet_wrap(~spar) + 
  geom_pointrange(aes(ymin = low_auc, ymax =  up_auc), color = "steelblue") + geom_hline(yintercept = 0.8, col = "red") +
  theme_bw() + labs(x = "Effect Size", y = "AUC")

plt <- (auc_plot + power_plot)/fdr_plot
#ggsave(plot = plt, "docs/fdr_power_auc_plot.png", dpi = 300, width = 12, height = 10)

