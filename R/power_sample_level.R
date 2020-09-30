# Script to get power at the single sample level. 
library(tidyverse)
library(furrr)
library(tictoc)
library(MASS)
library(ggsci)
library(progressr)
library(patchwork)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")


pwr_sim <- create_parameters(list(
  b_spar = c(0.2, 0.4, 0.6, 0.8),
  b_rho = c(0.1, 0.2, 0.5),
  eff_size = c(2,3,4,5)
))

plan(multiprocess, workers = round(availableCores()/2,0))
tic()
opt <- furrr_options(seed = T)
with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$sim <- furrr::future_map(pwr_sim$param, ~{
    p()
    zinb_simulation(n_samp = 10000, b_spar = .x$b_spar, b_rho = .x$b_rho, 
                    eff_size = .x$eff_size, n_inflate = 50, rho_ratio = 1, samp_prop = 1)
  }, .options = opt)
})
toc()
plan(sequential)

# generating scores 
tic()
plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$scores_cilr <- future_map(pwr_sim$sim, .f = ~{
    p()
    simple_cilr(X = .x$X, A = .x$A, preprocess = T, pcount = 1, transform = "prop", abs = F)
  })
})
with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$label_wc <- map(pwr_sim$sim,.f = ~{
    p()
    wc_test(X = .x$X, A = .x$A, thresh = 0.05, preprocess = T, pcount = 1, alt = "greater")
  })
})
plan(sequential)
toc()

pwr_sim$label_cilr_raw <- map(pwr_sim$scores_cilr, .f = ~cilr_eval(scores = .x, alt = "greater", thresh = 0.05, resample = F))
pwr_sim$label_cilr_norm <- map2(pwr_sim$scores_cilr, pwr_sim$sim, .f = ~cilr_eval(scores = .x, distr = "norm", alt = "greater", 
                                                                                  thresh = 0.05, resample = T, 
                                                                                  X = .y$X, A = .y$A))
pwr_sim$label_cilr_t <- map2(pwr_sim$scores_cilr, pwr_sim$sim, .f = ~cilr_eval(scores = .x, distr = "t", alt = "greater", 
                                                                               thresh = 0.05, resample = T, 
                                                                               X = .y$X, A = .y$A))

# generating fdr comparable statistics 
pwr_sim$pwr_cilr_raw <- map2(pwr_sim$label_cilr_raw, 
                               pwr_sim$sim, .f = ~calculate_statistic(eval = "pwr", pred = .x, true = .y$label))
pwr_sim$pwr_cilr_norm <- map2(pwr_sim$label_cilr_norm, 
                               pwr_sim$sim, .f = ~calculate_statistic(eval = "pwr", pred = .x,true = .y$label))
pwr_sim$pwr_cilr_t <- map2(pwr_sim$label_cilr_t, 
                             pwr_sim$sim, .f = ~calculate_statistic(eval = "pwr", pred = .x, true = .y$label))
pwr_sim$pwr_wc <- map2(pwr_sim$label_wc, 
                             pwr_sim$sim, .f = ~calculate_statistic(eval = "pwr", pred = .x, true = .y$label))


plot_df <- pwr_sim %>% unnest(starts_with("pwr")) %>% dplyr::select(!starts_with(c("label", "scores", "sim"))) %>%
    pivot_longer(starts_with("pwr")) %>% unnest(param) %>% group_by(name, b_spar, b_rho, eff_size) %>% 
    summarise(min = min(value, na.rm = T), max = max(value, na.rm = T), median = median(value, na.rm = T))
plot_df <- plot_df %>% dplyr::rename("Inter-taxa Correlation" = "b_rho", 
                                       "Effect Size" = "eff_size")
plt <- ggplot(plot_df, aes(y = median, x = b_spar, col = name)) + geom_point(size = 2) + 
    scale_color_nejm(labels = c("cILR norm", "cILR raw", "cILR t", "Wilcoxon Rank Sum Test")) + 
    geom_line() + geom_hline(yintercept = 0.8, col = "red") +  
    geom_errorbar(aes(ymin = min, ymax = max), width = 0.05) + 
    facet_grid(`Inter-taxa Correlation` ~ `Effect Size`, scales = "free_y", labeller = label_both) +
    labs(col = "Evaluation Models", x = "Sparsity", y = "Power (N = 10000, alpha = 0.05)") + 
    theme_bw()
plt

ggsave(plt, filename = "docs/manuscript/figures/power_single_sample.png", dpi = 300, width = 10, height = 10)


                                                                