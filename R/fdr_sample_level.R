# Sample level inference comparing to standard wilcoxon rank-sum test 
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

#TODO Wrap data generating in a nice function with a progress bar!


# False Discovery Rate
fdr_sim <- cross_df(list(
  b_spar = c(0.2, 0.4, 0.6, 0.8),
  b_rho = c(0.1, 0.2, 0.5),
  n_inflate = c(50,100,150,200)
))

fdr_sim <- fdr_sim %>% mutate(id = seq(1:nrow(fdr_sim))) %>% group_by(id) %>% nest()


plan(multiprocess, workers = round(availableCores()/2,0))
tic()
opt <- furrr_options(seed = T)
with_progress({
  p <- progressor(steps = nrow(fdr_sim))
  fdr_sim$sim <- furrr::future_map(fdr_sim$data, ~{
    p()
    zinb_simulation(n_samp = 10000, b_spar = .x$b_spar, b_rho = .x$b_rho, 
                    eff_size = 1, n_inflate = .x$n_inflate, rho_ratio = 1)
  }, .options = opt)
})
toc()
plan(sequential)

# generating scores 
tic()
plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(fdr_sim))
  fdr_sim$scores_cilr <- future_map(fdr_sim$sim, .f = ~{
    p()
    simple_cilr(X = .x$X, A = .x$A, preprocess = T, pcount = 1, transform = "prop", abs = F)
  })
})
plan(sequential)
toc()
tic()
plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(fdr_sim))
  fdr_sim$label_wc <- map(fdr_sim$sim,.f = ~{
    p()
    wc_test(X = .x$X, A = .x$A, thresh = 0.05, preprocess = T, pcount = 1, alt = "greater")
  })
})
plan(sequential)
toc()
# generating labels with threshold 
fdr_sim$label_cilr_raw <- map(fdr_sim$scores_cilr, .f = ~cilr_eval(scores = .x, alt = "greater", thresh = 0.05, resample = F))
fdr_sim$label_cilr_norm <- map2(fdr_sim$scores_cilr, fdr_sim$sim, .f = ~cilr_eval(scores = .x, distr = "norm", alt = "greater", thresh = 0.05, resample = T, 
                                                                        X = .y$X, A = .y$A))
fdr_sim$label_cilr_t <- map2(fdr_sim$scores_cilr, fdr_sim$sim, .f = ~cilr_eval(scores = .x, distr = "t", alt = "greater", thresh = 0.05, resample = T, 
                                                                               X = .y$X, A = .y$A))

# generating fdr comparable statistics 
fdr_sim$fdr_cilr_raw <- map(fdr_sim$label_cilr_raw, .f = ~calculate_statistic(eval = "fdr", pred = .x))
fdr_sim$fdr_cilr_norm <- map(fdr_sim$label_cilr_norm, .f = ~calculate_statistic(eval = "fdr", pred = .x))
fdr_sim$fdr_wc <- map(fdr_sim$label_wc, .f = ~calculate_statistic(eval = "fdr", pred = .x))
fdr_sim$fdr_cilr_t <- map(fdr_sim$label_cilr_t, .f = ~calculate_statistic(eval = "fdr", pred = .x))


# PLOTTING 
plot_df <- fdr_sim %>% unnest(starts_with("fdr")) %>% dplyr::select(!starts_with(c("label", "scores", "sim"))) %>%
  pivot_longer(starts_with("fdr")) %>% unnest(data) 
plot_df <- plot_df %>% dplyr::rename("Set Size" = "n_inflate", "Inter-taxa Correlation" = "b_rho")
print(head(plot_df))
#plot_df <- plot_df %>% filter(name != "fdr_cilr_resample")
plt <- ggplot(plot_df, aes(y = value, x = b_spar, col = name)) + geom_point(size = 3) +
  scale_color_nejm(labels = c("cILR normal", "cILR raw", "cILR t", "Wilcoxon Rank Sum")) + 
  geom_line() + geom_hline(yintercept = 0.05, col = "red") +  
  facet_grid(`Inter-taxa Correlation` ~ `Set Size`, scales = "free_y", labeller = label_both) +
  labs(col = "Evaluation Models", x = "Sparsity", y = "Type I error rate (N = 10,000, alpha = 0.05)") + 
  theme_bw()
plt
ggsave(plt, filename = "docs/manuscript/figures/fdr_single_sample.png", 
       dpi = 300, width = 10, height = 10)




