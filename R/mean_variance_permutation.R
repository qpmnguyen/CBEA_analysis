# Get distribution of mean, variance 
estim_meanvar <- map2(fdr_sim$scores_cilr, fdr_sim$sim, .f = ~cilr_eval(scores = .x, alt = "greater", thresh = 0.05, resample = T, 
                                                                        X = .y$X, A = .y$A, return = "resample"))
estim_meanvar <- do.call(rbind,estim_meanvar) %>% as.data.frame()
colnames(estim_meanvar) <- c("mu", "sd")

fdr_eval <- fdr_sim %>% dplyr::select(!starts_with(c("label", "scores", "sim")))
fdr_eval <- cbind(fdr_eval, estim_meanvar)
meanvar <- fdr_eval %>% unnest(data)
mean <- meanvar %>% pivot_longer(c(mu,sd)) %>% group_by(b_spar, b_rho, n_inflate, name) %>% 
  summarise(median = median(value), max = max(value), min = min(value)) %>% 
  filter(name == "mu") %>% 
  ggplot(aes(x = b_spar, y = median)) + geom_pointrange(aes(ymax = max, ymin = min)) + geom_hline(yintercept  = 0, col = "red") +
  labs(y = "Mean parameter", x = "Sparsity") + theme_bw() + 
  facet_grid(b_rho ~ n_inflate, labeller = label_both)

sd <- meanvar %>% pivot_longer(c(mu,sd)) %>% group_by(b_spar, b_rho, n_inflate, name) %>% 
  summarise(median = median(value), max = max(value), min = min(value)) %>% 
  filter(name == "sd") %>% 
  ggplot(aes(x = b_spar, y = median)) + geom_pointrange(aes(ymax = max, ymin = min)) + geom_hline(yintercept  = 1, col = "red") +
  labs(y = "Standard Deviation parameter", x = "Sparsity") + theme_bw() + 
  facet_grid(b_rho ~ n_inflate, labeller = label_both)
ggsave(mean + sd, filename = "null_resampling_effect.png", dpi = 300, width = 10, height = 8)