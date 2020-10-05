# Distribution of the cilr test statistic under the null 
# Quang Nguyen
# Last updated 9/24

library(tidyverse)
library(GSVA)
library(ROCR)
library(ggsci)
library(fitdistrplus)
library(glue)
library(patchwork)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

sim <- zinb_simulation(n_samp = 100000, b_spar = 0.2, b_rho = 0.2, eff_size = 1)
A <- sim$A
X <- sim$X

scores <- as.numeric(simple_cilr(X, A, pcount = 1, transform = NULL, preprocess = T))

norm_fit <- fitdist(scores, "norm")
norm_dist <- rnorm(100000, mean = norm_fit$estimate['mean'], sd = norm_fit$estimate['sd'])
std_norm <- rnorm(100000)

summary <- cbind(scores, norm_dist, std_norm) %>% as.data.frame() %>% pivot_longer(everything())

null_plot <- ggplot(summary, aes(x = value, fill = name)) + geom_density(alpha = 0.3) + theme_bw() + 
  scale_fill_nejm(labels = c("Fitted Normal", "Null cILR scores", "Standard Normal")) + 
  labs(fill = "Distribution", x = "Values", y = "Density")

qq_norm_std <- ggplot(as.data.frame(scores), aes(sample = scores)) + 
  stat_qq(col = "steelblue") + 
  stat_qq_line(size = 1.2, alpha = 0.7, col = "orange") + 
  theme_bw() + labs(x = "Theoretical", y = "Sample", title = "QQ-plot with Standard Normal")

qq_norm_fit <- ggplot(as.data.frame(scores), aes(sample = scores)) + 
  stat_qq(col = "steelblue", distribution = qnorm, dparams = as.list(norm_fit$estimate)) + 
  stat_qq_line(size = 1.2, alpha = 0.7, col = "orange", distribution = qnorm, dparams = as.list(norm_fit$estimate)) + 
  theme_bw() + labs(x = "Theoretical", y = "Sample", 
                    title = glue("QQ-plot with Fitted Normal (mean = {mean}, sd = {sd})",
                                 mean = round(norm_fit$estimate['mean'],2), 
                                 sd = round(norm_fit$estimate['sd'],2)))

null_distr <- (qq_norm_std /qq_norm_fit)|null_plot
ggsave(null_distr, filename = "docs/manuscript/figures/null_distribution.png", dpi = 300, 
       width = 12, height = 7)


fit_test <- fitdist(scores, "t", start = list(df = 1))
ggplot(as.data.frame(scores), aes(sample = scores)) + 
  stat_qq(distribution = qt, dparams = as.list(fit_test$estimate)) + 
  stat_qq_line(distribution = qt, dparams = as.list(fit_test$estimate))


cau <- rcauchy(1000, location = fit_test$estimate['location'], scale = fit_test$estimate['scale'])

plot(density(cau), xlim = c(-5,5))
lines(density(scores))
