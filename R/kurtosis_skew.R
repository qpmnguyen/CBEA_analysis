library(tidyverse)
library(fitdistrplus)
library(glue)
library(goftest)
library(ggsci)
library(parameters)
library(mixtools)
library(goftest)
library(fitdistrplus)
library(patchwork)
library(CBEA)
library(Rfast)
#source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")
source("R/functions_data_ss.R")

set.seed(1020)

if(Sys.info()["sysname"] == "Darwin"){
  save_dir <- "../cilr_manuscript/figures"
} else {
  save_dir <- "../teailr_manuscript/manuscript/figures"
}


parameters <- create_parameters(list(
  rep = 1, 
  spar = c(0.2, 0.6, 0.8),
  n_inflate = c(100),
  s_rho = c(0.0, 0.2, 0.6)
))

# Generating data 
parameters$data <- map(parameters$param, .f = ~{
  zinb_simulation(n_samp = 1e3, spar = .x$spar, s_rho = .x$s_rho, b_rho = 0, eff_size = 1, 
                  vary_params = FALSE, n_tax = 2000, 
                  n_inflate = .x$n_inflate, n_sets = 1, samp_prop = 1, method = "normal")
})


# Compute scores 
parameters$scores <- map(parameters$data, .f = ~ {
    enrichment_analysis(physeq = .x$obj, set = .x$set, method = "cbea", 
                        metric = "auc", abund_values = "Counts", output = "raw", 
                        parametric = FALSE)
})



# Compute kurtosis and skewness with bootstrapped resampling
parameters$shape <- map(parameters$scores, .f = ~{
  kur <- vector(length = 100)
  skew <- vector(length = 100)
  for (i in 1:100){
    samp <- sample(dplyr::pull(.x,2), replace = T, size = nrow(.x))
    kur[i] <- Rfast::kurt(as.vector(samp))
    skew[i] <- Rfast::skew(as.vector(samp))
  }
  tibble(kurtosis = c(Rfast::kurt(dplyr::pull(.x,2)), kur), 
         skewness = c(Rfast::skew(dplyr::pull(.x,2)), skew),
         label = c("obs",rep("boot",100)))
})

# Generating plots  
shape_plot <- parameters %>% unnest(c(param, shape)) %>% 
  dplyr::rename("Set Size" = n_inflate, "Sparsity" = spar, "Corr:" = s_rho)
shpp_plot<- ggplot(shape_plot %>% arrange(label), 
                   aes(x = skewness, y = kurtosis, col = label, alpha = factor(label), size = factor(label))) + 
  geom_point() + 
  facet_grid(`Correlation` ~ `Sparsity`, labeller = label_both) + theme_bw() + 
  geom_hline(yintercept = 0, col = 'red') + 
  geom_vline(xintercept = 0, col = "red") + 
  scale_color_d3(name = "Label", labels = c("Bootstrapped", "Observed")) + 
  scale_alpha_manual(guide = "none", values = c(0.3, 1)) +
  scale_size_manual(guide = "none", values  = c(1,2.5)) + 
  labs(x = "Skewness", y = "Kurtosis") + 
  theme(legend.position = "bottom", legend.margin = margin(), axis.text.x = element_text(angle = 45))



# GoF with regards to the simulations  

generate_values <- function(distr, n=1e3){
  if (length(distr) == 2){
    values <- rnorm(n = n, mean = distr$mean, sd = distr$sd)
  } else {
    values <- rnormmix(n = n, lambda = distr$lambda, sigma = distr$sigma, mu = distr$mu)
  }
  return(values)
}


values <- parameters

values <- values %>% unnest(param) %>% filter(spar %in% c(0.2,0.6), s_rho %in% c(0.2, 0.6)) 
values$norm_param <- map(values$scores, ~estimate_distr(dplyr::pull(.x,2), distr = "norm", init = NULL))
values$mnorm_param <- map(values$scores, ~estimate_distr(dplyr::pull(.x,2), distr = "mnorm", init = NULL, maxrestarts=1000, epsilon = 1e-06, maxit= 1e5))
values$norm <- map(values$norm_param, generate_values) 
values$mnorm <- map(values$mnorm_param, generate_values)
values <- values %>% mutate(raw = map(scores, ~pull(.x,2))) %>% unnest(c(norm, mnorm, raw)) %>% 
  pivot_longer(c(norm, mnorm, raw), names_to = "distr", values_to = "val") %>% 
  dplyr::rename("Corr" = "s_rho", "Sparsity" = "spar") %>% 
  mutate(distr = case_when(
    distr == "mnorm" ~ "Mixture Normal", 
    distr == "norm" ~ "Normal", 
    TRUE ~ "Raw CBEA"
  ))

dens_plot <- ggplot(values, aes(x = val, fill = distr)) + 
  geom_density(alpha = 0.4, aes(col = distr)) + theme_bw() + 
  scale_fill_d3() + 
  scale_color_d3() + 
  facet_grid(Corr ~ Sparsity, scales = "free_y", labeller = label_both) + 
  labs(x = "Values", y = "Density", fill = "Distribution", col = "Distribution") +
  theme(legend.position = "bottom")

# KS statistic 
parameters <- parameters %>% unnest(param)

parameters$KS <- map(parameters$scores, ~{
  distr_norm <- estimate_distr(pull(.x,2), distr = "norm", init = NULL)
  distr_mnorm <- estimate_distr(pull(.x,2), distr = "mnorm", init = NULL, maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
  norm <- rnorm(1e5, mean = distr_norm$mean, sd = distr_norm$sd)
  mnorm <- mixtools::rnormmix(n = 1e5, lambda = distr_mnorm$lambda, sigma = distr_mnorm$sigma, mu = distr_mnorm$mu)
  ks_norm <- ks.test(pull(.x,2), norm)$statistic
  ks_mnorm <- ks.test(pull(.x,2), mnorm)$statistic
  tibble(Distribution = c("Normal", "Mixture Normal"), KS = c(ks_norm, ks_mnorm))
})

parameters <- parameters %>% unnest(KS)

gof_plot_sim <- ggplot(parameters %>% dplyr::rename("Corr" = "s_rho"), 
       aes(x = spar, y = KS, col = Distribution)) + 
  geom_point(size = 2) + geom_line() + 
  facet_wrap(~Corr, labeller = label_both) + 
  theme_bw() + scale_color_d3() + 
  labs(y = "Kolmogorov-Smirnov D Statistic", x = "Sparsity") + scale_fill_d3() + 
  theme(legend.position = "bottom", legend.margin = margin(), axis.text.x = element_text(angle = 45))



# Applying function to normal data sets  
# normal_data <- readRDS(file = "data/hmp_stool_16S.rds")
# X <- otu_table(normal_data) %>% t()
# A <- taxtab2A(tax_table(normal_data), "GENUS", full = FALSE)
# scores <- cilr(X = X, A = A, resample = F)
# 
# distr_norm <- estimate_distr(scores[,1], distr = "norm", init = NULL)
# distr_mnorm <- estimate_distr(scores[,1], distr = "mnorm", init = NULL)
# 
# norm <- rnorm(1e4, mean = distr_norm$mean, sd = distr_norm$sd)
# mnorm <- mixtools::rnormmix(n = 1e4, lambda = distr_mnorm$lambda, sigma = distr_mnorm$sigma, mu = distr_mnorm$mu)
# df <- bind_rows(tibble(Distribution = rep("Raw", nrow(scores)), value = scores[,1]), 
#                 tibble(Distribution = rep("Normal", 1e4), value = norm), 
#                 tibble(Distribution = rep("Mixture Normal", 1e4), value = mnorm))
# distr_plot <- ggplot(df, aes(x = value, fill = Distribution, col = Distribution)) + geom_density(alpha = 0.3, size = 1.25) + 
#   scale_fill_d3() + scale_color_d3() + labs(x = "Value", y = "Density") + theme_bw() 
# 
# test <- goftest::ad.test(x = scores[,1], null = "pnorm", estimated = T, distr_norm$mean, distr_norm$sd)
# 
# ks_results <- map_dfr(scores, ~{
#   dist_norm <- estimate_distr(.x, distr = "norm", init = NULL)
#   distr_mnorm <- estimate_distr(.x, distr = "mnorm", init = NULL)
#   norm <- rnorm(1e5, mean = distr_norm$mean, sd = distr_norm$sd)
#   mnorm <- mixtools::rnormmix(n = 1e5, lambda = distr_mnorm$lambda, sigma = distr_mnorm$sigma, mu = distr_mnorm$mu)
#   ks_norm <- ks.test(.x, norm)$statistic
#   ks_mnorm <- ks.test(.x, mnorm)$statistic
#   tibble(Distribution = c("Normal", "Mixture Normal"), KS = c(ks_norm, ks_mnorm))
# })
# 
# gof_plot <- ggplot(ks_results, aes(x = Distribution, y = KS, fill = Distribution)) + geom_boxplot() + theme_bw() + 
#   labs(y = "Kolmogrov D Statistic") + scale_fill_d3()
#combined_plt <- shpp_plot + (distr_plot/gof_plot) + plot_annotation(tag_levels = "A")

combined_plt <- (shpp_plot + gof_plot_sim)/dens_plot + plot_annotation(tag_levels = "A") 

combined_plt

ggsave(combined_plt, filename = "figures/kurtosis_skewness_gof.png", dpi = 300, 
       width = 10, height = 9)
ggsave(combined_plt, filename = "figures/kurtosis_skewness_gof.eps", 
       dpi = 300, width = 10, height = 9, device = cairo_ps)

file.copy(from = Sys.glob("figures/*.png"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.eps"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
