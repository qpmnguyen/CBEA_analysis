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
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

set.seed(1020)

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
  cilr(X = .x$X, A = .x$A, resample = F, pcount = 1, transform = "prop", preprocess = T)
})



# Compute kurtosis and skewness with bootstrapped resampling
parameters$shape <- map(parameters$scores, .f = ~{
  kur <- vector(length = 100)
  skew <- vector(length = 100)
  for (i in 1:100){
    samp <- sample(.x[,1], replace = T, size = length(.x[,1]))
    kur[i] <- kurtosis(as.vector(samp))$Kurtosis
    skew[i] <- skewness(as.vector(samp))$Skewness
  }
  tibble(kurtosis = c(kurtosis(as.vector(.x))$Kurtosis, kur), 
         skewness = c(skewness(as.vector(.x))$Skewness, skew),
         label = c("obs",rep("boot",100)))
})

# Generating plots  
shape_plot <- parameters %>% unnest(c(param, shape)) %>% 
  rename("Set Size" = n_inflate, "Sparsity" = spar, "Correlation" = s_rho)
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
  theme(legend.position = "bottom", legend.margin = margin())



# GoF with regards to the simulations  

parameters <- parameters %>% unnest(param)

parameters$KS <- map(parameters$scores, ~{
  distr_norm <- estimate_distr(.x[,1], distr = "norm", init = NULL)
  distr_mnorm <- estimate_distr(.x[,1], distr = "mnorm", init = NULL, maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
  norm <- rnorm(1e5, mean = distr_norm$mean, sd = distr_norm$sd)
  mnorm <- mixtools::rnormmix(n = 1e5, lambda = distr_mnorm$lambda, sigma = distr_mnorm$sigma, mu = distr_mnorm$mu)
  ks_norm <- ks.test(.x, norm)$statistic
  ks_mnorm <- ks.test(.x, mnorm)$statistic
  tibble(Distribution = c("Normal", "Mixture Normal"), KS = c(ks_norm, ks_mnorm))
})

parameters <- parameters %>% unnest(KS)

gof_plot_sim <- ggplot(parameters %>% rename("Correlation" = "s_rho"), 
       aes(x = spar, y = KS, col = Distribution)) + 
  geom_point(size = 2) + geom_line() + 
  facet_wrap(~Correlation, labeller = label_both) + 
  theme_bw() + scale_color_d3() + 
  labs(y = "Kolmogrov D Statistic", x = "Sparsity") + scale_fill_d3() + 
  theme(legend.position = "bottom", legend.margin = margin())



# Applying function to normal data sets  
normal_data <- readRDS(file = "data/hmp_stool_16S.rds")
X <- otu_table(normal_data) %>% t()
A <- taxtab2A(tax_table(normal_data), "GENUS", full = FALSE)
scores <- cilr(X = X, A = A, resample = F)

distr_norm <- estimate_distr(scores[,1], distr = "norm", init = NULL)
distr_mnorm <- estimate_distr(scores[,1], distr = "mnorm", init = NULL)

norm <- rnorm(1e4, mean = distr_norm$mean, sd = distr_norm$sd)
mnorm <- mixtools::rnormmix(n = 1e4, lambda = distr_mnorm$lambda, sigma = distr_mnorm$sigma, mu = distr_mnorm$mu)
df <- bind_rows(tibble(Distribution = rep("Raw", nrow(scores)), value = scores[,1]), 
                tibble(Distribution = rep("Normal", 1e4), value = norm), 
                tibble(Distribution = rep("Mixture Normal", 1e4), value = mnorm))
distr_plot <- ggplot(df, aes(x = value, fill = Distribution, col = Distribution)) + geom_density(alpha = 0.3, size = 1.25) + 
  scale_fill_d3() + scale_color_d3() + labs(x = "Value", y = "Density") + theme_bw() 

test <- goftest::ad.test(x = scores[,1], null = "pnorm", estimated = T, distr_norm$mean, distr_norm$sd)

ks_results <- map_dfr(scores, ~{
  dist_norm <- estimate_distr(.x, distr = "norm", init = NULL)
  distr_mnorm <- estimate_distr(.x, distr = "mnorm", init = NULL)
  norm <- rnorm(1e5, mean = distr_norm$mean, sd = distr_norm$sd)
  mnorm <- mixtools::rnormmix(n = 1e5, lambda = distr_mnorm$lambda, sigma = distr_mnorm$sigma, mu = distr_mnorm$mu)
  ks_norm <- ks.test(.x, norm)$statistic
  ks_mnorm <- ks.test(.x, mnorm)$statistic
  tibble(Distribution = c("Normal", "Mixture Normal"), KS = c(ks_norm, ks_mnorm))
})

gof_plot <- ggplot(ks_results, aes(x = Distribution, y = KS, fill = Distribution)) + geom_boxplot() + theme_bw() + 
  labs(y = "Kolmogrov D Statistic") + scale_fill_d3()


#combined_plt <- shpp_plot + (distr_plot/gof_plot) + plot_annotation(tag_levels = "A")
combined_plt <- gof_plot_sim + shpp_plot + plot_annotation(tag_levels = "A") 


ggsave(combined_plt, filename = "figures/kurtosis_skewness_gof.png", dpi = 300, width = 8, height = 5)
file.copy("figures/kurtosis_skewness_gof.png", 
          "../teailr_manuscript/manuscript/figures/kurtosis_skewness_gof.png", overwrite = TRUE)
