library(tidyverse)
library(fitdistrplus)
library(glue)
library(goftest)
library(ggsci)
library(parameters)
library(mixtools)
library(goftest)

source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

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


map(parameters$scores, .f = ~{
  kur <- vector(length = 100)
  skew <- vector(length = 100)
  for (i in seq_len(100)){
    samp <- sample(.x, size = length(.x))
  }
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




#ggsave(shpp_plot, filename = "docs/manuscript/figures/kurtosis_skewness_sim.png", dpi = 300, width = 8, height = 5)
