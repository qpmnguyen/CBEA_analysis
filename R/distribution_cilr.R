# Distribution of the cilr test statistic under the null 
# Quang Nguyen
# Last updated 9/24

library(tidyverse)
library(ggsci)
library(fitdistrplus)
library(glue)
library(patchwork)
library(propagate)
library(gnorm)
library(SuppDists)
library(VGAM)
library(goftest)
library(furrr)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")



# Function to get fit and estimate fit from data, calculate ks test, aic and bic. If using fitdistplus, obtain 
# aic and bic directly. Else, use custon formula 
get_fit <- function(scores, thresh = 10){
  dist_list <- c("norm", "t", "gnorm", "laplace", "cauchy", "logistic", "JSU")
  result <- tibble(dist = dist_list, AIC = rep(0,length(dist_list)), BIC = rep(0, length(dist_list)),
                   AD = rep(0, length(dist_list)))
  for (i in 1:nrow(result)){
    if (result[i,]$dist %in% c("norm", "t", "gnorm", "laplace", "cauchy", "logistic", "JSU")){
      if(result[i,]$dist == "norm"){
        fit <- fitdist(scores, distr = "norm")
      } else if (result[i,]$dist == "t"){
        fit <- try(fitdist(scores, distr = "t", start = list(df = 1)), silent = T)
      } else if (result[i,]$dist == "gnorm"){
        fit <- fitdist(scores, distr = "gnorm", start = list(mu = 0, alpha = 3, beta = 2))
      } else if (result[i,]$dist == "laplace"){
        fit <- fitdist(scores, distr = "laplace", start = list(location = 0, scale = 1))
      } else if (result[i,]$dist == "cauchy"){
        fit <- fitdist(scores, distr = "cauchy", start = list(location = 0, scale = 1))
      } else if (result[i,]$dist == "logistic"){
        fit <- fitdist(scores, distr = "logis", start = list(location = 0, scale = 1))
      } else if (result[i,]$dist == "JSU"){
        fit <- fitdist(scores, distr = "JSU", start = list(mu = 0, sigma = 3, nu = 2, tau = 2), 
                       control = list(maxit = 1000))
      }
      try(stat <- gofstat(fit), silent = T)
      result[i,]$AIC <- stat$aic
      result[i,]$BIC <- stat$bic
      result[i,]$AD <- stat$ad
    } else {
      if (result[i,]$dist == "Johnson"){
        fit <- JohnsonFit(scores)
        loglikelihood <- sum(dJohnson(x = scores, parms = fit, log = T))
        dist_name <- "pJohnson"
      }
      result[i,]$AIC <- 2*4 - 2*loglikelihood
      result[i,]$BIC <- 4*log(length(scores)) - 2*loglikelihood 
      try(result[i,]$AD <- ad.test(scores, null = dist_name, fit, estimated = T)$statistic, silent = F)
    }
  }
  return(result)
}

parameters <- create_parameters(list(
  rep = seq(1,10),
  spar = c(0.2, 0.4, 0.6, 0.8),
  n_inflate = c(50,100,150),
  s_rho = c(0.1, 0.3, 0.5)
))

# Generating data 
plan(multisession, workers = 3)
opt <- furrr_options(seed = TRUE)
parameters$data <- future_map(parameters$param, .f = ~{
  zinb_simulation(n_samp = 1e3, spar = .x$spar, s_rho = .x$s_rho, b_rho = 0, eff_size = 1, 
                  vary_params = FALSE, n_tax = 1000, 
                  n_inflate = .x$n_inflate, n_sets = 1, samp_prop = 1, method = "normal")
},.options = opt, .progress = TRUE)
plan(sequential)

# cILR scores  
parameters$scores <- map(parameters$data, .f = ~ {
  simple_cilr(X = .x$X, A = .x$A, pcount = 1, transform = NULL, preprocess = T, method = "raw")
})


# Evaluate fit using custom function
parameters$fit_eval <- map(parameters$scores, .f = ~{
  get_fit(as.vector(.x))
})

fit_plot <- parameters %>% unnest(c(param, fit_eval)) %>% group_by(spar, n_inflate, s_rho, dist) %>% 
  mutate(meanAIC = mean(AIC), meanBIC = mean(BIC), meanAD = mean(AD), 
         upperAIC = mean(AIC) + sd(AIC), upperBIC = mean(BIC) + sd(BIC), upperAD = mean(AD) + sd(AD),
         lowerAIC = mean(AIC) - sd(AIC), lowerBIC = mean(BIC) - sd(BIC), lowerAD = mean(AD) - sd(AD)) %>%
  rename("Set Size" = n_inflate, "Inter-taxa correlation" = s_rho)

aic_plt <- ggplot(fit_plot, aes(x = spar, y = meanAIC, col = dist)) + geom_point() + 
  geom_errorbar(aes(ymax = upperAIC, ymin = lowerAIC), width = 0.05) + 
  geom_line() + facet_grid(`Set Size` ~ `Inter-taxa correlation`) + scale_color_lancet() + theme_bw()
aic_plt
bic_plt <- ggplot(fit_plot, aes(x = spar, y = meanBIC, col = dist)) + geom_point() + 
  geom_errorbar(aes(ymax = upperBIC, ymin = lowerBIC), width = 0.05) + 
  geom_line() + facet_grid(`Inter-taxa correlation` ~ `Set Size`, labeller = label_both) + 
  scale_color_npg(name = "Distributions", labels = c("Cauchy", "Generalized Normal", "Johnson SU", 
                                                     "Laplace", "Logistic", "Normal", "Student's t")) + 
  theme_bw() +
  labs(y = "Bayesian Information Criterion (10 data sets, N = 1000)", x = "Sparsity")


ad_plt <- ggplot(fit_plot, aes(x = spar, y = meanAD, col = dist)) + geom_point() + 
  geom_errorbar(aes(ymax = upperAD, ymin = lowerAD), width = 0.05) + 
  geom_line() + facet_grid(`Inter-taxa correlation` ~ `Set Size`, labeller = label_both) + 
  scale_color_npg(name = "Distributions", labels = c("Cauchy", "Generalized Normal", "Johnson SU", 
                                                     "Laplace", "Logistic", "Normal", "Student's t")) + theme_bw() + 
  labs(y = "Anderson-Darling GoF (10 data sets, N = 1000)", x = "Sparsity")


ftplot <- bic_plt + ad_plt + plot_layout(guides = "collect", tag_level = "new") 
ggsave(ftplot, filename = "docs/manuscript/figures/fit_evaluation.png", dpi = 300, width = 15, height = 8)


# TODO: Create bootstrapped estimates for kurtosis and skewness
parameters$kurtosis <- map_dbl(parameters$scores, .f = ~{
  boot <- vector(length = 500)
  for (i in 1:length(boot)){
    sample(.x, replace = T, size = length(.x))
  }
  kurtosis(as.vector(.x))
  return(boot)
})

parameters$skewness <- map_dbl(parameters$scores, .f = ~{
  skewness(as.vector(.x))
})
