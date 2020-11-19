# Distribution of the cilr test statistic under the null 
# Quang Nguyen
# Last updated 9/24

library(tidyverse)
library(fitdistrplus)
library(glue)
library(gnorm)
library(gamlss.dist)
library(VGAM)
library(goftest)
library(furrr)
library(parameters)
library(mixtools)
library(qs)
library(sn)
library(goftest)

source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")

data <- qread("objects/fdr_sim/simulation_21.qs")
sim <- qread("objects/fdr_sim/parameters.qs")
sim
cilr_adj <- cilr(X = data$X, A = data$A, resample = T, distr = "mnorm", adj = T,
                 output = "pval",
                 maxit = 1e6, maxrestarts = 1e3, epsilon = 1e-6)
cilr_unadj <- cilr(X = data$X, A = data$A, resample = T, distr = "mnorm", adj = F,
                 output = "pval",
                 maxit = 1e6, maxrestarts = 1e3, epsilon = 1e-6)

hist(cilr_adj, main = "Histogram of p-values", xlab = "Adjusted cILR w/ mixture normal distribution", col = "steelblue")
hist(cilr_unadj, main = "Histogram of p-values", xlab = "Unadjusted cILR w/ mixture normal distribution", col = "steelblue")

qqplot(x = qunif(ppoints(2e4), 0,1), cilr_adj, pch = 18)
qqline(cilr_adj, distribution = function(p) qunif(p, 0, 1), col = "salmon", lwd = 2)





# Fitting distributions to data 
# Function to get fit and estimate fit from data, calculate ks test, aic and bic. If using fitdistplus, obtain 
# aic and bic directly. Else, use custon formula 
get_fit <- function(scores, thresh = 10){
  #dist_list <- c("norm", "t", "gnorm", "laplace", "cauchy", "logistic", "mnorm", "st", "JSU")
  dist_list <- c("norm", "mnorm", "t")
  result <- tibble(dist = dist_list, AIC = rep(0,length(dist_list)), BIC = rep(0, length(dist_list)),
                   AD = rep(0, length(dist_list)))
  for (i in 1:nrow(result)){
    if (result[i,]$dist %in% c("norm", "t", "gnorm", "laplace", "cauchy", "logistic", "JSU", "st")){
      if(result[i,]$dist == "norm"){
        fit <- fitdist(scores, distr = "norm")
      } else if (result[i,]$dist == "t"){
        fit <- try(fitdist(scores, distr = "t", start = list(df = 1)), silent = T)
      } else if (result[i,]$dist == "gnorm"){
        fit <- fitdist(scores, distr = "gnorm", start = list(mu = 0, alpha = 1, beta = 2), control = list(maxit = 10000))
      } else if (result[i,]$dist == "laplace"){
        fit <- fitdist(scores, distr = "laplace", start = list(location = 0, scale = 1))
      } else if (result[i,]$dist == "cauchy"){
        fit <- fitdist(scores, distr = "cauchy", start = list(location = 0, scale = 1))
      } else if (result[i,]$dist == "logistic"){
        fit <- fitdist(scores, distr = "logis", start = list(location = 0, scale = 1))
      } else if (result[i,]$dist == "JSU"){
        fit <- fitdist(scores, distr = "JSU", start = list(mu = 0, sigma = 3, nu = 2, tau = 2), 
                       control = list(maxit = 10000))
      } else if (result[i,]$dist == "st"){
        fit <- fitdist(scores, distr = "st", start = list(xi = 0, omega = 1, alpha = 1, nu = 1), 
                       control = list(maxit = 10000))
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
      } else if (result[i,]$dist == "mnorm"){
        fit <- normalmixEM(x = scores)
        loglikelihood <- fit$loglik
        dist_name <- "pmnorm"
        fit <- list(mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda)
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
  spar = c(0.6),
  n_inflate = c(100),
  s_rho = c(0.0, 0.2, 0.4, 0.6)
))

# Generating data 
plan(multisession, workers = 3)
opt <- furrr_options(seed = TRUE)
parameters$data <- future_map(parameters$param, .f = ~{
  zinb_simulation(n_samp = 1e3, spar = .x$spar, s_rho = .x$s_rho, b_rho = 0, eff_size = 1, 
                  vary_params = FALSE, n_tax = 2000, 
                  n_inflate = .x$n_inflate, n_sets = 1, samp_prop = 1, method = "normal")
},.options = opt, .progress = TRUE)
plan(sequential)

# cILR scores  
parameters$scores <- map(parameters$data, .f = ~ {
  cilr(X = .x$X, A = .x$A, resample = F, pcount = 1, transform = "prop", preprocess = T)
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

# Remove t distribution

fit_plot <- fit_plot %>% filter(dist != "t")

bic_plt <- ggplot(fit_plot, aes(x = spar, y = meanBIC, col = dist)) + geom_point() + 
  geom_errorbar(aes(ymax = upperBIC, ymin = lowerBIC), width = 0.05) + 
  geom_line() + facet_grid(`Inter-taxa correlation` ~ `Set Size`, labeller = label_both) + 
  scale_color_npg(name = "Distributions", labels = c("Cauchy", "Generalized Normal", "Johnson SU", 
                                                     "Laplace", "Logistic", "Mixture Normal", "Normal", "Student's t")) + 
  theme_bw() +
  labs(y = "Bayesian Information Criterion (10 data sets, N = 1000)", x = "Sparsity")


ad_plt <- ggplot(fit_plot, aes(x = spar, y = meanAD, col = dist)) + geom_point() + 
  geom_errorbar(aes(ymax = upperAD, ymin = lowerAD), width = 0.05) + 
  geom_line() + facet_grid(`Inter-taxa correlation` ~ `Set Size`, labeller = label_both) + 
  scale_color_npg(name = "Distributions", labels = c("Cauchy", "Generalized Normal", "Johnson SU", 
                                                       "Laplace", "Logistic", "Mixture Normal", "Normal", "Student's t")) + theme_bw() + 
  labs(y = "Anderson-Darling GoF (10 data sets, N = 1000)", x = "Sparsity")


ftplot <- bic_plt + ad_plt + plot_layout(guides = "collect", tag_level = "new") 
ggsave(ftplot, filename = "docs/manuscript/figures/fit_evaluation_without_t.png", dpi = 300, width = 15, height = 8)


# Compute kurtosis and skewness with bootstrapped resampling
parameters$shape <- map(parameters$scores, .f = ~{
  kur <- vector(length = 100)
  skew <- vector(length = 100)
  for (i in 1:100){
    samp <- sample(.x, replace = T, size = length(.x))
    kur[i] <- kurtosis(as.vector(samp))$Kurtosis
    skew[i] <- skewness(as.vector(samp))$Skewness
  }
  tibble(kurtosis = c(kurtosis(as.vector(.x))$Kurtosis, kur), 
         skewness = c(skewness(as.vector(.x))$Skewness, skew),
         label = c("obs",rep("boot",100)))
})


shape_plot <- parameters %>% unnest(c(param, shape)) %>% rename("Set Size" = n_inflate, "Sparsity" = spar)



low_corr <- ggplot(shape_plot %>% filter(s_rho == 0.1) %>% arrange(label), 
                   aes(x = skewness, y = kurtosis, col = label, alpha = factor(label), size = factor(label))) + 
            geom_point() +
            facet_grid(`Set Size` ~ `Sparsity`, labeller = label_both) + 
            scale_color_npg(name = "Label", labels = c("Bootstrapped", "Observed")) + 
            scale_alpha_manual(guide = "none", values = c(0.3, 1)) + 
            scale_size_manual(guide = "none", values  = c(1,2.5)) +
            labs(x = "Skewness", y = "Kurtosis", subtitle = "Low Inter-taxa Correlation") + theme_bw()

med_corr <- ggplot(shape_plot %>% filter(s_rho == 0.3) %>% arrange(label), 
                   aes(x = skewness, y = kurtosis, col = label, alpha = factor(label), size = factor(label))) + 
  geom_point() +
  facet_grid(`Set Size` ~ `Sparsity`, labeller = label_both) + 
  scale_color_npg(name = "Label", labels = c("Bootstrapped", "Observed")) + 
  scale_alpha_manual(guide = "none", values = c(0.3, 1)) + 
  scale_size_manual(guide = "none", values  = c(1,2.5)) +
  labs(x = "Skewness", y = "Kurtosis", subtitle = "Medium Inter-taxa Correlation") + theme_bw()

high_corr <- ggplot(shape_plot %>% filter(s_rho == 0.5) %>% arrange(label), 
                    aes(x = skewness, y = kurtosis, col = label, alpha = factor(label), size = factor(label))) + 
  geom_point() +
  facet_grid(`Set Size` ~ `Sparsity`, labeller = label_both) + 
  scale_color_npg(name = "Label", labels = c("Bootstrapped", "Observed")) + 
  scale_alpha_manual(guide = "none", values = c(0.3, 1)) + 
  scale_size_manual(guide = "none", values  = c(1,2.5)) +
  labs(x = "Skewness", y = "Kurtosis", subtitle = "High Inter-taxa Correlation") + theme_bw()

shpplot <- low_corr + med_corr + high_corr +  plot_layout(guides = "collect", tag_level = "new") 
ggsave(shpplot, filename = "docs/manuscript/figures/kurtosis_skewness_sim.png", dpi = 300, width = 20, height = 10)



# Real data  
data <- qread(file = "data/hmp_stool_16S.qs")
otu_tab <- otu_table(data)

A <- taxtab2A(tax_table(data), agg_level = "GENUS")
X <- unclass(t(otu_table(data)))

real <- simple_cilr(X = X, A = A, preprocess = T, pcount = 1, resample = F, method = "raw")

otu_tab <- otu_tab[sample(1:nrow(otu_tab), replace = F),]
otu_table(data) <- otu_table(otu_tab, taxa_are_rows = T)
A <- taxtab2A(tax_table(data), agg_level = "GENUS")
X <- unclass(t(otu_table(data)))

perm <- simple_cilr(X = X, A = A, preprocess = T, pcount = 1, resample = F, method = "raw")


fitted <- apply(cilr_scores, 2, function(.x){
  get_fit(scores = as.vector(.x))
})





real_data_fit <- as_tibble_col(fitted) %>% mutate(label = names(fitted)) %>% unnest(value)

bic_plot <- ggplot(real_data_fit, aes(y = BIC, x = dist)) + geom_boxplot(fill = "steelblue", alpha = 0.8) + 
  scale_x_discrete(labels = c("cauchy" = "Cauchy", "gnorm" = "Generalized \n Normal", "JSU" = "Johnson \n SU", 
                              "laplace" = "Laplace", "logistic" = "Logistic", "mnorm" = "Mixture \n Normal", 
                              "norm" = "Normal")) +
  labs(x = "Distribution", y = "Bayesian Information Criterion (40 sets)") + theme_bw() 
ad_plot <- ggplot(real_data_fit, aes(y = AD, x = dist)) + geom_boxplot(fill = "steelblue", alpha = 0.8) + 
  scale_x_discrete(labels = c("cauchy" = "Cauchy", "gnorm" = "Generalized \n Normal", "JSU" = "Johnson \n SU", 
                              "laplace" = "Laplace", "logistic" = "Logistic", "mnorm" = "Mixture \n Normal", 
                              "norm" = "Normal")) +
  labs(x = "Distribution", y = "Anderson Darling GoF (40 sets)") + theme_bw() 

realdatplot <- bic_plot + ad_plot
ggsave(realdatplot, filename = "docs/manuscript/figures/real_data_fit.png", dpi = 300, width = 10, height = 10)

# Generating simulated values based on fitted distributions
generate_based_fit <- function(scores, n_sim = 1e4, distr = "t"){
  scores <- as.vector(scores)
  if (distr == "t"){
    fit <- fitdist(scores, distr = "t", method = "mle", start = list(df = 2), control = list(maxit = 1e4))
    data <- rt(n_sim, df = fit$estimate['df'])
    parms <- list(df = fit$estimate['df'])
  } else if (distr == "norm"){
    fit <- fitdist(scores, distr = "norm", method = "mle")
    data <- rnorm(n_sim, mean = fit$estimate['mean'], sd = fit$estimate['sd'])
    parms <- list(mean = fit$estimate['mean'], sd = fit$estimate['sd'])
  } else if (distr == "mnorm"){
    fit <-  normalmixEM(scores)
    data <- rnormmix(n_sim, lambda = fit$lambda, mu = fit$mu, sigma = fit$sigma)
    parms <- list(lambda = fit$lambda, mu = fit$mu, sigma = fit$sigma)
  }
  return(data)
}

# Plot for density
colors <- pal_npg('nrc')(3)
plot_fit <- qplot(x = cilr_scores[,1], geom = "blank", xlim = c(min(cilr_scores[,1]), max(cilr_scores[,1]))) + 
  geom_histogram(aes(y = stat(density)), fill = "steelblue", alpha = 0.8) + 
  geom_density(data = data.frame(fit = generate_based_fit(cilr_scores[,1], distr = 't')), 
               aes(x = fit, col = colors[1]), size = 2) +
  geom_density(data = data.frame(fit = generate_based_fit(cilr_scores[,1], distr = 'norm')), 
               aes(x = fit, col = colors[2]), size = 2) +
  geom_density(data = data.frame(fit = generate_based_fit(cilr_scores[,1], distr = 'mnorm')), 
               aes(x = fit, col = c('test' = colors[3])), size = 2) +
  scale_color_manual(name = "Distributions", labels = c("Mixed Normal", "Normal", "t"), values = colors) + 
  theme_bw() + labs(x = "cILR scores from null HMP data (1 set)", y = "Density")

ggsave(plot_fit, filename = "docs/manuscript/figures/distr_density_plot_real.png", dpi = 300, width = 5, height = 5)

