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
  rep = seq(1,10),
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
    zinb_simulation(n_samp = 300, b_spar = .x$b_spar, b_rho = .x$b_rho, 
                    eff_size = .x$eff_size, n_inflate = 50, rho_ratio = 1)
  }, .options = opt)
})
toc()
plan(sequential)


# getting scores 
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
  pwr_sim$scores_gsva_pois <- future_map(pwr_sim$sim, .f = ~{
    p()
    generate_alt_scores(X = .x$X, A = .x$A, method = "gsva", preprocess = T, transform=NULL, pcount=1)
  })
})

with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$scores_gsva_gauss <- future_map(pwr_sim$sim, .f = ~{
    p()
    generate_alt_scores(X = .x$X, A = .x$A, method = "gsva", preprocess = T, transform="clr", pcount=1)
  })
})

with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$scores_ssgsea_pois <- future_map(pwr_sim$sim, .f = ~{
    p()
    generate_alt_scores(X = .x$X, A = .x$A, method = "ssgsea", preprocess = T, transform=NULL, pcount=1)
  })
})

with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$scores_ssgsea_gauss <- future_map(pwr_sim$sim, .f = ~{
    p()
    generate_alt_scores(X = .x$X, A = .x$A, method = "ssgsea", preprocess = T, transform="clr", pcount=1)
  })
})

plan(sequential)

# get auc scores 
scores <- pwr_sim %>% dplyr::select(c(sim, starts_with("scores")))
auc <- vector(mode = "list", length = length(3:ncol(scores))-1)
for (i in 3:ncol(scores)){
  new_name <- glue("auc_{name}", name = colnames(scores)[i])
  values <- map_dbl(1:nrow(scores), .f = function(.x){
    stat <- calculate_statistic(eval = "auc", pred = scores[.x,i][[1]], 
                                true = scores$sim[[.x]]$label)
    return(stat)
  })
  auc[[new_name]] <- values
}
auc <- do.call(cbind, auc)
auc <- as_tibble(auc)
pwr_sim <- cbind(pwr_sim, auc)

plotting_data <- pwr_sim %>% dplyr::select(c(id,param, starts_with("auc"))) %>% pivot_longer(starts_with("auc")) %>% 
  unnest(param) %>% group_by(b_spar, b_rho, eff_size, name) %>% summarise(min = quantile(value, 0.05, na.rm = T), 
                                                                         median = median(value, na.rm = T),
                                                                         max = max(value, 0.95, na.rm = T), 
                                                                         hpdi_lower = rethinking::HPDI(value, prob = 0.95)[1],
                                                                         hpdi_upper = rethinking::HPDI(value, prob = 0.95)[2])
plotting_data <- plotting_data %>% filter(!name == "auc_scores_ssgsea_gauss")
plotting_data <- plotting_data %>% dplyr::rename("Effect Size" = "eff_size", "Inter-taxa correlation" = "b_rho")

auc_plot <- ggplot(plotting_data, aes(x = b_spar, y = median, col = name)) + geom_hline(yintercept = 0.8, col = "red") + 
  geom_point(size = 2) + facet_grid(`Inter-taxa correlation` ~ `Effect Size`, labeller = label_both) + theme_bw() +
  geom_errorbar(aes(ymin = hpdi_lower, ymax = hpdi_upper),width = 0.04) + 
  geom_line() + scale_color_nejm(labels = c("cILR", "GSVA (Gaussian Kernel)", "GSVA (Poisson Kernel)", "ssGSEA")) + 
  labs(col = "Models", x = "Sparsity", y = "Mediaan AUC (10 data sets) with HD intervals")
auc_plot

ggsave(auc_plot, filename = "docs/manuscript/figures/auc_plots.png", dpi = 300, width = 10, height = 10)


