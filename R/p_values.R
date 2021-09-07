library(tidyverse)
library(ggsci)
library(glue)
source("R/simulations.R")
source("R/utils.R")
source("R/cilr.R")

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}


if (!file.exists("additional_files/pval_distr.rds")){
    grid <- create_parameters(list(
        spar = c(0, 0.6),
        s_rho = c(0, 0.5)
    ))
    data <- zinb_simulation(n_samp = 1e4, spar = 0.2, s_rho = 0, eff_size = 1, 
                            n_tax = 1000, n_sets = 1, n_inflate = 100, vary_params = TRUE)
    grid <- grid %>% unnest(param) %>% rowwise() %>% mutate(data = list(
        zinb_simulation(n_samp = 1e4, spar = spar, s_rho = s_rho, eff_size = 1, 
                        n_tax = 500, n_sets = 1, n_inflate = 50, vary_params = TRUE)
    ))
    grid_eval <- crossing(grid, distr = c("norm", "mnorm")) %>% 
        mutate(adj = if_else(s_rho > 0, TRUE, FALSE))

    grid_eval <- grid_eval %>% rowwise() %>% mutate(model = list(
        cilr(X = data$X, A = data$A, resample = TRUE, output = "pval", 
             distr = distr, adj = adj, maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
    ))
} else {
    grid_eval <- readRDS(file = "additional_files/pval_distr.rds")
}



grid <- grid_eval %>% unnest(model) %>% rename(pval = Set1ss) %>% 
    dplyr::select(-data)

grid <- grid %>% rename(Sparsity = spar, Correlation = s_rho) %>% 
    mutate(distr = recode(distr, "mnorm" = "Mixture Normal", "norm" = "Normal"))

pval_plot <- ggplot(grid, aes(col = distr, sample = pval)) + 
    facet_grid(Sparsity ~ Correlation, labeller = label_both) +
    stat_qq(distribution = qunif) + 
    stat_qq_line(distribution = qunif, linetype = "dashed") + 
    theme_bw() + scale_color_d3() + 
    guides(color = guide_legend(override.aes = list(alpha = 1))) + 
    labs(x = "Theoretical", y = "Sample", col = "Distribution")
ggsave(pval_plot, filename = "figures/pval_distr.png", dpi = 300, width = 8, height = 6)
ggsave(pval_plot, filename = "figures/pval_distr.eps", 
       dpi = 300, width = 8, height = 6)

file.copy(from = Sys.glob("figures/*.png"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.eps"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
