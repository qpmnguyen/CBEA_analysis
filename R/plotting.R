library(tidyverse)
library(patchwork)
library(glue)
library(ggsci)

# FDR
plotting <- readRDS(file = "objects/fdr_sim/fdr_eval.rds")
plotting <- plotting %>% unnest(param) %>% mutate(eval = map(eval, as_tibble)) %>% unnest(eval)

plotting <- plotting %>% rename("Set Size" = "n_inflate", 
"Correlation" = "s_rho")
fdr_plot <- ggplot(plotting, aes(y = mean, x = spar, col = distr)) + geom_point() +
    geom_line(aes(linetype = adj)) + 
    facet_grid(`Correlation` ~ `Set Size`, labeller = "label_both", 
    scales = "free_y") + 
    geom_hline(yintercept = 0.05, col = "red") + 
    scale_color_npg(labels = c("Mixture normal", "Normal", "Wilcoxon Rank Sum Test")) + theme_bw() + 
    labs(x = "Sparsity", y = "Type I error (N = 2,000)", col = "Method", linetype = "Correlation adjusted") +
    geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.01)

ggsave(fdr_plot, file = "docs/manuscript/figures/ss_fdr_plot.png", dpi = 300, width = 8, height = 5)

# Plotting PWR
plotting <- readRDS(file = "objects/pwr_sim/pwr_eval.rds")
plotting <- plotting %>% unnest(param) %>% mutate(eval = map(eval, as_tibble)) %>% unnest(eval)

plotting <- plotting %>% rename("Set Size" = "n_inflate", 
"Correlation" = "s_rho", `Effect Size` = eff_size)
pwr_plot <- ggplot(plotting, aes(y = mean, x = spar, col = distr)) + geom_point() +
    geom_line(aes(linetype = adj)) + 
    facet_grid(`Correlation` ~ `Effect Size`, labeller = "label_both", 
    scales = "free_y") + 
    geom_hline(yintercept = 0.8, col = "red") + 
    scale_color_npg(labels = c("Mixture normal", "Normal", "Wilcoxon Rank Sum Test")) + theme_bw() + 
    labs(x = "Sparsity", y = "Power (N = 2,000)", col = "Method", linetype = "Correlation adjusted") +
    geom_errorbar(aes(ymax = upper, ymin = lower), width = 0.01)

ggsave(pwr_plot, file = "docs/manuscript/figures/ss_pwr_plot.png", dpi = 300, width = 8, height = 5)
