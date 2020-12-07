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

# Plotting AUC 
parameters <- readRDS(file = "objects/auc_sim/parameters.rds")
results <- readRDS(file = "objects/auc_sim/auc_eval.rds")
results <- results %>% unnest(eval)
results <- results %>% filter(!distr %in% c("gsva", "ssgsea"))
results <- results %>% group_by(spar, s_rho, eff_size, method, distr, output, adj) %>% summarise(auc = mean(eval), upper = mean(eval) + sd(eval), lower = mean(eval) - sd(eval))

ggplot(results, aes(x = spar, y = auc, color = distr, shape = output, linetype = adj)) + geom_line() + 
    facet_grid(s_rho ~ eff_size, labeller = label_both, scales = "free_y") + 
    geom_pointrange(aes(ymax = upper, ymin = lower)) + 
    scale_color_npg() + theme_bw() + 
    labs(y = "Mean AUC (10 data sets)", x = "Sparsity", shape = "Output type", linetype = "Correlation Adjustment", 
         color = "Model type")
#labels = c("GSVA", "Mixture Normal", "Raw", "Normal", "ssGSEA")