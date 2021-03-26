# Make sure to set working directory to be the main working directory as teailr
library(tidyverse)
library(phyloseq)
library(ggsci)
library(patchwork)

fdr <- readRDS(file = "analyses/simulations_single_sample_fdr/output/sim_ss_fdr.rds")
fdr_grid <- readRDS(file = "analyses/simulations_single_sample_fdr/output/simulation_grid_fdr.rds") 


fdr_grid <- fdr_grid %>% dplyr::select(spar, s_rho, n_inflate, id)
fdr_df <- inner_join(fdr, fdr_grid, by = "id") %>% 
    mutate(adj = replace_na(adj, "Not Applicable"), distr = replace_na(distr, "Not Applicable")) %>% 
    mutate(distr = recode(distr, mnorm = "Mixture Normal", norm = "Normal"), 
           model = recode(model, cilr = "cILR", wilcox = "Wilcoxon Rank Sum"), 
           adj = recode(adj, "FALSE" = "No", "TRUE" = "Yes")) %>% 
    rename(c("Set Size" = "n_inflate", "Correlation" = "s_rho"))

fdr_plt <- ggplot(fdr_df, aes(x = spar, y = est, col = model, linetype = adj, shape = distr)) + 
    geom_linerange(aes(ymax = upper, ymin = lower)) + 
    geom_point() + 
    geom_line() + 
    geom_hline(yintercept = 0.05, col = "red") + 
    scale_color_d3() +
    theme_bw() + 
    labs(x = "Sparsity", y = "Type I error", col = "Model type", shape = "Distribution", 
         linetype = "Correlation Adjusted") + 
    facet_grid(`Correlation` ~ `Set Size`, scales = "free", labeller = label_both)

pwr <- readRDS(file = "analyses/simulations_single_sample_pwr/output/sim_ss_pwr.rds")
pwr_grid <- readRDS(file = "analyses/simulations_single_sample_pwr/output/simulation_grid_pwr.rds")
pwr_grid <- pwr_grid %>% select(spar, s_rho, eff_size, id)
pwr_df <- inner_join(pwr, pwr_grid, by = "id") %>% 
    mutate(adj = replace_na(adj, "Not Applicable"), distr = replace_na(distr, "Not Applicable")) %>%
    mutate(distr = recode(distr, mnorm = "Mixture Normal", norm = "Normal"), 
           model = recode(model, cilr = "cILR", wilcox = "Wilcoxon Rank Sum"), 
           adj = recode(adj, "FALSE" = "No", "TRUE" = "Yes")) %>% 
    rename(c("Effect Size" = "eff_size", "Correlation" = "s_rho"))

pwr_plt <- ggplot(pwr_df, aes(x = spar, y = est, col = model, linetype = adj, shape = distr)) + 
    geom_linerange(aes(ymax = upper, ymin = lower)) + 
    geom_point() + 
    geom_line() + 
    geom_hline(yintercept = 0.8, col = "red") + 
    scale_color_d3() +
    theme_bw() + 
    labs(x = "Sparsity", y = "Power", col = "Model type", shape = "Distribution", 
         linetype = "Correlation Adjusted") + 
    facet_grid(`Correlation`~`Effect Size`, scales = "free", labeller = label_both)

hypo_plt <- fdr_plt + pwr_plt + plot_annotation(tag_levels = "A") + 
    plot_layout(guide = "collect") & 
    theme(legend.position = "bottom")

ggsave(hypo_plt, filename = "figures/sim_ss_hypo.png", dpi = 300, width = 12, height = 8)
file.copy(from = "figures/sim_ss_hypo.png", to = "../teailr_manuscript/manuscript/figures/sim_ss_hypo.png", overwrite = T)
