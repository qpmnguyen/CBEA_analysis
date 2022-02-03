# Make sure to set working directory to be the main working directory as teailr
library(tidyverse)
library(phyloseq)
library(ggsci)
library(patchwork)
library(glue)
source("R/plot_utils.R")


if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}

fdr <- readRDS(file = "output/sim_ss_fdr.rds")
fdr_grid <- readRDS(file = "output/simulation_grid_fdr.rds") 

fdr_grid <- fdr_grid %>% dplyr::select(spar, s_rho, n_inflate, id)

fdr_df <- inner_join(fdr, fdr_grid) %>% 
    unite("models", models, distr) %>% 
    mutate(models = case_when(
        models == "cbea_mnorm" ~ "CBEA Mixture Normal", 
        models == "cbea_norm" ~ "CBEA Normal",
        models == "wilcoxon_NA" ~ "Wilcoxon Rank Sum Test"
    )) %>% 
    mutate(adj = replace_na(adj, "Not Applicable")) %>% 
    dplyr::rename("Set Size" = "n_inflate", "Correlation" = "s_rho") 


fdr_plt <- ggplot(fdr_df, aes(x = spar, y = estimate, col = models, shape = adj)) + 
    geom_line() +
    geom_hline(yintercept = 0.05, col = "red", size = 1.25) +
    geom_pointrange(aes(ymax = upper, ymin = lower), size = 0.75) +
    facet_grid(Correlation ~ `Set Size`, labeller = label_both, scales = "free") +
    ylim(c(0,0.8)) +
    scale_color_d3() +
    labs(x = "Sparsity", y = "Type I error", col = "Models", shape = "Correlation Adjusted") +
    my_pretty_theme

fdr_plt
ggsave(fdr_plt, filename = "figures/sim_ss_fdr.png", width = 12, height = 10, units = "in")
ggsave(fdr_plt, filename = "figures/sim_ss_fdr.eps", width = 12, height = 10, units = "in", device = cairo_ps)


pwr <- readRDS(file = "output/sim_ss_pwr.rds")
pwr_grid <- readRDS(file = "output/simulation_grid_pwr.rds")
pwr_grid <- pwr_grid %>% dplyr::select(spar, s_rho, eff_size, id)
pwr_df <- inner_join(pwr, pwr_grid, by = "id") %>% 
    unite("models", models, distr) %>% 
    mutate(models = case_when(
        models == "cbea_mnorm" ~ "CBEA Mixture Normal", 
        models == "cbea_norm" ~ "CBEA Normal",
        models == "wilcoxon_NA" ~ "Wilcoxon Rank Sum Test"
    )) %>% 
    mutate(adj = replace_na(adj, "Not Applicable")) %>% 
    dplyr::rename("Effect Size" = "eff_size", "Correlation" = "s_rho") 

pwr_plt <- ggplot(pwr_df, aes(x = spar, y = estimate, col = models, shape = adj)) +
    geom_line() + 
    geom_pointrange(aes(ymax = upper, ymin = lower), size = 0.75) + 
    geom_hline(yintercept = 0.8, col = "red") + 
    scale_color_d3() +
    labs(x = "Sparsity", y = "Power", col = "Models", shape = "Correlation Adjusted", 
         linetype = "Correlation Adjusted") + 
    facet_grid(`Correlation`~`Effect Size`, scales = "free", labeller = label_both) + 
    my_pretty_theme + 
    theme(legend.position = "bottom", legend.box = "vertical")



auc <- readRDS(file = "output/sim_ss_auc.rds")
auc_grid <- readRDS(file = "output/simulation_grid_auc.rds")

auc_df <- inner_join(auc, auc_grid, by = "id") %>% 
    select(models, distr, adj, output, estimate, lower, upper, spar, s_rho, eff_size) %>% 
    unite("models", c(models, distr)) %>% 
    mutate(adj = replace_na(adj, "Not Applicable"), output = replace_na(output, "Not Applicable")) %>% 
    mutate(output = recode(output, zscore = "Z-score", cdf = "CDF")) %>% 
    mutate(models = recode(models, cbea_mnorm = "CBEA Mixture Normal", cbea_norm = "CBEA Normal", 
                           ssgsea_NA = "ssGSEA", gsva_NA = "GSVA", wilcoxon_NA = "Wilcoxon W-Statistic")) %>% 
    rename("Correlation" = "s_rho", "Effect Size" = "eff_size")
     
auc_plt <- ggplot(auc_df, aes(x = spar, y = estimate, col = models, shape = adj)) + 
    geom_line() +
    geom_pointrange(aes(ymax = upper, ymin = lower), size = 0.75) + 
    geom_hline(yintercept = 0.8, col = "red") + 
    scale_color_d3() +
    theme_bw() + 
    labs(x = "Sparsity", y = "AUC", col = "Models", shape = "Correlation Adjusted") + 
    facet_grid(`Correlation`~`Effect Size`, scales = "fixed", labeller = label_both) + 
    my_pretty_theme +
    theme(legend.position = "bottom", legend.box = "vertical")

sim_pheno_relevance <- pwr_plt + auc_plt + plot_annotation(tag_levels = "A")
ggsave(sim_pheno_relevance, filename = "figures/sim_ss_pwr.png", dpi = 300, width = 15, height = 10)
ggsave(sim_pheno_relevance, filename = "figures/sim_ss_pwr.eps", dpi = 300, width = 15, height = 10, device = cairo_ps)

file.copy(from = Sys.glob("figures/*.eps"), 
          to = glue("{save_dir}", save_dir = save_dir), recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.png"), 
          to = glue("{save_dir}", save_dir = save_dir), recursive = TRUE, overwrite = TRUE)
