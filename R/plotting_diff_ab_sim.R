library(tidyverse)
library(ggsci)
library(patchwork)
library(glue)
source("R/plot_utils.R")

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}

results <- readRDS(file = "output/sim_diff_ab.rds")
grid <- readRDS(file = "output/sim_diff_ab_grid.rds")

plot_df <- results %>% mutate(id = map_dbl(id, function(x) as.numeric(tail(strsplit(x, "_")[[1]], n = 1)))) %>%
    left_join(grid) %>% group_by(models, distr, output, adj, eff_size, spar, s_rho) %>% 
    mutate(value = if_else(eff_size > 1, 1 - value, value)) %>%
    summarise(estimate = mean(value), stderr = sd(value)/sqrt(10)) %>%
    mutate(lower = estimate - stderr, upper = estimate + stderr) %>% 
    unite(col = models, c(models, distr)) %>% 
    mutate(models = case_when(
        models == "cbea_mnorm" ~ "CBEA Mixture Normal",
        models == "cbea_norm" ~ "CBEA Normal", 
        models == "deseq2_NA" ~ "DESeq2",
        models == "corncob_NA" ~ "corncob"
    )) %>% 
    mutate(output = recode(output, "cdf" = "CDF", "zscore" = "Z-scores")) %>%
    mutate(output = replace_na(output, "Not applicable")) %>%
    mutate(adj = replace_na(adj, "Not applicable")) %>% 
    rename("Sparsity" = "spar", "Correlation" = "s_rho", "Effect Size" = "eff_size")


type_i_error <- plot_df %>% filter(`Effect Size` == 1)
power <- plot_df %>% filter(`Effect Size` > 1)

# type_i_plt <- ggplot(type_i_error, aes(x = model, y = est, shape = adj, linetype = output)) + 
#     geom_pointrange(aes(ymin = lower, ymax = upper, col = model), show.legend = FALSE, position = position_dodge(width = 0.5)) +  
#     facet_grid(Sparsity ~ Correlation, labeller = label_both) + geom_hline(aes(yintercept = 0.05), col = "red") + 
#     scale_color_d3() + theme_bw() + 
#     theme(axis.title.x = element_blank()) + 
#     labs(y = "Type I error", linetype = "Output type", shape = "Correlation adjustment", col = "Model")
(type_i_plt <- ggplot(type_i_error, aes(x = Sparsity, y = estimate, 
                                        col = models, shape = adj)) + 
        geom_point() +
        geom_errorbar(aes(ymin = lower, ymax = upper, linetype = output), width = 0.05) + 
        geom_line(aes(linetype = output), show.legend = FALSE) + 
        facet_grid(~ Correlation, labeller = label_both) + 
        scale_color_d3() + my_pretty_theme + 
        guides(linetype = guide_legend(override.aes = list(shape = NA))) +
        labs(y = "Type I error", linetype = "Output type", shape = "Correlation adjustment", col = "Model") + 
        geom_hline(yintercept = 0.05, col = "red"))
ggsave(type_i_plt, filename = "figures/sim_diffab_fdr.png", 
       dpi = 300, width = 10, height = 6)
ggsave(type_i_plt, filename = "figures/sim_diffab_fdr.eps", 
       dpi = 300, width = 10, height = 6, device = cairo_ps)
# accidentally coded true positive rate as 1 - sensitivity
(pwr_plot <- ggplot(power, aes(x = Sparsity, y = estimate,  col = models, shape = adj)) +
        geom_point() +
        geom_errorbar(aes(ymin = lower, ymax = upper, linetype = output), width = 0.05) + 
        geom_line(aes(linetype = output), show.legend = FALSE) + 
        facet_grid(`Effect Size` ~ Correlation, labeller = label_both) + 
        scale_color_d3() + my_pretty_theme + 
        guides(linetype = guide_legend(override.aes = list(shape = NA))) + 
        geom_hline(yintercept = 0.8, col = "red") + 
        labs(y = "Power", linetype = "Output type", shape = "Correlation adjustment", col = "Model"))
ggsave(pwr_plot, filename = "figures/sim_diffab_pwr.png", 
       dpi = 300, width = 10, height = 6)
ggsave(pwr_plot, filename = "figures/sim_diffab_pwr.eps", 
       dpi = 300, width = 10, height = 6, device = cairo_ps)



file.copy(from = Sys.glob("figures/*.png"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.eps"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
