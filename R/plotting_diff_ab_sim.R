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


combo_plot <- (type_i_plt/pwr_plot) + 
    plot_annotation(tag_levels = list(c("A. Type I error", "B. Power"))) + 
    plot_layout(guides = "collect", heights = c(1,2)) & 
    guides(linetype = guide_legend(override.aes = list(shape = NA)), 
           color = guide_legend(keyheight = 0.4, keywidth = 0.2, 
                                default.unit = "inch", 
                                override.aes = list(linetype = 0)), 
           shape = guide_legend(override.aes = list(linetype = 0))) &
    theme(plot.tag = element_text(face = "bold", size = 16, hjust = 0), 
          plot.tag.position = c(0,1.05), 
          plot.margin = margin(t = 22, b = 15, l = 10))


ggsave(combo_plot, filename = "figures/sim_diff_ab_comb.png", 
       dpi = 300, width = 12, height = 12)
ggsave(combo_plot, filename = "figures/sim_diff_ab_comb.eps", 
       dpi = 300, width = 15, height = 12, device = cairo_ps)

file.copy(from = Sys.glob("figures/*.png"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.eps"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
