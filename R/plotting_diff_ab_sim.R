library(tidyverse)
library(ggsci)
library(patchwork)
library(glue)

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}

results <- readRDS(file = "analyses/simulations_diff_ab/output/sim_diff_ab.rds")
grid <- readRDS(file = "analyses/simulations_diff_ab/output/sim_diff_ab_grid.rds")

results <- results %>% mutate(id = map(id, ~{str_split(.x, "_")[[1]][3]})) %>% 
    unnest(id) %>% mutate(id = as.numeric(id))

combined_df <- left_join(results, grid)
combined_df <- combined_df %>% group_by(model, distr, adj, output, spar, eff_size, s_rho) %>% 
    summarise(est = mean(res), 
              upper = mean(res) + (sd(res)/sqrt(10)), 
              lower = mean(res) - (sd(res)/sqrt(10))) %>% 
    unite("model", model, distr)  %>% 
    mutate(model = case_when(
        model == "cilr_wilcox_mnorm" ~ "cILR Mixture Normal w/ Wilcox Test", 
        model == "cilr_wilcox_norm" ~ "cILR Normal w/ Wilcox Test",
        model == "cilr_welch_mnorm" ~ "cILR Mixture Normal w/ Welch Test", 
        model == "cilr_welch_norm" ~ "cILR Normal w/ Welch Test",
        model == "deseq2_NA" ~ "DESeq2", 
        model == "corncob_NA" ~ "corncob"
    )) %>% 
    mutate(adj = replace_na(adj, "Not Applicable")) %>% 
    mutate(output = case_when(
        output == "cdf" ~ "CDF", 
        output == "zscore" ~ "z-score", 
        TRUE ~ "Not Applicable"
    )) %>%
    mutate(model = str_wrap(model, width = 20)) %>% 
    rename("Sparsity" = "spar", "Correlation" = "s_rho", "Effect Size" = "eff_size")

type_i_error <- combined_df %>% filter(`Effect Size` == 1)
power <- combined_df %>% filter(`Effect Size` > 1)

# type_i_plt <- ggplot(type_i_error, aes(x = model, y = est, shape = adj, linetype = output)) + 
#     geom_pointrange(aes(ymin = lower, ymax = upper, col = model), show.legend = FALSE, position = position_dodge(width = 0.5)) +  
#     facet_grid(Sparsity ~ Correlation, labeller = label_both) + geom_hline(aes(yintercept = 0.05), col = "red") + 
#     scale_color_d3() + theme_bw() + 
#     theme(axis.title.x = element_blank()) + 
#     labs(y = "Type I error", linetype = "Output type", shape = "Correlation adjustment", col = "Model")
type_i_plt <- ggplot(type_i_error, aes(x = Sparsity, y = est, col = model, shape = adj)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper, linetype = output)) + 
    geom_line(aes(linetype = output), show.legend = FALSE) + 
    facet_grid(~ Correlation, labeller = label_both) + 
    scale_color_d3() + theme_bw() + 
    ylim(c(0.03, 0.15)) + 
    labs(y = "Type I error", linetype = "Output type", shape = "Correlation adjustment", col = "Model") + 
    geom_hline(yintercept = 0.05, col = "red")


ggsave(type_i_plt, filename = "figures/sim_diff_ab_type_i_error.png", dpi = 300, width = 15, height = 8)
ggsave(type_i_plt, filename = "figures/sim_diff_ab_type_i_error.eps", dpi = 300, width = 15, height = 8)

pwr_plot <- ggplot(power, aes(x = Sparsity, y = est,  col = model, shape = adj)) +
    geom_pointrange(aes(ymin = lower, ymax = upper, linetype = output)) + 
    geom_line(aes(linetype = output), show.legend = FALSE) + 
    facet_grid(`Effect Size` ~ Correlation, labeller = label_both) + 
    scale_color_d3() + theme_bw() + geom_hline(yintercept = 0.8, col = "red") + 
    labs(y = "Power", linetype = "Output type", shape = "Correlation adjustment", col = "Model")
ggsave(pwr_plot, filename = "figures/sim_diff_ab_pwr.png", dpi = 300, width = 10, height = 8)
ggsave(pwr_plot, filename = "figures/sim_diff_ab_pwr.eps", dpi = 300, width = 10, height = 8)

# this plot has legend going horizontal at the bottom
# combo_plot <- (type_i_plt / pwr_plot) + 
#     plot_annotation(tag_levels = list(c("A. Type I error", "B. Power"))) + 
#     plot_layout(guides = "collect", heights = c(1,2)) & 
#     guides(linetype = guide_legend(override.aes = list(shape = NA)), 
#            color = guide_legend(override.aes = list(linetype = 0)), 
#            shape = guide_legend(override.aes = list(linetype = 0))) &
#     theme(plot.tag = element_text(face = "bold", size = 16, hjust = 0), 
#           legend.position = "bottom", legend.box = "horizontal", 
#           legend.direction = "horizontal", legend.justification = "center", 
#           legend.box.just = "center", 
#           legend.box.margin = margin(r = -12, l = -12, unit = "cm"), 
#           plot.tag.position = c(0,1.05), plot.margin = margin(t = 22, b = 15, l = 10))

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
       dpi = 300, width = 15, height = 12)

file.copy(from = Sys.glob("figures/*.png"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.eps"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
