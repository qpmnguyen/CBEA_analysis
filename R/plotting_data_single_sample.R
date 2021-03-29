library(tidyverse)
library(ggsci)
library(patchwork)

process_results <- function(df){
    df %>% mutate(distr = ifelse(is.na(distr), "None", distr), 
                  adj = ifelse(is.na(adj), FALSE, adj),
                  output = ifelse(is.na(output), "Raw Scores", output)) %>% 
        mutate(models = ifelse(models == "cilr", "cILR", models))
}

fdr <- readRDS(file = "analyses/data_single_sample/output/fdr_comparison.rds") %>%
    process_results()
fdr <- fdr %>% mutate(models = recode(models, "cILR" = "cILR", "wilcox" = "Wilcoxon Rank Sum"),
                      distr = recode(distr, "mnorm" = "Mixture Normal", "norm" = "Normal", "None" = "Not Applicable"))

 

fdr_plot <- ggplot(fdr, aes(x = models, y = est, color = adj, shape = distr)) + 
        geom_pointrange(aes(ymax = upper, ymin = lower), size = 1, position = position_dodge(width = 1)) +
        theme_bw(base_size = 15) + ggsci::scale_color_d3(labels = c("Unadjusted", "Adjusted")) + 
        geom_hline(yintercept = 0.05, col = "red") + theme(legend.position = "None") +
        labs(x = "Methods", y = "Type I error rate", shape = "Distribution", color = "Adjusted")




pwr <- readRDS(file = "analyses/data_single_sample/output/pwr_comparison.rds") %>% process_results()
pwr <- pwr %>% mutate(models = recode(models, "cILR" = "cILR", "wilcox" = "Wilcoxon Rank Sum"),
                      distr = recode(distr, "mnorm" = "Mixture Normal", "norm" = "Normal", "None" = "Not Applicable")) 
pwr_plot <- ggplot(pwr, aes(x = models, y = est, color = adj, shape = distr)) + 
    geom_pointrange(aes(ymax = upper, ymin = lower), size = 1, position = position_dodge(width = 1)) +
    theme_bw(base_size = 15) + ggsci::scale_color_d3(labels = c("Unadjusted", "Adjusted")) + 
    geom_hline(yintercept = 0.8, col = "red") + 
    labs(x = "Methods", y = "Power", shape = "Distribution", color = "Adjusted")


auc <- readRDS(file = "analyses/data_single_sample/output/auc_comparison.rds") %>% process_results()
auc <- auc %>% mutate(models = recode(models, "cILR" = "cILR", "wilcox" = "Wilcoxon Rank Sum"),
                      distr = recode(distr, "mnorm" = "Mixture Normal", "norm" = "Normal", "None" = "Not Applicable"))
auc_plot <- auc %>% group_by(models) %>% filter(est == max(est)) %>% 
    ggplot(aes(x = models, y = est, color = models)) + 
    scale_color_d3() + 
    geom_pointrange(aes(ymax = upper, ymin = lower), size = 1, position = position_dodge(width = 1)) + 
    theme_bw(base_size = 15) + geom_hline(yintercept = 0.8, col = "red") + 
    theme(legend.position = "None") +
    labs(x = "Methods", y = "AUC", col = "Models") + ylim(c(0.6,1))


layout <- "
AAAA
AAAA
#CC#
#CC#
"

top_plots <- fdr_plot + pwr_plot + plot_layout(guide = "collect") & 
    theme(legend.position = "bottom", legend.margin = margin())
comb_plot <- (top_plots / auc_plot) + 
    plot_layout(design = layout) + 
    plot_annotation(tag_levels = "A") & 
    theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin())


ggsave(comb_plot, filename = "figures/data_ss_comb.png", dpi = 300, width = 8, height = 8)
file.copy("figures/data_ss_comb.png", 
          "../teailr_manuscript/manuscript/figures/data_ss_comb.png", 
          overwrite = TRUE)


   
