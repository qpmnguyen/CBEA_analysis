library(tidyverse)
library(ggsci)
library(patchwork)

process_results <- function(df){
    df %>% mutate(distr = ifelse(is.na(distr), "None", distr), 
                  adj = ifelse(is.na(adj), FALSE, adj),
                  output = ifelse(is.na(output), "Raw Scores", output)) %>% 
        mutate(models = ifelse(models == "cilr", "cILR", models))
}

fdr <- readRDS(file = "output/fdr_comparison.rds") %>%
    process_results()


 

fdr_plot <- ggplot(fdr, aes(x = models, y = est, color = adj, shape = distr)) + 
        geom_pointrange(aes(ymax = upper, ymin = lower), size = 1, position = position_dodge(width = 1)) +
        theme_bw(base_size = 15) + ggsci::scale_color_npg(labels = c("Unadjusted", "Adjusted")) + 
        scale_shape_discrete(labels = c("Mixture Normal", "None", "Normal")) + 
        geom_hline(yintercept = 0.05, col = "red") + theme(legend.position = "None") +
        labs(x = "Methods", y = "False Discovery Rate", shape = "Distribution", color = "Adjusted")




pwr <- readRDS(file = "output/pwr_comparison.rds") %>% process_results()

pwr_plot <- ggplot(pwr, aes(x = models, y = est, color = adj, shape = distr)) + 
    geom_pointrange(aes(ymax = upper, ymin = lower), size = 1, position = position_dodge(width = 1)) +
    theme_bw(base_size = 15) + ggsci::scale_color_npg(labels = c("Unadjusted", "Adjusted")) + 
    scale_shape_discrete(labels = c("Mixture Normal", "None", "Normal")) + 
    geom_hline(yintercept = 0.8, col = "red") + 
    labs(x = "Methods", y = "Power", shape = "Distribution", color = "Adjusted")



auc <- readRDS(file = "output/auc_comparison.rds") %>% process_results()

auc_plot <- auc %>% group_by(models) %>% filter(est == max(est)) %>% 
    ggplot(aes(x = models, y = est, color = models)) + 
    geom_pointrange(aes(ymax = upper, ymin = lower), size = 1, position = position_dodge(width = 1)) + 
    theme_bw(base_size = 15) + geom_hline(yintercept = 0.8, col = "red") + 
    theme(legend.position = "None") +
    labs(x = "Methods", y = "AUC")

comb_plot <- (fdr_plot + pwr_plot) / (plot_spacer() + auc_plot + plot_spacer())

saveRDS(comb_plot, file = "output/data_single_sample.rds")
ggsave(comb_plot, filename = "output/data_single_sample.png", dpi = 300, width = 10, height = 8)    
