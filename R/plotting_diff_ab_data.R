library(tidyverse)
library(ggsci)

fdr_plot <- function(dir){
    data <- readRDS(file = dir)
    data <- data %>% 
        group_by(methods, output, distr, adj) %>% 
        summarise(est = mean(eval), sd = sd(eval), n = length(eval)) %>% 
        mutate(upper = est + sd/sqrt(n), lower = est - sd/sqrt(n)) %>%
        unite("model", methods, distr)  %>% 
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
        ))
    plt <- ggplot(data, aes(x = str_wrap(model, width = 15), y = est, shape = output, linetype = adj, 
                               col = str_wrap(model, width = 15))) + 
        geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 0.5)) + 
        geom_hline(yintercept = 0.05, col = "red") + scale_color_d3() + coord_flip() + 
        theme_bw() + 
        labs(y = "Type I error", x = "Methods", linetype = "Adjusted", 
             shape = "Output") + 
        guides(color = FALSE, 
               linetype = guide_legend(override.aes = list(shape = NA)))
    return(plt)
}

fdr_files <- list("analyses/data_diff_ab/output/stool_16S_fdr.rds", 
                  "analyses/data_diff_ab/output/stool_wgs_fdr.rds")

fdr_plots <- map(fdr_files, fdr_plot)
pwr <- readRDS(file = "analyses/data_diff_ab/output/gingival_16S_pwr.rds")

pwr <- pwr %>% mutate(adj = replace_na(adj, "Not Applicable"), 
                      output = case_when(
                          output == "cdf" ~ "CDF", 
                          output == "zscore" ~ "z-score", 
                          TRUE ~ "Not Applicable"
                      ))
# simulations  
results <- readRDS(file = "analyses/simulations_diff_ab/output/sim_diff_ab.rds")
grid <- readRDS(file = "analyses/simulations_diff_ab/output/sim_diff_ab_grid.rds")

results <- results %>% rowwise() %>% mutate(id = as.numeric(tail(strsplit(id, "_")[[1]],1)))

combined <- left_join(results,grid)

combined %>% group_by(model, distr, adj, output, spar, s_rho, eff_size) %>% 
    summarise(est = mean(est), upper = mean(est) + sd(est), lower = mean(est) - sd(est)) %>% 
    ggplot(aes(x = spar, y = est, col = model)) + geom_line() + facet_grid(eff_size~s_rho)
