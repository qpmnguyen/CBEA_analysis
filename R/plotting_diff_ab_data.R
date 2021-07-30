library(tidyverse)
library(ggsci)
library(glue)

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}


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
        )) %>%
        mutate(model = str_wrap(model, width = 20))
    
    plt <- ggplot(data, aes(x = model, y = est, linetype = output, shape = adj, col = model)) + 
        geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 0.75)) + 
        geom_hline(yintercept = 0.05, col = "red") + 
        scale_color_manual(values = pal_d3("category10")(6)) + coord_flip() + 
        theme_bw() + 
        labs(y = "Type I error", shape = "Correlation adjustment", 
             linetype = "Output type", subtitle = "16S rRNA sequencing of stool samples") + 
        guides(color = FALSE, 
               linetype = guide_legend(override.aes = list(shape = NA))) +
        theme(axis.title.y = element_blank()) + 
        ylim(c(0,0.3)) + 
        annotate(geom = "text", x = 0.90, y = 0.15, label = expression(paste(alpha,"= 0.05"))) +
        annotate(geom = "curve", x = 1, xend = 1, y = 0.15, yend = 0.05, 
                 curvature = 0.2, arrow = arrow(length = unit(2, "mm")))
    
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
                   TRUE ~ "Not Applicable")) %>% 
    unite("model", methods, distr)  %>% 
    mutate(model = case_when(
        model == "cilr_wilcox_mnorm" ~ "cILR Mixture Normal w/ Wilcox Test", 
        model == "cilr_wilcox_norm" ~ "cILR Normal w/ Wilcox Test",
        model == "cilr_welch_mnorm" ~ "cILR Mixture Normal w/ Welch Test", 
        model == "cilr_welch_norm" ~ "cILR Normal w/ Welch Test",
        model == "deseq2_NA" ~ "DESeq2", 
        model == "corncob_NA" ~ "corncob"
    )) %>% 
    mutate(model = str_wrap(model, width = 20))

pwr_plots <- ggplot(pwr, aes(x = model, y = est, linetype = output, shape = adj, col = model)) + 
    geom_point(position = position_dodge(width = 0.75), show.legend = FALSE) +
    geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 0.75)) + 
    scale_color_manual(values = pal_d3("category10")(6), guide = FALSE) + theme_bw() +
    guides(linetype = guide_legend(override.aes = list(shape = NA))) + coord_flip() + 
    theme(axis.title.y = element_blank()) + 
    labs(y = "True Positive Rate", shape = "Correlation adjustment", linetype = "Output type", 
         subtitle = "16S rRNA sequencing of the gingival site") +
    ylim(y = c(0,1))


combo_diff_ab <- fdr_plots[[1]] + (pwr_plots & theme(axis.text.y = element_blank())) + 
    plot_layout(guides = "collect")

ggsave(combo_diff_ab, filename = "figures/data_diff_ab.png", dpi = 300, width = 10, height = 5.5)
ggsave(combo_diff_ab, filename = "figures/data_diff_ab.eps", dpi = 300, width = 10, height = 5.5)

file.copy(from = Sys.glob("figures/*.png"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.eps"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)

