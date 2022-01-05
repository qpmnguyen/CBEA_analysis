# Plotting results for prediction analysis with real data  

library(tidyverse)
library(ggsci)
library(glue)

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}

data <- readRDS(file = "analyses/data_prediction/output/data_prediction.rds")
data <- data %>% mutate(output = replace_na(output, "Raw Scores"), 
                              distr = replace_na(distr, "Not Applicable"),
                              adj = replace_na(adj, "Not Applicable")) %>% 
    mutate(dset = case_when(
        dset == "nielsen_ibd_wgs" ~ "IBD vs Control - WGS",
        TRUE ~ "IBD vs Control - 16S"
    )) %>% 
    mutate(distr = case_when(
        distr == "mnorm" ~ "Mixture Normal", 
        distr == "norm" ~ "Normal",
        TRUE ~ "Not Applicable"
    )) %>% 
    mutate(method = case_when(
        method == "cilr" ~ "CBEA",
        method == "ssgsea" ~ "ssGSEA",
        method == "gsva" ~ "GSVA", 
        method == "clr" ~ "CLR"
    )) %>% 
    mutate(output = case_when(
        output == "cdf" ~ "CDF",
        output == "zscore" ~ "z-score",
        TRUE ~ output
    ))

plt <- ggplot(data, aes(x = method, y = auc, col = distr, shape = adj, linetype = output)) + 
    geom_pointrange(aes(ymax = upper, ymin = lower),position = position_dodge(width = 1)) + 
    facet_grid(dset~.) + theme_bw(base_size = 15) + scale_color_d3() + 
    labs(x = "Models", y = "AUC", col = "Distribution", shape = "Correlation Adjustment", 
         linetype = "Output type") +
    scale_linetype_manual(values = c("solid","longdash","dotted"), 
                          guide = guide_legend(override.aes = list(shape = c(NA,NA,NA)))) + 
    geom_hline(yintercept = 0.8, color = "red") + theme(legend.position = "bottom", legend.box = "vertical", 
                                                        legend.margin = margin())

#saveRDS(plt, file = "output/data_prediction_plot.rds")
ggsave(plt, filename = "figures/data_prediction_plot.png", dpi = 300, width = 8, height = 8)
ggsave(plt, filename = "figures/data_prediction_plot.eps", dpi = 300, width = 8, height = 8)

file.copy(from = "figures/data_prediction_plot.png", 
         to = glue("{save_dir}/data_prediction_plot.png", save_dir = save_dir), overwrite = TRUE)

file.copy(from = Sys.glob("figures/*.eps"), 
          to = glue("{save_dir}", save_dir = save_dir), recursive = TRUE, overwrite = TRUE)
