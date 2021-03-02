library(tidyverse)
library(ggsci)

data <- readRDS(file = "output/data_prediction.rds")
data <- data %>% mutate(output = replace_na(output, "Raw Scores"), 
                              distr = replace_na(distr, "Not applicable"),
                              adj = replace_na(adj, "Not Applicable")) %>% 
    mutate(dset = case_when(
        dset == "nielsen_ibd_wgs" ~ "IBD vs Control - WGS",
        TRUE ~ "UC vs Control - 16S"
    )) %>% 
    mutate(method = case_when(
        method == "cilr" ~ "cILR",
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
    geom_pointrange(aes(ymax = upper, ymin = lower), size = 1, position = position_dodge(width = 1)) + 
    facet_wrap(~dset) + theme_bw(base_size = 15) + scale_color_npg() + 
    labs(x = "Models", y = "AUC", col = "Distribution", shape = "Correlation Adjustment", 
         linetype = "Output type")
saveRDS(plt, file = "output/data_prediction_plot.rds")
ggsave(plt, filename = "output/data_prediction_plot.png", dpi = 300, width = 13, height = 8)
