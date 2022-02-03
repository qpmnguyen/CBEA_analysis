# Plotting results for prediction analysis with real data  
library(tidyverse)
library(ggsci)
library(glue)
source("R/plot_utils.R")



if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}

data <- readRDS(file = "output/data_pred.rds")
data <- data %>% unite("models", c(models, distr)) %>% 
    mutate(upper = estimate + std_err, lower = estimate - std_err) %>% 
    mutate(output = replace_na(output, "Not Applicable"), adj = replace_na(adj, "Not Applicable")) %>%
    mutate(models = recode(models, cbea_mnorm = "CBEA Mixture Normal", cbea_norm = "CBEA Normal",
                           ssgsea_NA = "ssGSEA", gsva_NA = "GSVA", clr_NA = "CLR")) %>% 
    mutate(models = if_else(output == "raw", "CBEA Raw Scores", models)) %>% 
    mutate(output = recode(output, raw = "Raw Scores", zscore = "Z-score", cdf = "CDF")) %>%
    mutate(data = recode(data, `16s` = "16S rRNA gene sequencing", wgs = "Whole genome sequencing")) 
    


plt <- ggplot(data, aes(x = models, y = estimate, col = models, shape = adj, linetype = output)) + 
    geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 1)) + 
    facet_grid(data~.) + theme_bw(base_size = 15) + scale_color_d3() + 
    labs(x = "Models", y = "AUROC", shape = "Correlation Adjustment", 
         linetype = "Output type") +
    guides(col = "none") +
    #scale_linetype_manual(values = c("solid","longdash","dotted"), 
    #                     guide = guide_legend(override.aes = list(shape = c(NA,NA,NA)))) + 
    geom_hline(yintercept = 0.8, color = "red") + 
    my_pretty_theme +
    theme(legend.position = "bottom", legend.box = "vertical", 
          legend.margin = margin()) 

plt

#saveRDS(plt, file = "output/data_prediction_plot.rds")
ggsave(plt, filename = "figures/data_pred.png", dpi = 300, width = 9, height = 8)
ggsave(plt, filename = "figures/data_pred.eps", dpi = 300, width = 9, height = 8, device = cairo_ps)


file.copy(from = Sys.glob("figures/*.eps"), 
          to = glue("{save_dir}", save_dir = save_dir), recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.png"), 
          to = glue("{save_dir}", save_dir = save_dir), recursive = TRUE, overwrite = TRUE)
