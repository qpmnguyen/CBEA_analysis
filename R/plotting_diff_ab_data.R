library(tidyverse)
library(ggsci)
library(glue)
library(patchwork)

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}


# FALSE DISCOVERY RATE ####
fdr <- readRDS(file = "output/data_diffab_fdr.rds")
fdr_df <- fdr %>% group_by(models, distr, adj, output) %>% 
    summarise(mean = mean(res, na.rm = TRUE), stderr = sd(res, na.rm = TRUE)/sqrt(500)) %>%
    unite(col = "models", c(models, distr)) %>% 
    mutate(models = case_when(
        models == "cbea_mnorm" ~ "CBEA Mixture Normal",
        models == "cbea_norm" ~ "CBEA Normal",
        models == "deseq2_NA" ~ "DESeq2",
        models == "corncob_NA" ~ "corncob"
    )) %>% 
    mutate(adj = replace_na(adj, "Not Applicable")) %>%
    mutate(output = case_when(
        output == "raw" ~ "Raw scores",
        output == "cdf" ~ "CDF", 
        output == "zscore" ~ "Z-scores",
        TRUE ~ "Not Applicable"
    )) %>% 
    mutate(models = if_else(output == "Raw scores", "CBEA Raw Scores", models)) %>% 
    mutate(output = if_else(output == "Raw scores", "Not Applicable", output)) %>% 
    mutate(upper = mean + stderr, lower = mean - stderr) 

fdr_plt <- ggplot(fdr_df, aes(x = str_wrap(models, width = 10), y = mean, col = output, shape = adj)) + 
    geom_point(position = position_dodge(width = 1)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper), size = 1, 
                    position = position_dodge(width = 1)) +
    geom_hline(yintercept = 0.05, col = "red") + 
    labs(x = "Type I error", y = "Models", shape = "Correlation adjusted", col = "Output type") +
    coord_flip() + 
    scale_color_d3() +
    my_pretty_theme +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical",
          legend.margin = margin(0,0,0,0), plot.margin = unit(c(0,0,30,0), "pt"))

rset <- readRDS(file = "output/data_diffab_rset.rds")

rset_df <- rset %>% unite("models", c(models, distr)) %>% 
    mutate(models = case_when(
        models == "cbea_mnorm" ~ "CBEA Mixture Normal",
        models == "cbea_norm" ~ "CBEA Normal",
        models == "deseq2_NA" ~ "DESeq2",
        models == "corncob_NA" ~ "corncob"
    )) %>% 
    mutate(adj = replace_na(adj, "Not Applicable")) %>%
    mutate(output = case_when(
        output == "raw" ~ "Raw scores",
        output == "cdf" ~ "CDF", 
        output == "zscore" ~ "Z-scores",
        TRUE ~ "Not Applicable"
    )) %>% 
    mutate(models = if_else(output == "Raw scores", "CBEA Raw Scores", models)) %>% 
    mutate(output = if_else(output == "Raw scores", "Not Applicable", output)) %>% 
    mutate(method = if_else(str_detect(models, "CBEA"), "CBEA", models))

rset_plt <- ggplot(rset_df, aes(x = size, y = estimate, col = models, shape = output, 
                                 linetype = adj)) + 
        geom_point(size = 4) +
        geom_pointrange(aes(ymin = lower, ymax = upper), show.legend = FALSE) +
        geom_line(size = 1) +
        guides(linetype = guide_legend(override.aes = list(shape = NA))) +
        facet_wrap(~method) +
        labs(x = "Set size", y = "Proportion of significant sets", 
             col = "Models", shape = "Output type", linetype = "Correlation adjusted") +
        ylim(c(0,1)) +
        my_pretty_theme + 
        theme(legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical", 
              legend.margin = margin(0,0,0,0), plot.margin = unit(c(30,0,0,0), "pt"))

combo_diff_ab <- fdr_plt/rset_plt + plot_annotation(tag_levels = "A") + plot_layout(heights = c(0.6,0.4))

ggsave(combo_diff_ab, filename = "figures/data_diffab_fdr.png", dpi = 300, width = 10, height = 10)
ggsave(combo_diff_ab, filename = "figures/data_diffab_fdr.eps", dpi = 300, width = 10, height = 10, device = cairo_ps)

# PWR ANALYSES ####  
pwr <- readRDS(file = "output/data_diffab_pwr.rds")
pwr_plotting <- pwr %>% unite("models", c(models, distr)) %>% 
    mutate(models = case_when(
        models == "cbea_mnorm" ~ "CBEA Mixture Normal",
        models == "cbea_norm" ~ "CBEA Normal",
        models == "deseq2_NA" ~ "DESeq2",
        models == "corncob_NA" ~ "corncob"
    )) %>% 
    mutate(adj = replace_na(adj, "Not Applicable")) %>%
    mutate(output = case_when(
        output == "raw" ~ "Raw scores",
        output == "cdf" ~ "CDF", 
        output == "zscore" ~ "Z-scores",
        TRUE ~ "Not Applicable"
    )) %>% 
    mutate(models = if_else(output == "Raw scores", "CBEA Raw Scores", models)) %>% 
    mutate(output = if_else(output == "Raw scores", "Not Applicable", output)) %>% 
    mutate(method = if_else(str_detect(models, "CBEA"), "CBEA", models))


pwr_plt <- ggplot(pwr_plotting, aes(y = estimate, x = str_wrap(models,10), col = output, shape = adj)) + 
    geom_pointrange(aes(ymin = lower, ymax = upper), size = 1, 
                    position = position_dodge(width = 1)) + 
    labs(y = "Power", x = "Models", col = "Output type", shape = "Correlation adjusted") + 
    coord_flip() +
    scale_color_d3() +
    ylim(c(0,1)) +
    my_pretty_theme

ggsave(pwr_plt, filename = "figures/data_diffab_pwr.png", dpi = 300, width = 10, height = 5.5)
ggsave(pwr_plt, filename = "figures/data_diffab_pwr.eps", dpi = 300, width = 10, height = 5.5, device = cairo_ps)

file.copy(from = Sys.glob("figures/*.png"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.eps"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)

