library(tidyverse)
library(ggsci)
library(patchwork)
library(glue)
library(ggridges)

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}

### FDR PLOTS NEW #### 
df_fdr_new <- readRDS(file = "output/fdr_ss_randset.rds")

# diagnostic distribution plots 
df_fdr_new %>%  mutate(adj = if_else(adj, "Adjusted", "Not Adjusted")) %>% 
    mutate(adj = if_else(is.na(adj), "Not applicable", adj), 
           models = if_else(is.na(models), "Not applicable", models)) %>%
    unite("models", c(models, distr, adj)) %>% 
    mutate(models = if_else(str_detect(models, "_NA"), 
                            str_remove_all(models, "_NA"), models)) %>% 
    ggplot(aes(x = res, y = models, fill = models)) + 
    geom_density_ridges() + 
    facet_wrap(~size, labeller = label_both) + 
    labs(x = "Type I error", y = "Counts") + theme_bw() +
    scale_fill_d3()

df_fdr_new <- df_fdr_new %>% mutate(adj = if_else(adj, "Yes", "No")) %>% 
    mutate(adj = if_else(is.na(adj), "Not applicable", adj), 
           models = if_else(is.na(models), "Not applicable", models)) %>%
    unite("models", models:distr) %>% 
    mutate(models = if_else(str_detect(models, "_NA"), 
                            str_remove_all(models, "_NA"), models)) %>% 
    group_by(models, adj, size) %>% 
    summarise(estimate = mean(res), 
              se = sd(res)/sqrt(1000)) %>% 
    mutate(upper = estimate + se, lower = estimate - se) %>% ungroup() %>% 
    mutate(models = case_when(
        models == "cbea_mnorm" ~ str_wrap("CBEA Mixture Normal", width = 70),
        models == "cbea_norm" ~ str_wrap("CBEA Normal", width = 70),
        models == "wilcoxon" ~ str_wrap("Wilcoxon Rank Sum Test", width = 70)
    ))



fdr_new <- ggplot(df_fdr_new, aes(x = models, y = estimate, col = models, shape = adj)) +
    geom_point(position = position_dodge(width = 1)) +
    scale_color_d3() +
    geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 1)) + 
    facet_wrap(~ size, dir = "v", labeller = label_both) + 
    theme_bw() + labs(x = "Models", y = "False Discovery Rate", 
                      col = "Models", shape = "Correlation adjusted")

fdr_wgs_new <- ggplot(df_fdr_new %>% filter(type == "wgs"), aes(x = models, y = estimate, col = models, shape = adj)) +
    geom_point(position = position_dodge(width = 1)) +
    scale_color_d3() +
    geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 1)) + 
    facet_wrap(~ size, dir = "v", labeller = label_both) + 
    theme_bw() + labs(x = "Models", y = "False Discovery Rate", 
                      col = "Models", shape = "Correlation adjusted", 
                      title = "WGS Dataset")

fdr_new <- fdr_16s_new + fdr_wgs_new + plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A")

df_pwr_new <- readRDS(file = "output/pwr_ss_pheno.rds") %>%
    mutate(adj = if_else(adj, "Yes", "No")) %>% 
    mutate(adj = if_else(is.na(adj), "Not applicable", adj), 
           models = if_else(is.na(models), "Not applicable", models)) %>%
    unite("models", models:distr) %>% 
    mutate(models = if_else(str_detect(models, "_NA"), 
                            str_remove_all(models, "_NA"), models)) %>% 
    mutate(models = case_when(
        models == "cbea_mnorm" ~ str_wrap("CBEA Mixture Normal", width = 70),
        models == "cbea_norm" ~ str_wrap("CBEA Normal", width = 70),
        models == "wilcoxon" ~ str_wrap("Wilcoxon Rank Sum Test", width = 70)
    ))

pwr_new <- ggplot(df_pwr_new, aes(x = models, y = estimate, col = models, shape = adj)) + 
    geom_point(position = position_dodge(width = 1)) + 
    geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 1)) + 
    theme_bw() +
    labs(x = "Models", y = "Power", shape = "Correlation adjusted", col = "Models") + 
    scale_color_d3()


df_auc_new <- readRDS(file = "output/auc_ss_pheno.rds") %>% 
    mutate(adj = if_else(adj, "Yes", "No")) %>% 
    mutate(adj = if_else(is.na(adj), "Not applicable", adj), 
           models = if_else(is.na(models), "Not applicable", models)) %>%
    unite("models", models:distr) %>% 
    mutate(models = if_else(str_detect(models, "_NA"), 
                            str_remove_all(models, "_NA"), models)) %>% 
    mutate(models = case_when(
        models == "cbea_mnorm" ~ str_wrap("CBEA Mixture Normal", width = 70),
        models == "cbea_norm" ~ str_wrap("CBEA Normal", width = 70),
        models == "wilcoxon" ~ str_wrap("Wilcoxon Rank Sum Test", width = 70),
        models == "ssgsea" ~ str_wrap("ssGSEA", width = 90),
        models == "gsva" ~ str_wrap("GSVA", width = 70) 
    ))

auc_new <- ggplot(df_auc_new, aes(x = models, y = estimate, col = models, shape = adj)) + 
    geom_point(position = position_dodge(width = 1)) + 
    geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 1)) + 
    theme_bw() +
    labs(x = "Models", y = "AUC", shape = "Correlation adjusted", col = "Models") + 
    scale_color_d3()

pwr_new + auc_new
