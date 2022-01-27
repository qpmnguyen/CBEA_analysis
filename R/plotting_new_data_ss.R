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

my_pretty_theme <- theme_minimal(base_family = "Lato", base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          # Bold, bigger title
          #plot.title = element_text(face = "bold", size = rel(1.7)),
          # Plain, slightly bigger subtitle that is grey
          plot.subtitle = element_text(face = "plain", size = rel(1.3), color = "grey70"),
          # Italic, smaller, grey caption that is left-aligned
          plot.caption = element_text(face = "italic", size = rel(0.7), 
                                      color = "grey70", hjust = 0),
          # Bold legend titles
          legend.title = element_text(face = "bold"),
          # Bold, slightly larger facet titles that are left-aligned for the sake of repetition
          strip.text = element_text(face = "bold", size = rel(1.1), hjust = 0),
          # Add a light grey background to the facet titles, with no borders
          strip.background = element_rect(fill = "grey90", color = NA),
          # Add a thin grey border around all the plots to tie in the facet titles
          panel.border = element_rect(color = "grey90", fill = NA))


### FDR PLOTS NEW #### 
df_fdr_new <- readRDS(file = "output/fdr_ss_randset.rds")

# diagnostic distribution plots 
# df_fdr_new %>%  mutate(adj = if_else(adj, "Adjusted", "Not Adjusted")) %>% 
#     mutate(adj = if_else(is.na(adj), "Not applicable", adj), 
#            models = if_else(is.na(models), "Not applicable", models)) %>%
#     unite("models", c(models, distr, adj)) %>% 
#     mutate(models = if_else(str_detect(models, "_NA"), 
#                             str_remove_all(models, "_NA"), models)) %>% 
#     ggplot(aes(x = res, y = models, fill = models)) + 
#     geom_density_ridges() + 
#     facet_wrap(~size, labeller = label_both) + 
#     labs(x = "Type I error", y = "Counts") + theme_bw() +
#     scale_fill_d3()

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



fdr_new <- ggplot(df_fdr_new, aes(x = size, y = estimate, col = models, shape = adj, fill = models)) +
    geom_point(size = 3) +
    geom_line(alpha = 0.5) + 
    geom_pointrange(aes(ymax = upper, ymin = lower)) + 
    scale_color_d3() +
    scale_fill_d3() +
    labs(x = "Set size", y = "Type I error", 
         col = "Models", shape = "Correlation adjusted", fill = "Models") + 
    geom_hline(yintercept = 0.05, col = "red") + my_pretty_theme


ggsave(fdr_new, filename = "figures/data_ss_fdr_new.png", dpi = 300, width = 10, height = 6)
ggsave(fdr_new, filename = "figures/data_ss_fdr_new.eps", dpi = 300, width = 10, height = 6, device = cairo_ps)
file.copy(from = Sys.glob("figures/*.png"), 
          to = glue("{save_dir}", save_dir = save_dir), overwrite = TRUE)

file.copy(from = Sys.glob("figures/*.eps"), 
          to = glue("{save_dir}", save_dir = save_dir), recursive = TRUE, overwrite = TRUE)




# PWR PLOTS NEW ####
col_values <- pal_d3(palette = "category10")(3)
names(col_values) <- c("CBEA Mixture Normal", "CBEA Normal", "Wilcoxon Rank Sum Test")

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

pwr_new <- ggplot(df_pwr_new, aes(x = str_wrap(models, width = 20), y = estimate, col = models, shape = adj)) + 
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 0.5)) + 
    labs(x = "Models", y = "Power", shape = "Correlation adjusted", col = "Models") + 
    scale_color_manual(values = col_values[1:3]) + guides(color = "none") + 
    my_pretty_theme


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
        models == "wilcoxon" ~ str_wrap("Wilcoxon Rank Sum W statistic", width = 50),
        models == "ssgsea" ~ str_wrap("ssGSEA", width = 90),
        models == "gsva" ~ str_wrap("GSVA", width = 70) 
    ))


col_values <- pal_d3(palette = "category10")(5)
names(col_values) <- c("CBEA Mixture Normal", "CBEA Normal", "Wilcoxon Rank Sum W statistic", "GSVA", "ssGSEA")

auc_new <- ggplot(df_auc_new, aes(x = str_wrap(models, width = 20), y = estimate, col = models, shape = adj)) + 
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 0.5)) + 
    labs(x = "Models", y = "AUROC", shape = "Correlation adjusted", col = "Models") + 
    scale_color_manual(values = col_values) + guides(color = "none") + 
    my_pretty_theme + theme(axis.title.y = element_blank())

combo_rel <- pwr_new + auc_new + plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = "A") & 
    coord_flip() &
    theme(plot.tag = element_text(face = "bold"))

ggsave(combo_rel, filename = "figures/data_ss_pwr_new.png", dpi = 300, width = 10, height = 6)
ggsave(combo_rel, filename = "figures/data_ss_pwr_new.eps", dpi = 300, width = 10, height = 6, device = cairo_ps)
file.copy(from = Sys.glob("figures/*.png"), 
          to = glue("{save_dir}", save_dir = save_dir), overwrite = TRUE)

file.copy(from = Sys.glob("figures/*.eps"), 
          to = glue("{save_dir}", save_dir = save_dir), recursive = TRUE, overwrite = TRUE)
