library(tidyverse)
library(ggsci)
library(patchwork)
library(glue)
library(bench)
library(ggrepel)
source("R/plot_utils.R")

results <- readRDS(file = "output/runtime.rds")

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}

results <- results %>% mutate(distr = recode(distr, "mnorm" = "Mixture Normal", "norm" = "Normal")) %>% 
    mutate(n_perm = as.factor(n_perm)) %>% 
    dplyr::rename("Permutations" = "n_perm")


perf_plot <- ggplot(results, aes(x = distr, y = time, col = distr, shape = adj, label = as.character(time))) + 
    geom_point(size = 3, position = position_dodge(width = 0.5)) + 
    geom_pointrange(ymin = 0, aes(ymax = time), position = position_dodge(width = 0.5)) + 
    geom_label_repel(position = position_dodge(width = 0.5)) +
    facet_wrap(~Permutations, labeller = label_both) +
    coord_flip() +
    labs(x = "Distribution", y = "Total runtime", col = "Distribution", shape = "Correlation Adjusted") +
    my_pretty_theme

ggsave(perf_plot, filename = "figures/performance.png", dpi = 300, width = 8, height = 5)
ggsave(perf_plot, filename = "figures/performance.eps", dpi = 300, width = 8, height = 5, device = cairo_ps)

file.copy(from = Sys.glob("figures/*.png"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.eps"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
