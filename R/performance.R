library(tidyverse)
library(ggsci)
library(patchwork)
library(glue)
results <- readRDS(file = "analyses/performance/output/results.rds")
grid <- readRDS(file = "analyses/performance/output/grid_sim.rds")

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}

data <- left_join(grid, results)

data <- data %>% mutate()

samp_dat <- data %>% filter(eval == "samp")

samp_plot <- ggplot(samp_dat, aes(x = n_samp, y = time, col = distr, shape = adj)) + 
    geom_point(size = 3) + geom_line() + scale_color_d3() + theme_bw() + 
    labs(x = "Number of samples", y = "Time (seconds)", color = "Distribution", 
         shape = "Correlation adjustment", 
         title = "Varying number of samples \n(10 sets) ") + 
    ylim(c(0,310))

set_dat <- data %>% filter(eval == "tax") %>% mutate(n_sets = n_tax/50)
set_plot <- ggplot(set_dat, aes(x = n_sets, y = time, col = distr, shape = adj)) + 
    geom_point(size = 3) + 
    geom_line() + scale_color_d3() + theme_bw() + 
    labs(x = "Number of sets", y = "Time (seconds)", color = "Distribution", 
         shape = "Correlation adjustment", 
         title = "Varying number of sets \n(N = 1000)") + 
    theme(axis.title.y = element_blank()) + ylim(c(0,310))

combo_plot <- samp_plot + set_plot + plot_layout(guides = "collect")
ggsave(combo_plot, filename = "figures/performance.png", dpi = 300, width = 10, height = 5)
ggsave(combo_plot, filename = "figures/performance.eps", dpi = 300, width = 10, height = 5)

file.copy(from = Sys.glob("figures/*.png"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.eps"), to = glue("{save_dir}", dir = save_dir), 
          recursive = TRUE, overwrite = TRUE)
