library(ggsci)
library(patchwork)
library(tidyverse)
library(glue)
source("R/plot_utils.R")

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}



# define some functions ####
plot_func <- function(df, sat_value, type){
    type <- match.arg(type, c("classif", "regr"))
    if (type == "classif"){
        ylab <- "AUROC"
    } else {
        ylab <- "RMSE"
    }
    
    title <- paste("Effect Saturation:", sat_value)
    
    df <- df %>% filter(sat == sat_value)
    plt <- ggplot(df, aes(x = snr, y = estimate, col = model)) + 
        geom_linerange(aes(ymax = upper, ymin = lower), show.legend = FALSE) + 
        geom_point(aes(shape = adj)) + 
        geom_line(aes(linetype = output)) + 
        facet_grid(`Corr` ~ `Sparsity`, labeller = label_both) + 
        scale_color_d3() + my_pretty_theme + 
        guides(linetype = guide_legend(override.aes = list(shape = NA, col = "black")),
               shape = guide_legend(override.aes = list(size = 3)), 
               color = guide_legend(override.aes = list(shape = NA, size = 2)))
    
    plt <- plt + labs(x = "Signal-to-noise ratio", y = ylab, 
                      col = "Models", 
                      linetype = "Output type", 
                      shape = "Correlation adjustment",
                      subtitle = title)
    return(plt)
}

process_outputs <- function(grid, results){
    df <- left_join(grid, results, by = "id") %>% 
        dplyr::select(snr, sat, spar, s_rho, models, 
                      distr, adj, output, estimate, stderr) %>% 
        unite(model, c(models, distr), sep = "_") %>%
        mutate(model = str_replace(model, "_NA", ""), 
               adj = replace_na(adj, "Not applicable"), 
               output = case_when(
                   output == "zscore" ~ "z-score",
                   output == "cdf" ~ "CDF", 
                   TRUE ~ "Not applicable"
               )) %>% 
        mutate(model = case_when(
            model == "cbea_mnorm" ~ "CBEA Mixture Normal",
            model == "cbea_norm" ~ "CBEA Normal",
            model == "gsva" ~ "GSVA", 
            model == "ssgsea" ~ "ssGSEA",
            model == "clr" ~ "CLR"
        )) %>%
        mutate(lower = estimate - stderr, upper = estimate + stderr) %>% 
        dplyr::rename("Sparsity" = "spar",  "Corr" = "s_rho")
    return(df)
}


# Process results #### 

# classif ####
classif_grid <- readRDS(file = "output/simulation_grid_classif.rds")
classif_results <- readRDS(file = "output/sim_pred_classif.rds")
classif_df <- process_outputs(classif_grid, classif_results)

classif_auc <- plot_func(classif_df, 0.1, "classif") + plot_func(classif_df, 0.5, "classif") 

# Regr ####
regr_grid <- readRDS(file = "output/simulation_grid_regr.rds")
regr_results <- readRDS(file = "output/sim_pred_regr.rds")
regr_df <- process_outputs(regr_grid, regr_results)

regr_rmse <- plot_func(regr_df, 0.1, "regr") + plot_func(regr_df, 0.5, "regr") 

combined_sim <- regr_rmse / classif_auc + plot_layout(guides = "collect") + 
    plot_annotation(tag_levels = list(c("A. Regression", "", "B. Classification", ""))) &
    theme(plot.tag = element_text(face = "bold", size = 16, hjust = 0), plot.tag.position = c(0,1.05),
          legend.position = "right", plot.margin = margin(t = 15, b = 15, l = 10))


ggsave(combined_sim, filename = "figures/sim_pred_combined.png", dpi = 300, width = 12, height = 10)
ggsave(combined_sim, filename = "figures/sim_pred_combined.eps", dpi = 300, width = 12, height = 10, device = cairo_ps)


file.copy(from = Sys.glob("figures/*.eps"), 
          to = glue("{save_dir}", save_dir = save_dir), recursive = TRUE, overwrite = TRUE)
file.copy(from = Sys.glob("figures/*.png"),
          to = glue("{save_dir}", save_dir = save_dir), recursive = TRUE, overwrite = TRUE)





