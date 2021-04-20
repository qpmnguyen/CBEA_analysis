library(ggsci)
library(patchwork)
library(tidyverse)

# define some functions ####
plot_func <- function(metrics, sat_value, df_main){
    df <- df_main %>% filter(.metric == metrics) %>% filter(sat == sat_value)
    plt <- ggplot(df, aes(x = snr, y = mean, col = model)) + 
        geom_linerange(aes(ymax = upper, ymin = lower), show.legend = FALSE) + 
        geom_point(aes(shape = adj)) + 
        geom_line(aes(linetype = output)) + 
        facet_grid(`Correlation` ~ `Sparsity`, labeller = label_both) + 
        scale_color_d3() + theme_bw() + 
        guides(linetype = guide_legend(override.aes = list(shape = NA, col = "black")),
               shape = guide_legend(override.aes = list(size = 3)), 
               color = guide_legend(override.aes = list(shape = NA, size = 2)))
    if (metrics == "accuracy"){
        ylab <- "Accuracy"
    } else if (metrics == "roc_auc") {
        ylab <- "AUC"
    } else if (metrics == "rmse"){
        ylab <- "RMSE"
    } else if (metrics == "rsq"){
        ylab <- "R-squared"
    }
    
    if(sat_value == 0.1){
        title <- "Effect Saturation: 0.1"
    } else {
        title <- "Effect Saturation: 0.5"
    }
    
    plt <- plt + labs(x = "Signal-to-noise ratio", y = ylab, 
                      col = "Models", linetype = "Output type", shape = "Correlation adjustment",
                      title = title)
    return(plt)
}

process_outputs <- function(grid, results){
    df <- left_join(grid, results, by = "id") %>% 
        select(snr, sat, spar, s_rho, model, distr, adj, output, .metric, mean, std_err) %>% 
        unite(model, c(model, distr), sep = "_") %>%
        mutate(model = str_replace(model, "_NA", ""), 
               adj = replace_na(adj, "Not applicable"), 
               output = case_when(
                   output == "zscore" ~ "z-score",
                   output == "cdf" ~ "CDF", 
                   TRUE ~ "Not applicable"
               )) %>% 
        mutate(model = case_when(
            model == "cilr_mnorm" ~ "cILR Mixture Normal",
            model == "cilr_norm" ~ "cILR Normal",
            model == "gsva" ~ "GSVA", 
            model == "ssgsea" ~ "ssGSEA"
        )) %>%
        mutate(lower = mean - std_err, upper = mean + std_err) %>% 
        rename(c("Sparsity" = "spar",  "Correlation" = "s_rho"))
    return(df)
}


# Process results #### 

# classif ####
classif_grid <- readRDS(file = "analyses/simulations_prediction_classif/output/simulation_grid_classif.rds")
classif_results <- readRDS(file = "analyses/simulations_prediction_classif/output/sim_pred_classif.rds")
classif_df <- process_outputs(classif_grid, classif_results)
classif_names <- cross_df(list(
    metrics = c("accuracy", "roc_auc"),
    sat_value = c(0.1, 0.5)
))

classif_names <- classif_names %>% rowwise() %>% 
    mutate(plot = list(plot_func(metrics = metrics, sat_value = sat_value, df_main = classif_df)))

roc <- classif_names %>% filter(metrics == "roc_auc")
classif_auc <- roc$plot[[1]] / roc$plot[[2]] + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(plot.title = element_text(face  = "bold")) & ylim(c(0.7,1))
ggsave(classif_auc, filename = "figures/sim_pred_auc.png", dpi = 300, width = 8, height = 9)
file.copy(from = "figures/sim_pred_auc.png",
          to = "../teailr_manuscript/manuscript/figures/sim_pred_auc.png", overwrite = TRUE)

acc <- classif_names %>% filter(metrics == "accuracy")
classif_acc <- acc$plot[[1]] / acc$plot[[2]] + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(plot.title = element_text(face  = "bold")) & ylim(c(0.6,1))
ggsave(classif_auc, filename = "figures/sim_pred_acc.png", dpi = 300, width = 8, height = 9)
file.copy(from = "figures/sim_pred_acc.png",
          to = "../teailr_manuscript/manuscript/figures/sim_pred_acc.png", overwrite = TRUE)

# Regr ####
regr_grid <- readRDS(file = "analyses/simulations_prediction_regr/output/simulation_grid_regr.rds")
regr_results <- readRDS(file = "analyses/simulations_prediction_regr/output/sim_pred_regr.rds")
regr_df <- process_outputs(regr_grid, regr_results)
regr_names <- cross_df(list(
    metrics = c("rsq", "rmse"),
    sat_value = c(0.1, 0.5)
))

regr_names <- regr_names %>% rowwise() %>% 
    mutate(plot = list(plot_func(metrics = metrics, sat_value = sat_value, df_main = regr_df)))

rmse <- regr_names %>% filter(metrics == "rmse")
regr_rmse <- rmse$plot[[1]] / rmse$plot[[2]] + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(plot.title = element_text(face  = "bold"))
ggsave(regr_rmse, filename = "figures/sim_pred_rmse.png", dpi = 300, width = 8, height = 9)
file.copy(from = "figures/sim_pred_rmse.png",
          to = "../teailr_manuscript/manuscript/figures/sim_pred_rmse.png", overwrite = TRUE)


rsq <- regr_names %>% filter(metrics == "rsq")
regr_rsq <- rsq$plot[[1]] / rsq$plot[[2]] + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(plot.title = element_text(face  = "bold"))

ggsave(regr_rsq, filename = "figures/sim_pred_rsq.png", dpi = 300, width = 8, height = 9)
file.copy(from = "figures/sim_pred_rsq.png",
          to = "../teailr_manuscript/manuscript/figures/sim_pred_rsq.png", overwrite = TRUE)



