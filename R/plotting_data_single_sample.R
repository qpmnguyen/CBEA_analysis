library(tidyverse)
library(ggsci)
library(patchwork)
library(glue)

if(Sys.info()["sysname"] == "Darwin"){
    save_dir <- "../cilr_manuscript/figures"
} else {
    save_dir <- "../teailr_manuscript/manuscript/figures"
}

# Processing fdr and pwr ####  
load_and_process <- function(type = c("fdr", "pwr")){
    type <- match.arg(type)
    real_data <- readRDS(file = glue("analyses/data_single_sample/output/{type}_comparison.rds",
                                     type = type))
    sim_grid <- readRDS(file = glue( "analyses/simulations_single_sample_{type}/output/simulation_grid_{type}.rds", 
                                     type = type))
    sim_results <- readRDS(file = glue( "analyses/simulations_single_sample_{type}/output/sim_ss_{type}.rds", 
                                        type = type))
    full_data <- full_join(sim_grid, sim_results, by = "id")
    full_data <- full_data %>% unite(model, c("model", "distr")) %>% 
        mutate(model = case_when(
            model == "cilr_mnorm" ~ "cILR Mixture Normal",
            model == "cilr_norm" ~ "cILR Normal", 
            model == "wilcox_NA" ~ "Wilcoxon Rank Sum"
        )) %>% 
        mutate(adj = replace_na(adj, "Not applicable")) %>% 
        rename("Correlation" = "s_rho", "Effect Size" = "eff_size", "Set Size" = "n_inflate")
    
    real_data <- real_data %>% unite(models, c("models", "distr")) %>% 
        mutate(model = case_when(
            models == "cilr_mnorm" ~ "cILR Mixture Normal",
            models == "cilr_norm" ~ "cILR Normal", 
            models == "wilcox_NA" ~ "Wilcoxon Rank Sum",
        )) %>% 
        mutate(adj = replace_na(adj, "Not applicable")) 
    
    sim_plots <- ggplot(full_data, aes(x = spar, y = est, col = model, shape = adj)) + 
        geom_line() + 
        geom_point() + 
        geom_linerange(aes(ymax = upper, ymin = lower), show.legend = FALSE) + 
        theme_bw() + scale_color_d3() + 
        guides(linetype = guide_legend(override.aes = list(shape = NA)))
    
    if (type == "fdr"){
        sim_plots <- sim_plots + facet_grid(`Correlation` ~ `Set Size`, labeller = label_both) + 
            labs(x = "Sparsity", y = "Type I error", col = "Models", shape = "Correlation adjusted") +
            geom_hline(yintercept = 0.05, col = "red") 
    } else if (type == "pwr"){
        sim_plots <- sim_plots + facet_grid(`Correlation` ~ `Effect Size`, labeller = label_both) + 
            labs(x = "Sparsity", y = "Power", col = "Models", shape = "Correlation adjusted") + 
            geom_hline(yintercept = 0.8, col = "red") 
    }
    
    data_plot <- ggplot(real_data, aes(x = model, y = est, col = model, shape = adj)) + 
        geom_line() + 
        geom_point(position = position_dodge(width = 0.5)) + 
        geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), 
                      show.legend = FALSE) + 
        scale_color_d3() +
        guides(linetype = guide_legend(override.aes = list(shape = NA))) + theme_bw()
    
    if (type == "fdr"){
        data_plot <- data_plot + 
            labs(x = "Models", y = "False Positive Rate", color = "Models", shape = "Correlation adjusted") + 
            geom_hline(yintercept = 0.05, col = "red")
    } else if (type == "pwr"){
        data_plot <- data_plot + 
            labs(x = "Models", y = "True Positive Rate", color = "Models", shape = "Correlation adjusted") + 
            geom_hline(yintercept = 0.8, col = "red")
    }
    data_plot <- data_plot + coord_flip()
    return(list(data_plot = data_plot, sim_plot = sim_plots))
}

fdr <- load_and_process("fdr")
pwr <- load_and_process("pwr")

layout <- "
AA#
AAC
AAC
BBD
BBD
BB#
"
comb_plot <- fdr$sim_plot + pwr$sim_plot + fdr$data_plot + pwr$data_plot + 
    plot_layout(guides = "collect", design = layout, widths = c(2,1)) + 
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin())
saveRDS(comb_plot, file = "figures/sim_data_ss_hypo.rds")
ggsave(comb_plot, filename = "figures/sim_data_ss_hypo.png", dpi = 300, width = 10, height = 9)
ggsave(comb_plot, filename = "figures/sim_data_ss_hypo.eps", dpi = 300, width = 10, height = 9)

# Plotting AUC data ####
auc_data <- readRDS(file = "analyses/data_single_sample/output/auc_comparison.rds")
sim_auc_grid <- readRDS(file = "analyses/simulations_single_sample_auc/output/simulation_grid_auc.rds")
sim_auc_results <- readRDS(file = "analyses/simulations_single_sample_auc/output/sim_ss_auc.rds")
full_data <- full_join(sim_auc_grid, sim_auc_results, by = "id")
full_data <- full_data %>% unite(model, c("model", "distr")) %>% 
    mutate(model = case_when(
        model == "cilr_mnorm" ~ "cILR Mixture Normal",
        model == "cilr_norm" ~ "cILR Normal", 
        model == "gsva_NA" ~ "GSVA", 
        model == "ssgsea_NA" ~ "ssGSEA",
        model == "wilcox_NA" ~ "Wilcoxon W Statistic"
    )) %>% 
    mutate(output = case_when(
        output == "cdf" ~ "CDF",
        output == "zscore" ~ "z-scores",
        is.na(output) ~ "Not applicable"
    )) %>% 
    mutate(adj = replace_na(adj, "Not applicable")) %>% 
    rename("Correlation" = "s_rho", "Effect Size" = "eff_size")

auc_data <- auc_data %>% unite(models, c("models", "distr")) %>% 
    mutate(model = case_when(
        models == "cilr_mnorm" ~ "cILR Mixture Normal",
        models == "cilr_norm" ~ "cILR Normal", 
        models == "gsva_NA" ~ "GSVA", 
        models == "ssgsea_NA" ~ "ssGSEA",
        models == "wilcoxon_NA" ~ "Wilcoxon W Statistic"
    )) %>% 
    mutate(output = case_when(
        output == "cdf" ~ "CDF",
        output == "zscore" ~ "z-scores",
        is.na(output) ~ "Not applicable"
    )) %>% 
    mutate(adj = replace_na(adj, "Not applicable")) 

sim_plot <- ggplot(full_data, aes(x = spar, y = est, col = model, shape = output)) +
    geom_line(aes(linetype = adj)) + 
    geom_point() + 
    geom_linerange(aes(ymax = upper, ymin = lower), show.legend = FALSE) + 
    facet_grid(`Correlation` ~ `Effect Size`, labeller = label_both) + theme_bw() + scale_color_d3() + 
    labs(x = "Sparsity", y = "AUC", col = "Models", linetype = "Correlation adjusted", shape = "Output type") +
    guides(linetype = guide_legend(override.aes = list(shape = NA)))

data_plot <- ggplot(auc_data, aes(x = model, y = est, col = model, linetype = adj, shape = output)) + 
    geom_line() + 
    geom_point(position = position_dodge(width = 0.5)) + 
    geom_linerange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), 
                  show.legend = FALSE) + 
    labs(x = "Models", y = "AUC", color = "Models", 
         shape = "Output type", linetype = "Correlation adjusted") +
    scale_color_d3() +
    guides(linetype = guide_legend(override.aes = list(shape = NA))) + theme_bw() + 
    coord_flip() + 
    ylim(c(0.6,0.8))

layout <- "
AAB
AAB
"

comb_plot <- sim_plot + data_plot + plot_layout(guides = "collect", design = layout) + 
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin())
#saveRDS(comb_plot, file = "figures/sim_data_ss_auc.rds")
ggsave(comb_plot, filename = "figures/sim_data_ss_auc.png", dpi = 300, width = 10, height = 6)
ggsave(comb_plot, filename = "figures/sim_data_ss_auc.eps", dpi = 300, width = 10, height = 6)



file.copy(from = Sys.glob("figures/*.png"), 
          to = glue("{save_dir}", save_dir = save_dir), overwrite = TRUE)

file.copy(from = Sys.glob("figures/*.eps"), 
          to = glue("{save_dir}", save_dir = save_dir), recursive = TRUE, overwrite = TRUE)


### FDR PLOTS NEW #### 
df <- readRDS(file = "output/fdr_ss_randset.rds")
df <- df %>% mutate(adj = if_else(adj, "Yes", "No")) %>% 
    mutate(adj = if_else(is.na(adj), "Not applicable", adj), 
           models = if_else(is.na(models), "Not applicable", models)) %>%
    unite("models", models:distr) %>% 
    mutate(models = if_else(str_detect(models, "_NA"), 
                            str_remove_all(models, "_NA"), models)) %>% 
    group_by(type, models, adj, size) %>% 
    summarise(estimate = mean(res), 
              se = sd(res)/sqrt(500)) %>% 
    mutate(upper = estimate + se, lower = estimate - se) %>% ungroup()



ggplot(df %>% filter(type == "16s"), aes(x = models, y = estimate, col = models, shape = adj)) +
    geom_point(position = position_dodge(width = 1)) +
    scale_color_d3() +
    geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 1)) + 
    facet_grid(. ~ size) + theme_bw()
