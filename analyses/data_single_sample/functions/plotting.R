library(tidyverse)
library(ggsci)
library(showtext)
font_add_google("Lora")
showtext_auto()


generate_plots <- function(df, metric){
    df <- df %>% mutate(distr = ifelse(is.na(distr), "None", distr), 
                        adj = ifelse(is.na(adj), FALSE, adj))
    showtext_begin()
    if (metric == "fdr"){
        plt <- ggplot(df, aes(x = models, y = fdr, color = adj, shape = distr)) + 
            geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 1)) + 
            ggsci::scale_color_npg(labels = c("Adjusted", "Unadjusted")) + 
            theme_bw() + theme(text = element_text(family = "Roboto", size = 20)) +
            scale_shape_discrete(labels = c("Mixture Normal", "None", "Normal")) + 
            geom_hline(yintercept = 0.05, col = "red") + 
            labs(x = "Models", y = "Type I error", color = "Correlation adjustment", shape = "Distribution")
    } else if (metric == "pwr"){
        plt <- ggplot(df, aes(x = models, y = pwr, color = adj, shape = distr)) + 
            geom_pointrange(aes(ymax = upper, ymin = lower), position = position_dodge(width = 1)) + 
            ggsci::scale_color_npg(labels = c("Adjusted", "Unadjusted")) + 
            theme_bw() + theme(text = element_text(family = "Roboto", size = 20)) +
            scale_shape_discrete(labels = c("Mixture Normal", "None", "Normal")) + 
            geom_hline(yintercept = 0.05, col = "red") + 
            labs(x = "Models", y = "Power", color = "Correlation adjustment", shape = "Distribution")
    } else if (metric == "auc"){
        df <- df %>% group_by(models, fdr, adj, distr) %>% summarise()
    }
    return(plt)
    showtext_en()
}