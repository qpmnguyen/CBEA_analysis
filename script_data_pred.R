library(targets)
library(tarchetypes)
library(future)
library(tidyverse)
library(mia)
library(phyloseq)
library(BiocSet)
source("R/functions_data_pred.R")
source("R/utils.R")

plan(multisession)

# DEFINE SETTINGS ####
ibd_load <- function(type){
    type <- match.arg(type, c("16s", "wgs"))
    if (type == "16s"){
        data <- readRDS("data/ackerman_ibd_16S.rds")
        label <- if_else(sample_data(data)$diagnosis %in% c("CD", "UC"), "case", "control")
        sample_data(data)$outcome <- factor(label, c("case", "control"))
    } else {
        data <- readRDS("data/nielsen_ibd_wgs.rds")
        label <- if_else(sample_data(data)$disease == "IBD", "case", "control")
        sample_data(data)$outcome <- factor(label, c("case", "control"))
        sample_data(data)$sample_name <- sample_data(data)$subjectID
    }
    #data <- mia::makeTreeSummarizedExperimentFromPhyloseq(data)
    return(data)
}


# if not auc, then only care about significant outcomes 
get_settings <- function(){
    settings <- cross_df(list(
        models = c("cbea"),
        data = c("16s", "wgs"),
        distr = c("mnorm", "norm"),
        adj = c(TRUE, FALSE),
        output = c("zscore", "cdf")
    ))
    addition <- cross_df(list(
        data = c("16s", "wgs"),
        models = c("ssgsea", "gsva", "clr")
    ))
    addition_2 <- cross_df(list(
        data = c("16s", "wgs"), 
        models = "cbea", distr = "norm", adj = FALSE, output = "raw"))
    settings <- full_join(settings, addition)
    settings <- full_join(settings, addition_2)
    settings$id <- seq_len(nrow(settings))
    return(settings)
}

# FIT MODELS #### 
pred <- tar_map(unlist = FALSE, values = get_settings(), 
               tar_target(input_data, ibd_load(data)),
               tar_target(proc_df, {
                   generate_aggregation(input_data, method = models, output = output, distr = distr, 
                                        adj = adj)
               }),
               tar_target(fit, {
                   results <- fit_and_eval(proc_df$aug_df, nfolds = 10, task = "classification")
                   results
               }), 
               tar_target(eval, {
                   res <- fit %>% dplyr::filter(.metric == "roc_auc")
                   tibble(
                       data = data,
                       models = models, 
                       distr = distr, 
                       adj = adj, 
                       output = output,
                       estimate = res %>% dplyr::pull(mean),
                       std_err = res %>% dplyr::pull(std_err)
                   )
               })
)



pred_summary <- tar_combine(combine_pred, pred[[4]], command = dplyr::bind_rows(!!!.x))

pred_save <- tarchetypes::tar_rds(save_pred, saveRDS(combine_pred, 
                                                   "output/data_pred.rds"))

list(pred, pred_summary, pred_save)



