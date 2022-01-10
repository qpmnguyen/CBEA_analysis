options(tidyverse.quiet = TRUE)
library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
library(phyloseq)
# library(future.batchtools)
source("R/functions_data_ss.R")
tar_option_set(workspace_on_error = TRUE)

plan(multisession)
#plan(batchtools_torque, template = "batchtools.torque.tmpl")

set.seed(1020)

# DEFINE SETTINGS ####
# if not auc, then only care about significant outcomes 
get_settings <- function(mode){
    mode <- match.arg(mode, choices = c("sig","pheno", "auc"))
    if (mode == "sig"){
        settings <- cross_df(list(
            models = c("cbea"),
            distr = c("mnorm", "norm"),
            adj = c(TRUE, FALSE), 
            size = c(20,50,100,150,200)
        ))
        addition <- cross_df(list(
            models = c("wilcoxon"),
            size = c(20,50,100,150,200)
        ))
        settings <- full_join(settings, addition, by = c("models", "size"))
        settings$id <- seq_len(nrow(settings))
    } else if (mode == "auc"){
        settings <- cross_df(list(
            models = c("cbea"),
            distr = c("mnorm", "norm"),
            adj = c(TRUE, FALSE)
        ))
        addition <- cross_df(list(
            models = c("ssgsea", "gsva", "wilcoxon")
        
        ))
        settings <- full_join(settings, addition, by = c("models"))
        settings$id <- seq_len(nrow(settings))
    } else if (mode == "pheno"){
        settings <- cross_df(list(
            models = c("cbea"),
            distr = c("mnorm", "norm"),
            adj = c(TRUE, FALSE) 
        ))
        addition <- cross_df(list(
            models = c("wilcoxon")
        ))
        settings <- full_join(settings, addition, by = c("models"))
        settings$id <- seq_len(nrow(settings))
        
    }
    return(settings)
}

# ANALYSIS FOR RANDOM SET FDR ####
# data 2 and 3  
ibd_load <- function(type){
    type <- match.arg(type, c("16s", "wgs"))
    if (type == "16s"){
        return(readRDS("data/ackerman_ibd_16S.rds"))
    } else (
        return(readRDS("data/nielsen_ibd_wgs.rds"))
    )
}

#' @title Function that loads the data 
gingival_load <- function(){
    readRDS(file = "data/hmp_supergingival_supragingival_16S.rds")
}

fdr <- tar_map(unlist = FALSE, values = get_settings("sig"), 
               tar_target(index_batch, seq_len(10)),
               tar_target(index_rep, seq_len(100)),
               tar_rds(input_data, gingival_load()),
               tar_target(rand_set, {
                   purrr::map(index_rep, ~get_rand_sets(input_data, size = size, n_sets = 1))
               }, pattern = map(index_batch)),
               tar_target(enrich_test, {
                   purrr::map(rand_set, ~enrichment_analysis(input_data, 
                                                             set = .x, method = models, 
                                                             metric = "fdr", distr = distr, 
                                                             adj = adj, output = "sig"))
               }, pattern = map(rand_set)),
               tar_target(enrich_eval, {
                   purrr::map_dfr(enrich_test, ~{
                       tibble(
                           type = type,
                           models = models,
                           distr = distr,
                           adj = adj,
                           size = size,
                           res = sum(.x[,2])/nrow(.x)
                       )
                   })
               }, pattern = map(enrich_test))
)



fdr_summary <- tar_combine(combine_fdr, fdr[[6]], command = dplyr::bind_rows(!!!.x))

fdr_save <- tarchetypes::tar_rds(save_fdr, saveRDS(combine_fdr, 
                                                   "output/fdr_ss_randset.rds"))

# ACCURACY ANALYSES - PHENOTYPE RELEVANCE ####

phenotype_sig <- tar_map(unlist = FALSE, values = get_settings("pheno"), 
    tar_target(sig_data, gingival_load()),
    tar_target(proc_sig_data, gingival_processing(data = sig_data)),
    tar_target(sig_scores, enrichment_analysis(proc_sig_data$physeq, 
                                             set = proc_sig_data$sets, 
                                             method = models, 
                                             metric = "fdr", distr = distr, 
                                             adj = adj, output = "sig")),
    tar_target(sig_eval, gingival_evaluate(physeq = proc_sig_data$physeq, 
                                         results = sig_scores)),
    tar_target(sig_metric, {
        result <- calculate_statistic(eval = "pwr", 
                                      pred = sig_eval$Aerobic, 
                                      true = sig_eval$label)
        as_tibble(data.frame(
            models = models, 
            distr = distr, 
            adj = adj,
            result
        ))
    })
)



phenotype_sig_summary <- tar_combine(combine_pheno_sig, phenotype_sig[[5]], 
                                     command = dplyr::bind_rows(!!!.x))

phenotype_sig_save <- tarchetypes::tar_rds(save_pheno_sig, saveRDS(combine_pheno_sig, 
                                                   "output/pwr_ss_pheno.rds"))

# AUC ANALYSES - PHENOTYPE RELEVANCE ####  
phenotype_auc <- tar_map(unlist = FALSE, values = get_settings("auc"),
     tar_target(auc_data, gingival_load()),
     tar_target(proc_auc_data, gingival_processing(data = auc_data)),
     tar_target(auc_scores, enrichment_analysis(proc_auc_data$physeq, 
                                              set = proc_auc_data$sets, 
                                              method = models, 
                                              metric = "auc", distr = distr, 
                                              adj = adj, output = "sig")),
     tar_target(auc_eval, gingival_evaluate(physeq = proc_auc_data$physeq, 
                                          results = auc_scores)),
     tar_target(auc_metric, {
         result <- calculate_statistic(eval = "auc", 
                                       pred = auc_eval$Aerobic, 
                                       true = auc_eval$label)
         as_tibble(data.frame(
             models = models, 
             distr = distr, 
             adj = adj,
             result
         ))
     })
)

phenotype_auc_summary <- tar_combine(combine_pheno_auc, phenotype_auc[[5]], 
                                     command = dplyr::bind_rows(!!!.x))

pheno_auc_save <- tarchetypes::tar_rds(save_pheno_auc, saveRDS(combine_pheno_auc, 
                                                               "output/auc_ss_pheno.rds"))



list(
    fdr, fdr_summary, fdr_save,
    phenotype_sig, phenotype_sig_summary, phenotype_sig_save,
    phenotype_auc, phenotype_auc_summary, pheno_auc_save
)






