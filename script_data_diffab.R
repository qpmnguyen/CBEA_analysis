options(tidyverse.quiet = TRUE)
library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
library(phyloseq)
source("R/functions_data_diffab.R")

gingival_load <- function(){
    readRDS(file = "data/hmp_supergingival_supragingival_16S.rds")
}

# DEFINE SETTINGS ####
# if not auc, then only care about significant outcomes 
get_settings <- function(mode){
    mode <- match.arg(mode, choices = c("fdr","rset"))
    if (mode == "rset"){
        settings <- cross_df(list(
            models = c("cbea"),
            distr = c("mnorm", "norm"),
            adj = c(TRUE, FALSE),
            output = c("zscore", "cdf", "raw"),
            size = c(20,50,100,150,200)
        ))
        addition <- cross_df(list(
            models = c("corncob", "deseq2"),
            size = c(20,50,100,150,200)
        ))
        settings <- full_join(settings, addition, by = c("models", "size"))
        settings$id <- seq_len(nrow(settings))
    } else if (mode == "fdr"){
        settings <- cross_df(list(
            models = c("cbea"),
            distr = c("mnorm", "norm"),
            adj = c(TRUE, FALSE),
            output = c("zscore", "cdf", "raw")
        ))
        addition <- cross_df(list(
            models = c("corncob", "deseq2")
        ))
        settings <- full_join(settings, addition, by = c("models"))
        settings$id <- seq_len(nrow(settings))
    }
    return(settings)
}



fdr_analysis <- tar_map(unlist = FALSE, values = get_settings("fdr"), 
    tar_target(index_batch, seq_len(10)),
    tar_target(index_rep, seq_len(1)),
    tar_target(input_data, {
        gingival_process(gingival_load())
    }),
    tar_target(rand_label, {
        purrr::map(index_rep, ~{
            sample_data(input_data$physeq)[,"group"] <- rbinom(n = nsamples(input_data$physeq), 
                                                               size = 1, prob = 0.5)
            return(input_data$physeq)
        })
    }, pattern = map(index_batch)),
    tar_target(diff_analysis, {
        purrr::map(rand_label, ~{
            
        })
    }, pattern = map(rand_label)), 
    tar_target(eval_diff, {
        
    })
)


rset_analysis <- tar_map(unlist = FALSE, values = get_settings("sig"), 
               tar_target(index_batch, seq_len(10)),
               tar_target(index_rep, seq_len(1)),
               tar_target(input_data, {gingival_load()$physeq}),
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

