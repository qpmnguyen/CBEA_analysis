options(tidyverse.quiet = TRUE)
library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
library(CBEA)
library(phyloseq)
library(mia)
source("R/functions_data_diffab.R")

gingival_load <- function(){
    # readRDS(file = "data/hmp_supergingival_supragingival_16S.rds")
    data(hmp_gingival, package = "CBEA")
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
            output = c("zscore", "cdf"),
            size = c(20,50,100,150,200)
        ))
        addition <- cross_df(list(
            models = c("corncob", "deseq2"),
            size = c(20,50,100,150,200)
        ))
        addition_2 <- cross_df(list(
            models = "cbea", 
            distr = "norm", 
            output = "raw", 
            adj = FALSE, 
            size = c(20,50,100,150,200)
        ))
        settings <- full_join(settings, addition, by = c("models", "size"))
        settings <- full_join(settings, addition_2, by = c("models", "size", "adj", "distr", "output"))
        settings$id <- seq_len(nrow(settings))
    } else if (mode == "fdr"){
        settings <- cross_df(list(
            models = c("cbea"),
            distr = c("mnorm", "norm"),
            adj = c(TRUE, FALSE),
            output = c("zscore", "cdf")
        ))
        addition <- cross_df(list(
            models = c("corncob", "deseq2")
        )) 
        settings <- full_join(settings, addition, by = c("models"))
        settings <- bind_rows(settings, data.frame(models = "cbea", distr = "norm", adj = TRUE, output = "raw"))
        settings$id <- seq_len(nrow(settings))
    }
    return(settings)
}


# FDR ANALYSIS #### 
fdr_analysis <- tar_map(unlist = FALSE, values = get_settings(mode = "fdr"), 
    tar_target(index_batch, seq_len(1)),
    tar_target(index_rep, seq_len(1)),
    tar_target(input_data, gingival_load()),
    tar_target(rand_seq, {
        purrr::map(index_rep, ~{
            obj <- input_data$data
            colData(obj)$group <- factor(rbinom(n = nrow(colData(obj)), size = 1, prob = 0.5))
            return(obj)
        })
    }, pattern = map(index_batch)),
    tar_target(diff_analysis, {
        purrr::map(rand_seq, ~diff_ab(obj = .x, eval = "fdr", method = models, 
                                      thresh = 0.05, return = "sig", 
                                      distr = distr, adj = adj, output = output))
    }, pattern = map(rand_seq)), 
    tar_target(eval_diff, {
        purr::map(diff_analysis, ~{
            tibble(
                models = models, 
                distr = distr, 
                adj = adj,
                output = output, 
                res = sum(.x == 1)/length(.x)
            )
        })
    }, pattern = map(diff_analysis))
)

fdr_summary <- tar_combine(combine_fdr, fdr_analysis[[6]], command = dplyr::bind_rows(!!!.x))

fdr_save <- tarchetypes::tar_rds(save_fdr, saveRDS(combine_fdr, 
                                                   "output/data_diffab_fdr.rds"))


# RANDOM SET ANALYSIS ####
rset_analysis <- tar_map(unlist = FALSE, values = get_settings("rset"),
               tar_target(input_data, {
                   obj <- gingival_load()$data
                   colData(obj)$group <- factor(ifelse(colData(obj)$HMP_BODY_SUBSITE == "Supragingival Plaque", 1,0))
               }),
               tar_target(rand_set, {
                   purrr::map(index_rep, ~get_rand_sets(input_data, size = size, n_sets = 100))
               }, pattern = map(index_batch)),
               tar_target(enrich_test, {
                   purrr::map(rand_set, ~diff_ab(obj = input_data,
                                                 eval = "rset", 
                                                 return = "sig",
                                                 abund_values = "16SrRNA",
                                                 sets = .x, method = models,
                                                 distr = distr, adj = adj, output = output))
               }, pattern = map(rand_set)),
               tar_target(enrich_eval, {
                   purrr::map_dfr(enrich_test, ~{
                       tibble(
                           models = models,
                           distr = distr,
                           adj = adj,
                           size = size,
                           res = sum(.x == 1)/length(.x)
                       )
                   })
               }, pattern = map(enrich_test))
)

rset_summary <- tar_combine(combine_rset, rset_analysis[[4]], command = dplyr::bind_rows(!!!.x))

rset_save <- tarchetypes::tar_rds(save_rset, saveRDS(combine_rset, 
                                                   "output/data_diffab_rset.rds"))


list(
    fdr_analysis, fdr_summary, fdr_save,
    rset_analysis, rset_summary, rset_save
)


