options(tidyverse.quiet = TRUE)
library(targets)
library(tarchetypes)
library(tidyverse)
library(future)
library(phyloseq)
# library(future.batchtools)
source("functions/enrichment.R")
tar_option_set(error = "workspace")

plan(multisession)
#plan(batchtools_torque, template = "batchtools.torque.tmpl")

set.seed(1020)
# define some settings to map across 
# these are different settings for cilr method  
cilr_settings <- cross_df(list(
    models = c("cilr"),
    distr = c("mnorm", "norm"),
    adj = c(TRUE, FALSE),
    output = c("zscore", "cdf")
))

# if not auc, then only care about significant outcomes 
cilr_settings_sig <- cross_df(list(
    models = c("cilr"),
    distr = c("mnorm", "norm"),
    adj = c(TRUE, FALSE)
))

# these are different single sample enrichment models  
auc_models <- tibble(
    models = c("gsva", "ssgsea", "wilcoxon")
)

# data enrichment 
data_auc <- tar_rds(data_enrich, {
    readRDS("../../data/hmp_supergingival_supragingival_16S.rds") %>% enrichment_processing() -> data_enrich
})

# all cilr models under different evaluations  
auc_cilr_job <- tar_map(unlist = FALSE, values = cilr_settings, {
     tar_target(auc_cilr, {
         X <- data_enrich$X
         idx <- sample(1:nrow(X), size = nrow(X), replace = F)
         X_boot <- X[idx,]
         label_boot <- data_enrich$label[idx]
         auc <- enrichment_analysis(X = X_boot, A = data_enrich$A, method = models, metric = "auc", label = label_boot, 
                                              distr = distr, adj = adj, output = output)
         data.frame(est = auc$est, upper = auc$upper, lower = auc$lower, 
                    models = models, distr = distr, adj = adj, output = output)
     })
 })



fdr_cilr_job <- tar_map(unlist = FALSE, values = cilr_settings_sig,{
    tar_target(fdr_cilr, {
        X <- data_enrich$X
        idx <- sample(1:nrow(X), size = nrow(X), replace = F)
        X_boot <- X[idx,]
        label_boot <- data_enrich$label[idx]
        fdr <- enrichment_analysis(X = X_boot, A = data_enrich$A, method = models, label = label_boot, 
                                   distr = distr, adj = adj, metric = "fdr")
        data.frame(est = fdr$est, upper = fdr$upper, lower = fdr$lower, models = models, distr = distr, adj = adj, output = "sig")
    })
})

pwr_cilr_job <- tar_map(unlist = FALSE, values = cilr_settings_sig,{
    tar_target(pwr_cilr, {
        X <- data_enrich$X
        idx <- sample(1:nrow(X), size = nrow(X), replace = F)
        X_boot <- X[idx,]
        label_boot <- data_enrich$label[idx]
        pwr <- enrichment_analysis(X = X_boot, A = data_enrich$A, method = models, label = label_boot, 
                                   distr = distr, adj = adj, metric = "pwr")
        data.frame(est = pwr$est, upper = pwr$upper, lower = pwr$lower, models = models, distr = distr, adj = adj, output = "sig")
    })
})

# Models for comparison  
auc_other_job <- tar_map(unlist = FALSE, values = auc_models, {
    tar_target(auc_models, {
        X <- data_enrich$X
        idx <- sample(1:nrow(X), size = nrow(X), replace = F)
        X_boot <- X[idx, ]
        label_boot <- data_enrich$label[idx]
        auc <- enrichment_analysis(X = X_boot, A = data_enrich$A, method = models, label = label_boot, metric = "auc")
        data.frame(est = auc$est, upper = auc$upper, lower = auc$lower, models = models)
    })
})

pwr_other_job <- tar_target(pwr_models, {
    X <- data_enrich$X
    idx <- sample(1:nrow(X), size = nrow(X), replace = F)
    X_boot <- X[idx,]
    label_boot <- data_enrich$label[idx]
    pwr <- enrichment_analysis(X = X_boot, A = data_enrich$A, method = "wilcoxon", 
                               label = label_boot, metric = "pwr", output = "sig")
    data.frame(est = pwr$est, upper = pwr$upper, lower = pwr$lower, models = "wilcox")
}, error = "workspace")

fdr_other_job <- tar_target(fdr_models, {
    X <- data_enrich$X
    idx <- sample(1:nrow(X), size = nrow(X), replace = F)
    X_boot <- X[idx,]
    label_boot <- data_enrich$label[idx]
    fdr <- enrichment_analysis(X = X_boot, A = data_enrich$A, method = "wilcoxon", 
                               label = label_boot, metric = "fdr", output = "sig")
    data.frame(est = fdr$est, upper = fdr$upper, lower = fdr$lower, models = "wilcox")
}, error = "workspace")


auc_combined <- tar_combine(auc, auc_cilr_job, auc_other_job, command = dplyr::bind_rows(!!!.x))
pwr_combined <- tar_combine(pwr, pwr_cilr_job, pwr_other_job, command = dplyr::bind_rows(!!!.x))
fdr_combined <- tar_combine(fdr, fdr_cilr_job, fdr_other_job, command = dplyr::bind_rows(!!!.x))
save_auc <- tarchetypes::tar_rds(auc_save, saveRDS(auc, "output/auc_comparison.rds"))
save_fdr <- tarchetypes::tar_rds(fdr_save, saveRDS(fdr, "output/fdr_comparison.rds"))
save_pwr <- tarchetypes::tar_rds(pwr_save, saveRDS(pwr, "output/pwr_comparison.rds"))

#list(data_enrich, auc_cilr, auc_other, auc_combined)
#list(data_enrich, fdr_cilr, fdr_other, fdr_combined)
#list(data_enrich, pwr_cilr, pwr_other, pwr_combined)
list(data_auc, 
     auc_cilr_job, auc_other_job, auc_combined, 
     pwr_cilr_job, pwr_other_job, pwr_combined, 
     fdr_cilr_job, fdr_other_job, fdr_combined,
     save_auc, save_fdr, save_pwr)

