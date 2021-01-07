options(tidyverse.quiet = TRUE)
library(targets)
library(here)
library(tarchetypes)
library(tidyverse)
library(future)
# library(future.batchtools)
source(here("functions","utils.R"))
source(here("functions", "cilr.R"))
source(here("functions", "enrichment.R"))

#plan(multisession)
#plan(batchtools_torque, template = "batchtools.torque.tmpl")

set.seed(2105)
cilr_settings <- cross_df(list(
    models = c("cilr"),
    distr = c("mnorm", "norm"),
    adj = c(TRUE, FALSE),
    output = c("zscore", "cdf")
))

auc_models <- tibble(
    models = c("gsva", "ssgsea")
)


data <- tar_rds(data_enrich, {
    readRDS(here("data", "hmp_supergingival_supragingival_16S.rds")) %>% enrichment_processing()
})


auc_cilr <- tar_map(unlist = FALSE, values = cilr_settings, {
     tar_rep(auc_cilr, {
         X <- data_enrich$X
         idx <- sample(1:nrow(X), size = nrow(X), replace = F)
         X_boot <- X[idx,]
         label_boot <- data_enrich$label[idx]
         auc <- enrichment_analysis(X = X_boot, A = data_enrich$A, method = models, label = label_boot, 
                                              distr = distr, adj = adj, output = output)
         data.frame(auc = auc, models = models, distr = distr, adj = adj, output = output)
     }, batches = 1, reps = 1)
 })

auc_other <- tar_map(unlist = FALSE, values = auc_models, {
    tar_rep(auc_models, {
        X <- data_enrich$X
        idx <- sample(1:nrow(X), size = nrow(X), replace = F)
        X_boot <- X[idx, ]
        label_boot <- data_enrich$label[idx]
        auc <- enrichment_analysis(X = X_boot, A = data_enrich$A, method = models, label = label_boot)
        data.frame(auc = auc, models = models)
    }, batches = 1, reps = 1)
})

auc_combined <- tar_combine(auc, auc_cilr[[2]], auc_other[[2]], command = dplyr::bind_rows(!!!.x))

list(data, auc_cilr, auc_other, auc_combined)
