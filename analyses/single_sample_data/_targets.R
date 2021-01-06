library(targets)
library(here)
library(tarchetypes)
library(tidyverse)
library(future)
# library(future.batchtools)
source(here("functions","utils.R"))
source(here("functions", "cilr.R"))
source(here("functions", "enrichment.R"))

plan(multisession)
options(tidyverse.quite = TRUE)
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


data <- tar_rds(data, {
    readRDS(here("data", "hmp_supergingival_supragingival_16S.rds")) %>% enrichment_processing()
})

auc_cilr <- tar_map(unlist = FALSE, values = cilr_settings, {
    tar_rep(auc_cilr, {
        X <- data$X
        idx <- sample(1:nrow(X), size = nrow(X), replace = F)
        X_boot <- X[idx,]
        label_boot <- data$label[idx]
        data.frame(auc = enrichment_analysis(X = X_boot, A = data$A, method = models, label = label_boot, 
                                             distr = distr, adj = adj, output = output))
    }, batches = 2, reps = 2)
})
auc_other <- tar_map(unlist = FALSE, values = auc_models, {
    tar_rep(auc_models, {
        X <- data$X
        idx <- sample(1:nrow(X), size = nrow(X), replace = F)
        X_boot <- X[idx, ]
        label_boot <- data$label[idx]
        data.frame(auc = enrichment_analysis(X = X_boot, A = data$A, method = models, label = label_boot))
    }, batches = 2, reps = 2)
})

combined <- tar_combine(auc, auc_other[[1]], auc_cilr[[1]], command = dplyr::bind_rows(!!!.x))

list(data, auc_cilr, auc_other)
