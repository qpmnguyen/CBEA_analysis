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
# define some settings to map across 

# these are different settings for cilr method  
cilr_settings <- cross_df(list(
    models = c("cilr"),
    distr = c("mnorm", "norm"),
    adj = c(TRUE, FALSE),
    output = c("zscore", "cdf")
))

# these are different single sample enrichment models  
auc_models <- tibble(
    models = c("gsva", "ssgsea")
)

# these are different models to test for power  
pwr_models <- tibble(
    models = c("wc")
)


data_enrich <- tar_rds(data_enrich, {
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
     }, batches = 10, reps = 10)
 })

pwr_cilr <- tar_map(unlist = FALSE, values = cilr_settings,{
    tar_rep(pwr_cilr, {
        X <- data_enrich$X
        idx <- sample(1:nrow(X), size = nrow(X), replace = F)
        X_boot <- X[idx,]
        label_boot <- data_enrich$label[idx]
        pwr <- enrichment_analysis(X = X_boot, A = data_enrich$A, method = models, label = label_boot, 
                                   distr = distr, adj = adj, output = output)
        data.frame(pwr = pwr, models = models, distr = distr, adj = adj, output = output)
    }, batches = 10, reps = 10)
})


auc_other <- tar_map(unlist = FALSE, values = auc_models, {
    tar_rep(auc_models, {
        X <- data_enrich$X
        idx <- sample(1:nrow(X), size = nrow(X), replace = F)
        X_boot <- X[idx, ]
        label_boot <- data_enrich$label[idx]
        auc <- enrichment_analysis(X = X_boot, A = data_enrich$A, method = models, label = label_boot)
        data.frame(auc = auc, models = models)
    }, batches = 10, reps = 10)
})

auc_combined <- tar_combine(auc, auc_cilr[[2]], auc_other[[2]], command = dplyr::bind_rows(!!!.x))

save_auc <- tarchetypes::tar_rds(auc_save, saveRDS(auc, "output/auc_comparison.rds"))


list(data_enrich, auc_cilr, auc_other, auc_combined)
