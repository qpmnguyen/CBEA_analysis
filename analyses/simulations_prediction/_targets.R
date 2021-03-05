library(targets)
library(tarchetypes)

source("functions/pred_functions.R")

set.seed(1020)

tar_option_set(error = "workspace")

sim_grid <- 

sim_grid <- tar_target(auc_test_grid, {
    eval_settings <- cross_df(list(
        model = "cilr",
        distr = c("mnorm", "norm"),
        adj = c(TRUE, FALSE),
        output = c("zscore", "cdf")
    ))
    other <- tibble(model = c("ssgsea", "gsva"))
    eval_settings <- full_join(eval_settings, other, by = "model")
    eval_settings
})