library(tidyverse)
library(furrr)
library(optparse)

source("../R/cilr.R")
source("../R/utils.R")
source("../R/diff_ab_functions.R")

option_list <- list(
    make_option("--ncores", type="integer", help="Number of cores going to be used")
)

opt <- parse_args(OptionParser(option_list=option_list))

ncores <- opt$ncores
dir <- "diff_ab_sim"

sim <- readRDS(file = glue("{dir}/parameters.rds", dir = dir))

eval_settings <- cross_df(list(
    methods = c("cilr_welch", "cilr_wilcox", "deseq2","corncob"),
    distr = c("norm", "mnorm", "identity"),
    output = c("cdf", "zscore", "identity"),
    adj = c(TRUE, FALSE),
    id = sim$id
))

eval_settings <- eval_settings %>% 
    filter(!(methods %in% c("deseq2", "corncob") & adj == TRUE)) %>% 
    filter(!(methods %in% c("deseq2", "corncob") & distr %in% c("norm", "mnorm"))) %>% 
    filter(!(methods %in% c("deseq2", "corncob") & output %in% c("cdf", "zscore")))

sim <- left_join(sim, eval_settings)

plan(multisession, workers = ncores)
tic()
sim$eval <- future_map(1:nrow(sim), .f = ~{
    source("../R/cilr.R")
    data <- readRDS(file = glue("{dir}/simulation_{i}.rds", 
                            dir = dir, i = sim$id[.x]))
    phydat <- sim2phylo(data)
    if (sim$methods %in% c("corncob", "deseq2")){
        scores <- diff_ab(phydat, method = sim$methods[.x], agg_level = "GENUS",
                        thresh = 0.05, padj = FALSE, return = "sig")
    } else {
        scores <- diff_ab(phydat, method = sim$methods[.x], agg_level = "GENUS",
                            thresh = 0.05, padj = FALSE, return = "sig", 
                            distr = sim$distr[.x], adj = sim$adj[.x], 
                            output = sim$output[.x])
    }
    evaluate <- c("fdr", "pwr")
    out <- sapply(evaluate, calculate_statistic, pred = as.vector(scores))
    names(out) <- evaluate
    return(out)
}, .options = furrr_options(seed = TRUE), .progress = TRUE)
toc()
plan(sequential)

saveRDS(sim, file = glue("{dir}/diff_ab_eval.rds", dir = dir))