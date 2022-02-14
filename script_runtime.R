library(targets)
library(tarchetypes)
library(CBEA)
library(tidyverse)
library(bench)
library(mia)
library(future)
library(future.callr)
source("R/simulations.R")

set.seed(1020)
plan(multisession)

values <- cross_df(list(
    n_perm = c(50, 100), 
    adj = c(TRUE, FALSE), 
    distr = c("mnorm", "norm")
))


data <- tar_target(df, {
    sim_dat <- zinb_simulation(n_samp = 500, spar = 0, s_rho = 0, eff_size = 1, 
                    n_tax = 800, n_inflate = 20, n_sets = 40, vary_params = TRUE)
    assay(sim_dat$obj, "Counts") <- assay(sim_dat$obj, "Counts") + 1
    sim_dat$obj <- mia::transformSamples(sim_dat$obj, abund_values = "Counts", method = "relabundance")
    sim_dat
})
mapping <- tar_map(values = values, unlist = FALSE, 
    tar_target(benchmarking, {
	args <- list(obj = df$obj, set = df$set, abund_values = "relabundance", 
			    distr = distr, adj = adj, output = "sig", n_perm = n_perm)
	if (distr == "mnorm"){
	    args <- c(args, list(control = list(fix_comp = "large")))
	}
	print(head(args))
        results <- bench::mark(
	    do.call(CBEA::cbea, args),
            memory = FALSE, iterations = 1)
	results 
    }),
    tar_target(eval_bench, {
	tibble(
	    adj = adj, 
	    distr = distr, 
	    n_perm = n_perm, 
	    time = benchmarking %>% dplyr::pull(median)
	)
    })
) 
combine <- tar_combine(combine_runtime, mapping[[2]], command = dplyr::bind_rows(!!!.x))
save <- tar_rds(save_runtime, saveRDS(combine_runtime, file = "output/runtime.rds"))


list(data, mapping, combine, save)
