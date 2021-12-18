#library(teaR)
library(bench)
library(tidyverse)
library(targets)
library(tarchetypes)
library(future.callr)
source("functions.R")


set.seed(1020)
targets::tar_option_set(error = "workspace")

plan(callr)
grid <- tibble(
    n_samp = c(1000, 2000, 3000, 1000, 1000, 1000),
    n_tax = c(500, 500, 500, 500, 2500, 5000),
    eval = c("samp", "samp", "samp", "tax", "tax", "tax")
)
grid$id <- seq(1, nrow(grid))
saveRDS(grid, file = "output/grid_sim.rds")

evaluation_grid <- tar_target("eval_grid",
    purrr::cross_df(list(
        distr = c("mnorm", "norm"), 
        adj = c(TRUE, FALSE)))
)


sim_and_eval <- tar_map(
    values = grid, 
    unlist = FALSE,
    names = "id",
    tar_target("simulation", {
        zinb_simulation(n_samp = n_samp, n_tax = n_tax, n_sets = n_tax/50, 
                        s_rho = 0.5, eff_size = 1, spar = 0.2)
    }), 
    tar_target("transformation", {
        output <- vector(mode = "list", length = 2)
        output$X <- add_pseudocount(simulation$X)
        output$A <- A2list(simulation$A)
        return(output)
    }),
    tar_target("benchmark", {
        res <- bench::bench_time(cilr(ab_tab = transformation$X, set_list = transformation$A, 
                         output = "cdf", distr = eval_grid$distr, adj = eval_grid$adj))
        bind_cols(id = id, distr = eval_grid$distr, adj = eval_grid$adj, time = unname(res[2]))
    }, pattern = map(eval_grid))
)

combine_results <- tar_combine("combined_summaries", sim_and_eval[[3]], 
                               command = dplyr::bind_rows(!!!.x))

save_results <- tar_rds("output", saveRDS(combined_summaries, file = "output/results.rds"))

list(evaluation_grid, sim_and_eval, combine_results, save_results)