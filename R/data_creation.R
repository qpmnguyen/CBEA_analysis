library(tidyverse)
library(furrr)
library(progressr)
library(optparse)
library(tictoc)
source("simulations.R")

option_list <- list(
  make_option("--ncores", type = "integer", default=5, help="Number of cores"),
  make_option("--export", type = "string", help = "Filename for export")
)

opt <- parse_args(OptionParser(option_list = option_list))

sim <- list(
  rep = seq(1,100),
  b_spar = c(0.2, 0.4, 0.6, 0.8),
  b_rho = c(0.1, 0.2, 0.5),
  n_inflate = c(50,100,150,200)
)


sim <- create_parameters(list(sim))

tic()
plan(multicore, workers = opt$ncores)
opt <- furrr_options(seed = T)
with_progress({
  p <- progressor(steps = nrow(sim))
  sim$sim <- furrr::future_map(sim$param, ~{
    p()
    suppressMessages(zinb_simulation(n_samp = 1000, b_spar = .x$b_spar, b_rho = .x$b_rho, 
                    eff_size = 1, n_inflate = 50, rho_ratio = 1, n_tax = 1000))
  }, .options = opt)
})
plan(sequential)
toc()

saveRDS(sim, file = glue("parameters_{name}.rds", name = opt$export))