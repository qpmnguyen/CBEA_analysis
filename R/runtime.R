library(CBEA)
library(bench)
source("R/simulations.R")

nsets_grid <- cross_df(list(
    n_samp = c(100,500,1000),
    adj = c(TRUE, FALSE),
    distr = c("norm", "mnorm")
))

# benchmarking some results 
results <- bench::press(n_samp = c(100,500, 1000), adj = c(TRUE, FALSE), distr = c("norm", "mnorm"), {
    df <- zinb_simulation(n_samp = n_samp, spar = 0.1, s_rho = 0.2, eff_size = 1, b_rho = 0.1, n_tax = 800,
                            n_inflate = 20, n_sets = 40, prop_set_inflate = 0.5)
    bench::mark(
        cbea(df$obj, df$set, abund_values = "Counts", distr = "norm", adj = FALSE, output = "sig"), 
        max_iterations = 10
    )
})

