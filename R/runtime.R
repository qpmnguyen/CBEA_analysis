library(CBEA)
library(tidyverse)
library(bench)
library(furrr)
library(future)
source("R/simulations.R")
source("R/plot_utils.R")
data(hmp_gingival)
# benchmarking the same data set 
results <- bench::press(adj = c(TRUE, FALSE), distr = c("norm", "mnorm"), n_perm = c(50,100,200), {
    bench::mark(
        cbea(hmp_gingival$data, hmp_gingival$set, abund_values = "16SrRNA", distr = distr, adj = adj, output = "sig", 
             n_perm = n_perm)
    )    
})
res_plot <- results %>% dplyr::select(adj, distr, n_perm, median)

res_plot <- res_plot %>% 
    mutate(distr = recode(distr, "norm" = "Normal", "mnorm" = "Mixture Normal")) %>% 
    dplyr::rename("Distribution" = distr, "Correlation adjust" = adj) 
saveRDS(res_plot, "output/runtime_hmp_gingival.rds")

ggplot(res_plot, aes(x = factor(n_perm), y = median)) + geom_bar(stat = "identity") +
    facet_grid(Distribution~`Correlation adjust`, labeller = label_both) + 
    geom_text(aes(label = median, vjust = -0.5)) + 
    my_pretty_theme



nsamp_grid <- cross_df(list(
    n_samp = c(100, 500, 1000),
    adj = c(TRUE, FALSE), 
    distr = c("norm", "mnorm")
))


# benchmarking some results 
plan(multisession, workers = 3)
results_nsamp <- furrr::future_pmap(nsamp_grid, function(n_samp, adj, distr){
    df <- zinb_simulation(n_samp = n_samp, spar = 0.1, s_rho = 0.2, eff_size = 1, 
                          b_rho = 0.1, n_tax = 800,
                          n_inflate = 20, n_sets = 40, prop_set_inflate = 0.5)
    assay(df$obj, "Counts") <- assay(df$obj, "Counts") + 1
    bench::mark(
        cbea(df$obj, df$set, abund_values = "Counts", distr = distr, 
             adj = adj, output = "sig"), 
        iterations = 1, memory = FALSE,
    )
})

saveRDS(results_nsamp, "output/runtime_nsamp.rds")

results_nsets <- bench::press(nsets = c(50, 100, 150), adj = c(TRUE, FALSE), 
                              distr = c("norm", "mnorm"), {
    df <- zinb_simulation(n_samp = 300, spar = 0.1, s_rho = 0.2, eff_size = 1, b_rho = 0.1, 
                          n_tax = n_sets * 25, n_inflate = 25, n_sets = n_sets)
    bench::mark(
        cbea(df$obj, df$set, abund_values = "Counts", distr = distr, adj = adj, 
             output = "sig"), 
        iterations = 1, memory = FALSE
    )
})
saveRDS(results_nsets, "output/runtime_nsets.rds")

results_nperm <- bench::press(nperm = c(100,200,300), adj = c(TRUE, FALSE), 
                              distr = c("norm", "mnorm"), {
    df <- zinb_simulation(n_samp = 300, spar = 0.1, s_rho = 0.2, eff_size = 1, b_rho = 0.1, 
                          n_tax = 800, n_sets = 40, n_inflate = 20, prop_set_inflate = 0.5)
    bench::mark(
        cbea(df$obj, df$set, abund_values = "Counts", distr = distr, adj = adj, output = "sig", 
             n_perm = n_perm), 
        iterations =1, memory = FALSE
    )
})
saveRDS(results_nperm, "output/runtime_nperm.rds")
