
sim_grid <- sim_grid[1,] %>% dplyr::select(-rep)
sim_dat <- do.call(zinb_simulation, as.list(sim_grid))

sim <- proc_sim(sim_dat)

begin <- Sys.time()
diff_ab(obj = sim$obj, 
        sets = sim$set, 
        make_phylo_manual = TRUE,
        abund_values = "Counts", 
        method = "cbea",
        distr = "mnorm", 
        adj = TRUE, 
        output = "zscore", n_perm = 50, 
        eval = "rset", thresh = 0.05, return = "sig")
end <- Sys.time()
