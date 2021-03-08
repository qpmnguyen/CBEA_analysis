
grid <- tar_target(grid, {
    generate_grid(eval = "fdr")[1:3,]
})

sim <- tar_target(sim, {
    zinb_simulation(n_samp = grid$n_samp, spar = grid$spar, s_rho = grid$s_rho, eff_size = grid$eff_size, 
                    samp_prop = grid$samp_prop, vary_params = grid$vary_params, n_tax = grid$n_tax, 
                    n_inflate = grid$n_inflate, n_sets = grid$n_sets, prop_set_inflate = grid$prop_set_inflate, 
                    method = "normal")
}, pattern = map(grid), iteration = "list")

hypo_grid <- tar_target(hypo_grid, {
    eval_settings <- cross_df(list(
        model = "cilr",
        distr = c("mnorm", "norm"),
        adj = c(TRUE, FALSE)
    ))
    wlcx <- tibble(model = "wilcox")
    eval_settings <- full_join(eval_settings, wlcx, by = c("model"))
    eval_settings
})

eval <- tar_target(eval, {
    res <- analysis(sim = sim, model = hypo_grid$model, distr = hypo_grid$distr, 
                    eval = "fdr", adj = hypo_grid$adj)
    dplyr::bind_cols(grid, hypo_grid, res)
}, pattern = cross(map(sim), hypo_grid))

# res1 <- analysis(sim[[1]], model = "cilr", distr = "norm", eval = "fdr", adj = T)
# res2 <- analysis(sim[[1]], model = "cilr", distr = "norm", eval = "fdr", adj = T)