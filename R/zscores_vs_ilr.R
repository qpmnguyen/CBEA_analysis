### SIDE TRACKING ILR vs Z-score  
pwr_sim <- create_parameters(list(
  rep = seq(1,10),
  b_rho = c(0.1, 0.2, 0.5),
  eff_size = c(2,3,4,5),
  n_inflate = c(50,100,150,200)
))

plan(multiprocess, workers = round(availableCores()/2,0))
tic()
opt <- furrr_options(seed = T)
with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$sim <- furrr::future_map(pwr_sim$param, ~{
    p()
    zinb_simulation(n_samp = 300, b_spar = 0.4, b_rho = .x$b_rho, 
                    eff_size = .x$eff_size, n_inflate = .x$n_inflate, rho_ratio = 1)
  }, .options = opt)
})
toc()
plan(sequential)

plan(multiprocess, workers = round(availableCores()/2,0))
with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$scores_cilr <- future_map(pwr_sim$sim, .f = ~{
    p()
    simple_cilr(X = .x$X, A = .x$A, preprocess = T, pcount = 1, transform = "prop", abs = F)
  })
})
with_progress({
  p <- progressor(steps = nrow(pwr_sim))
  pwr_sim$scores_zscore <- future_map(pwr_sim$sim, .f = ~{
    p()
    generate_alt_scores(X = .x$X, A = .x$A, method = "zscore", preprocess = T, transform="clr", pcount=1)
  })
})
plan(sequential)

scores <- pwr_sim %>% dplyr::select(c(sim, starts_with("scores")))
auc <- vector(mode = "list", length = length(3:ncol(scores))-1)
for (i in 3:ncol(scores)){
  new_name <- glue("auc_{name}", name = colnames(scores)[i])
  values <- map_dbl(1:nrow(scores), .f = function(.x){
    stat <- calculate_statistic(eval = "auc", pred = scores[.x,i][[1]], 
                                true = scores$sim[[.x]]$label)
    return(stat)
  })
  auc[[new_name]] <- values
}
auc <- do.call(cbind, auc)
auc <- as_tibble(auc)
pwr_sim <- cbind(pwr_sim, auc)

plotting_data <- pwr_sim %>% dplyr::select(c(id,param, starts_with("auc"))) %>% pivot_longer(starts_with("auc")) %>% 
  unnest(param) %>% group_by(n_inflate, b_rho, eff_size, name) %>% summarise(min = min(value, na.rm = T), 
                                                                          mean = mean(value, na.rm = T),
                                                                          max = max(value, na.rm = T))
plotting_data <- plotting_data %>% dplyr::rename("Set Size" = "n_inflate", "Inter-taxa correlation" = "b_rho")

auc_plot <- ggplot(plotting_data, aes(x = eff_size, y = mean, col = name)) + geom_hline(yintercept = 0.8, col = "red") + 
  geom_point(size = 2) + facet_grid(`Set Size` ~ `Inter-taxa correlation`, labeller = label_both) + theme_bw() + 
  geom_line() + scale_color_nejm(labels = c("cILR", "z-score")) + 
  labs(col = "Models", x = "Effect Size", y = "Mean AUC (10 data sets)")
ggsave(auc_plot, filename = "docs/manuscript/figures/zscore_vs_cilr.png", dpi = 300, width = 8, height = 8)
