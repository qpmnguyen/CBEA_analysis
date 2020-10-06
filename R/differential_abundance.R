library(tidyverse)
library(furrr)
library(tictoc)
library(ggsci)
library(progressr)
library(patchwork)
source("R/cilr.R")
source("R/simulations.R")
source("R/utils.R")






data <- zinb_simulation(n_samp = 300, b_spar = 0.2, b_rho = 0.1, eff_size = 3, n_tax = 5000, n_sets = 100, 
                        prop_set_inflate = 0.1, n_inflate = 50, prop_inflate = 0.5, samp_prop = 0.5)


sum_agg <- aggregate(data$X, data$A)

gsva_scores <- generate_alt_scores(X = data$X, A = data$A, method = "gsva", preprocess = T, pcount = 1)

cilr_scores <- simple_cilr(data$X, A = data$A, preprocess = T, pcount = 1, method = "sample")

mean_values <- apply(data$X, 2, mean)

plot(mean_values)

heatmap(cilr_scores)
heatmap(unclass(acomp(sum_agg)))
heatmap(unclass(clr(sum_agg)))
heatmap(gsva_scores)

test_result <- c()
for (i in 1:ncol(gsva_scores)){
  test_result[i] <- t.test(gsva_scores[1:150,i], gsva_scores[151:300,i])$p.value
}
test_result

mod <- get_diff_ab(cilr_scores, labels = data$label, method = "welch")
length(which(mod < 0.05))
sum(p.adjust(mod,"BH") <= 0.05)
which(p.adjust(mod, "BH") <= 0.05)


