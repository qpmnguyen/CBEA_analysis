library(compositions)


source("R/simulations.R")
source("R/utils.R")

spar <- seq(0,0.9,by = 0.1)
set_size <- c(25,50,100,150)

data_list <- list()
for (i in 1:length(spar)){
  data_list[[i]] <- list()
  for (j in 1:10){
    abundance <- abundance_sim(n_tax = 1000, n_samp = 2000, sparsity = spar[i])
    tax_table <- 
    output <- 
  }
}
