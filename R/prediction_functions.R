library(caret)
library(rsample)
library(yardstick)

source("R/utils.R")
source("R/simulations.R")

dat <- sim_prediction(beta_eff = 1, snr = 2, sat = 0.2)




scores <- generate_alt_scores(dat$predictors$X, dat$predictors$A, method = "gsva")


prediction_eval <- function(outcome, predictors){
    
}