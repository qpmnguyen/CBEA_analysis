library(caret)
library(rsample)
library(yardstick)
library(mlr3)
library(mlr3learners)
library(paradox)

source("R/utils.R")
source("R/simulations.R")
source("R/cilr.R")

dat <- sim_prediction(beta_eff = 1, snr = 2, sat = 0.2)
# TODO: make sure output has row and column names 
scores <- generate_alt_scores(dat$predictors$X,
            dat$predictors$A, method = "gsva")
scores <- cilr(dat$predictors$X, dat$predictors$A,
                output = "zscore", distr = "norm",
                resample = T)

scores <- generate_alt_scores(dat$predictors$X,
                              dat$predictors$A, method = "prop")

methods <- c("cilr", "gsva", "clr")
for (i in 1:length(methods)){
    #scores <- 
}



comb <- cbind(dat$outcome, scores)
colnames(comb)[1] <- "outcome"
colnames(comb) <- c("outcome", paste0("set",1:40))
comb <- as.data.frame(comb)
task <- TaskRegr$new(id = "scores", backend = comb, target = "outcome")
learner <- mlr_learners$get("regr.ranger")
learner$param_set$values <- list(num.trees = 1000, mtry = round(ncol(comb-1)/3,0))



train_set <- sample(task$nrow, 0.8 * task$nrow)
test_set = setdiff(seq_len(task$nrow), train_set)
learner$train(task, row_ids = train_set)
measure <- msr("regr.rsq")
print(learner$model)
prediction <- learner$predict(task, row_ids = test_set)
prediction$score(measure)


prediction_eval <- function(outcome, predictors, method=c("xgboost","rf"), nfolds=10){
    method <- match.arg(method)
    if(method == "xgboost"){
        tune_grid <- expand.grid(
            nrounds = seq(from = 200, to = 1000, by = 100), # number of trees
            eta = c(0.05, 0.1, 0.2), # learning rate 
            gamma = 0, # min_split_loss - thresholding for tree pruning
            colsample_bytree = c(1.0), # column subsample 
            subsample = c(1.0),# row subsample
            max_depth = c(4, 6, 8, 10), # max tree depth,
            min_child_weight = 3
        )
    } else {
        tune_grid <- data.frame(mtry = sqrt(ncol(predictors)))

    }

}