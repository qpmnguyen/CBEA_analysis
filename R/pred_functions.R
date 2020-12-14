library(mlr3)
library(mlr3measures)
library(mlr3learners)
library(phyloseq)
source("R/utils.R")

#' @title Function to fit
#' @param object the simulation object or a phyloseq type object
#' @param type What is the learning task
#' @param ... Additional arguments to be passed to cilr.R
fit_prediction <- function(object, type=c("regr", "classif"), ...) {
    type <- match.arg(type)
    # making defaults for cilr
    def <- list(resample = T, distr = "norm", nperm = 5, adj = T, output = "zscore",
                    pcount = 1, transform = "prop",
                    maxrestarts=1000, epsilon = 1e-06, maxit= 1e5 )
    if (missing(...)){
        args <- def
    } else {
        sup <- list(...)
        args <- merge_lists(defaults = def, supplied = sup)
    }
    # return defaults
    if (class(object) == "phyloseq") {
        predictors <- as(object@otu_table, "matrix")
        response <- as.vector(object@sam_data$group)
        A <- taxtab2A(object@tax_table, "GENUS")
    } else {
        predictors <- object$predictors$X
        response <- object$outcome
        A <- object$predictors$A
    }
    comb <- as.data.frame(cbind(response, predictors))
    colnames(comb)[1] <- "outcome"
    # create task list for all outcomes
    methods <- c("cilr", "gsva", "clr")
    task_list <- vector(mode = "list", length = 3)
    # creating the list of tasks
    for (i in seq(length(methods)){
        if (methods[i] == "cilr"){
            scores <- do.call(predictors, args)
        } else if (methods[i] == "gsva"){
            scores <- generate_alt_scores(predictors, A, method = "gsva")
        } else if (methods[i] == "clr"){
            scores <- generate_alt_scores(predictors, A, method = "prop")
            scores <- unclass(clr(acomp(scores))) %>% as.data.frame()
        } else if (methods[i] == "everything"){
            scores <- unclass(clr(acomp(dat$predictors$X))) %>% as.data.frame()
        } 
        comb <- cbind(dat$outcome, scores)
        colnames(comb)[1] <- "outcome"
        task_list[[i]] <- create_task(type = type, id = methods[i], backend = comb, id = "outcome")
    }
    if (type == "classif"){
        learner <- lrn("classif.ranger")
        learner$predict_type <- "prob"
        mtry <- round(sqrt(ncol(predictors)),0)
        measure <- list(msr("classif.auc"), msr("classif.acc"))
    } else if (type == "regr"){
        learner <- mlr_learners$get("regr.ranger")
        mtry <- round(ncol(predictors)/3,0)
        measure <- list(msr("regr.rmse"), msr("regr.srho"))
    }
    learner$param_set$values <- list(num.trees = 1000, mtry = mtry)
    # getting values in a grid for benchmarking
    design <- benchmark_grid(
        tasks = task_list,
        learners = learner,
        resamplings = rsmp("cv", folds = 10)
    )
    # benchmarking values via 10 fold cross validation
    bmr <- benchmark(design)
    tab <- bmr$aggregate(measure)
    return(tab[,c(3,7,8)])
}

#' Function to create the tasks for mlr3
#' @param type What is the learning type
#' @param id A string indicating the identifier of the task
#' @param backend The dataframe backend of the task
#' @param target A strong for the target column of the task
create_task <- function(type = c("classif", "regr"), id, backend, target){
    type <- match.arg()
    if (type == "classif"){
        task <- mlr3::TaskClassif$new(id, backend, target)
    } else if (type == "regr"){
        task <- mlr3::TaskRegr$new(id, backend, target)
    }
    return(task)
}