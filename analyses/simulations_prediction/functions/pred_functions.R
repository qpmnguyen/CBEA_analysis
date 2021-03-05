library(tidymodels)
library(tidyverse)

source("../../R/cilr.R")
source("../../R/utils.R")
source("../../R/simulations.R")


# Function that takes a physeq object and return a data frame and a taxa table 
# data frame contains the column label, which has the label of the possible prediction task 
#' @param lab_col Name of the column detailing the label 
#' @param case_label Name of the label that should be indicated as a case (i.e. equals 1) for a binary prediction problem
#' @param control_label Name of the label that stands for control  
process_pred <- function(physeq, lab_col, case_label, data_type, control_label=NULL){
    label <- sample_data(physeq) %>% as("data.frame") %>% dplyr::select(!!lab_col)
    if(is.null(control_label)){
        label <- ifelse(label == case_label, 1, 0) 
    } else {
        label <- dplyr::case_when(label == case_label  ~ 1, 
                                  label == control_label ~ 0,
                                  TRUE ~ NA_real_) 
    }
    if (data_type == "wgs"){
        data <- otu_table(physeq) %>% as("matrix") %>% t() %>% 
            as.data.frame() %>% 
            dplyr::select(starts_with("s__")) 
        table <- tax_table(physeq) %>% unclass() %>% as.data.frame() %>% 
            rownames_to_column(var = "colnames") %>% filter(str_detect(colnames, "s__")) %>%
            column_to_rownames(var = "colnames")
    } else if (data_type == "16S"){
        data <- otu_table(physeq) %>% as("matrix") %>% t() %>% as.data.frame()
        table <- tax_table(physeq) %>% unclass() %>% as.data.frame()
    }
    na_idx <- which(is.na(label))
    if (length(na_idx) != 0){
        label <- label[-na_idx]
        data <- data[-na_idx,]
    }
    return(list(data = data, label = label %>% as.factor(), table = table))
}


# Generating a workflowr from a data object for the purposes of execution 
generate_wkflow <- function(data, task = c("regression", "classification")){
    task <- match.arg(task)
    if (task == "regr"){
        mtry <- round(sqrt(ncol(data) - 1),0)
    } else {
        mtry <- round((ncol(data)-1)/3,0)
    }
    rf <- rand_forest(mtry = mtry, trees = 1000) %>%
        set_engine("ranger", importance = "impurity") %>%
        set_mode(task)
    wkflow <- workflow() %>% add_model(rf) %>% add_formula(label ~ .)
    return(wkflow)
}

# Function to generate aggregation from a wgs processed object 
generate_aggregation <- function(proc_object, method = c("cilr", "gsva", "clr", "ssgsea"), ...){
    proper_names <- c("outcome", "predictors", "beta")
    # Checking sure that all elements are there 
    if (!is.list(proc_object)) rlang::abort(message = "Object needs to be a list")
    if (length(proc_object) != 3) rlang::abort(message = "Object needs all three elements")
    if (length(intersect(proper_names, names(proc_object))) != 3){
        rlang::abort(message = "Object needs to be properly labeled")
    }
    if (length(intersect(c("X", "A"), names(proc_object$predictors))) != 2){
        rlang::abort(message = "Predictor needs to have both X and A matrix")
    }
    
    # Get matrix from data 
    data <- proc_object$predictors$X
    A <- proc_object$predictors$A
    # drop single tons 
    counts <- apply(A, 2, sum)
    idx_singletons <- which(counts <= 1)
    if (length(idx_singletons) >= 1){
        A <- A[,-idx_singletons]    
    }
    
    # Then generate aggregations  
    if (method == "cilr"){
        new_dat <- cilr(X = data, A = A, resample = T, preprocess = T, pcount = 1, ..., 
                        maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
    } else if (method == "clr"){
        new_dat <- aggregate(X = data, A = A)
        new_dat[new_dat == 0] <- 1 
        new_dat <- unclass(clr(new_dat))
    } else {
        new_dat <- generate_alt_scores(X = data, A = A, method = method, preprocess = T, pcount = 1, ...)
    }
    
    new_dat <- as.data.frame(new_dat)
    new_dat <- new_dat %>% mutate(label = proc_object$outcome) %>% 
        dplyr::select(label, everything())
    return(new_dat)
}


fit_and_eval <- function(data, nfolds = 10, task = "classification"){
    wkflow <- generate_wkflow(data, task = task)
    folds <- vfold_cv(data, v = nfolds)
    eval <- wkflow %>% fit_resamples(folds)
    metrics <- collect_metrics(eval)
    return(metrics)
}

