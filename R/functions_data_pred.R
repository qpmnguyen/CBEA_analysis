library(tidymodels)
library(tidyverse)
library(phyloseq)
library(themis)
library(BiocSet)
library(CBEA)
source("R/functions_data_diffab.R")
source("R/utils.R")

# Function to generate aggregation from a wgs processed object 
generate_aggregation <- function(physeq, method, ...){
    method <- match.arg(method, c("cbea", "gsva", "ssgsea", "clr"))
    
    label <- sample_data(physeq) %>% as("data.frame") %>% 
        dplyr::select(sample_name, outcome) %>% as_tibble()
    
    # first let's get the set and filter out all singular sets 
    set <- const_set_taxtable(table = tax_table(physeq), rank = "Genus")
    filt_set <- set %>% filter_set(size >= 2) %>% es_set() %>% pull(set)
    set <- set %>% filter_elementset(set %in% filt_set)
    
    if (method == "cbea"){
        physeq <- phyloseq::transform_sample_counts(physeq, function(x) x + 1)
        physeq <- phyloseq::transform_sample_counts(physeq, function(x) x/sum(x))
        obj <- mia::makeTreeSummarizedExperimentFromPhyloseq(physeq)
        args <- list(...)
        if ("distr" %in% names(args)){
            if (args$distr == "mnorm"){
                args <- c(args, list(control = list(fix_comp = "large")))
            }
        }
        args <- c(args, list(obj = obj, set = set, parametric = TRUE, abund_values = "counts"))
        mod <- do.call(cbea, args)
        print(mod)
        df <- tidy(mod)
        if (!"sample_name" %in% colnames(df)){
            df <- df %>% tibble::add_column(sample_name = colnames(assay(obj, "counts")), 
                                                    .before = 1)
        }
    } else if (method %in% c("ssgsea", "gsva")){
        physeq <- transform_sample_counts(physeq, function(x) x + 1)
        obj <- mia::makeTreeSummarizedExperimentFromPhyloseq(physeq)
        df <- alt_scores_physeq(obj, set, method = method, preprocess = FALSE, 
                                    abund_values = "counts")
        df <- df %>% rownames_to_column(var = "sample_name") %>% as_tibble()
    } else if (method == "clr"){
        physeq <- transform_sample_counts(physeq, function(x) x + 1)
        physeq <- speedyseq::tax_glom(physeq, "Genus")
        clr_transformed <- compositions::clr(compositions::acomp(t(as(otu_table(physeq), "matrix"))))
        df <- as.data.frame(clr_transformed)
        df <- df %>% rownames_to_column("sample_name") %>% as_tibble()
    }
    aug_df <- df %>% mutate(outcome = label$outcome)
    out <- list(df = df, label = label, aug_df = aug_df)
    return(out)
}


# Generating a workflowr from a data object for the purposes of execution 
generate_wkflow <- function(data, task = c("regression", "classification")){
    task <- match.arg(task)
    if (task == "regression"){
        mtry <- round(sqrt(ncol(data) - 1),0)
    } else {
        mtry <- round((ncol(data)-1)/3,0)
    }
    
    proc_rec <- recipe(outcome ~ ., data = data) %>% 
        step_rm("sample_name") %>% 
        step_smote(outcome) %>%
        step_normalize(all_numeric_predictors())
    
    rf <- rand_forest(mtry = mtry, trees = 1000) %>%
        set_engine("ranger", importance = "impurity") %>%
        set_mode(task)
    wkflow <- workflow() %>% add_model(rf) %>% add_recipe(proc_rec)
    return(wkflow)
}




fit_and_eval <- function(data, nfolds = 10, task = "classification"){
    wkflow <- generate_wkflow(data, task = task)
    folds <- vfold_cv(data, v = nfolds)
    eval <- wkflow %>% fit_resamples(folds)
    metrics <- collect_metrics(eval)
    return(metrics)
}


