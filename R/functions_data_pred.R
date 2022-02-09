library(tidymodels)
library(tidyverse)
library(phyloseq)
library(themis)
library(BiocSet)
library(CBEA)
source("R/functions_data_diffab.R")
source("R/utils.R")

# Function to generate aggregation from a wgs processed object 
generate_aggregation <- function(physeq, method, set = NULL,...){
    method <- match.arg(method, c("cbea", "gsva", "ssgsea", "clr"))
    
    # extract label from sample data 
    label <- sample_data(physeq) %>% as("data.frame") %>% 
        dplyr::select(sample_name, outcome) %>% as_tibble()
    if (is.null(set)){
        # first let's get the set and filter out all singular sets 
        set <- const_set_taxtable(table = tax_table(physeq), rank = "Genus")
        filt_set <- set %>% filter_set(size >= 2) %>% es_set() %>% pull(set)
        set <- set %>% filter_elementset(set %in% filt_set)
    }
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
        if (is.null(set)){
            physeq <- speedyseq::tax_glom(physeq, "Genus")
        } else {
            grp_vec <- rep(NA, ntaxa(physeq))
            names(grp_vec) <- taxa_names(physeq)
            set_names <- set %>% es_set() %>% pull(set)
            for (i in seq_along(set_names)){
                el_names <- set %>% es_elementset() %>% filter(set == set_names[i]) %>% 
                    pull(element)
                idx <- which(names(grp_vec) %in% el_names)
                grp_vec[idx] <- set_names[i]
            }
            grp_vec <- as.vector(grp_vec)
            physeq <- speedyseq::merge_taxa_vec(physeq, grp_vec)
        }
        clr_transformed <- compositions::clr(compositions::acomp(t(as(otu_table(physeq), "matrix"))))
        df <- as.data.frame(clr_transformed)
        df <- df %>% rownames_to_column("sample_name") %>% as_tibble()
    }
    aug_df <- df %>% mutate(outcome = label$outcome)
    out <- list(df = df, label = label, aug_df = aug_df)
    return(out)
}


# Generating a workflowr from a data object for the purposes of execution 
generate_wkflow <- function(data, task = c("regression", "classification"), unbal_class = TRUE){
    task <- match.arg(task)
    if (task == "regression"){
        mtry <- round(sqrt(ncol(data) - 1),0)
    } else {
        mtry <- round((ncol(data)-1)/3,0)
    }
    
    proc_rec <- recipe(outcome ~ ., data = data) %>% 
        step_rm("sample_name") %>% 
        step_normalize(all_numeric_predictors())
    
    if (unbal_class == TRUE){
        proc_rec <- proc_rec %>% step_smote(outcome)
    }

    
    rf <- rand_forest(mtry = mtry, trees = 1000) %>%
        set_engine("ranger") %>%
        set_mode(task)
    wkflow <- workflow() %>% add_model(rf) %>% add_recipe(proc_rec)
    return(wkflow)
}


fit_and_eval <- function(data, nfolds = 10, task = c("regression", "classification"), unbal_class = TRUE){
    wkflow <- generate_wkflow(data, task = task, unbal_class = unbal_class)
    folds <- vfold_cv(data, v = nfolds)
    eval <- wkflow %>% fit_resamples(folds)
    metrics <- collect_metrics(eval)
    return(metrics)
}

proc_sim <- function(simulation_dat){
    physeq <- simulation_dat$predictors$X
    physeq <- mia::makePhyloseqFromTreeSummarizedExperiment(physeq, abund_values = "Counts")
    sample_names(physeq) <- paste0("samp", seq_len(nsamples(physeq)))
    sample_data(physeq)$sample_name <- sample_names(physeq)
    return(list(
        physeq = physeq, 
        set = simulation_dat$predictors$set
    ))
}


