library(phyloseq)
library(tidyverse)
library(CBEA)
library(yardstick)
source("R/utils.R")


#' Function to preprocess the gingival data set
gingival_processing <- function(data) {
    physeq <- data$physeq
    annotation <- data$annotation
    colnames(annotation) <- c("GENUS", "METABOLISM")
    tax_table <- as.data.frame(as(physeq@tax_table, "matrix")) %>% 
        rownames_to_column(var = "id") %>% as_tibble()
    full_table <- full_join(tax_table, annotation, by = "GENUS") %>% 
        filter(!is.na(id) & !is.na(METABOLISM))
    sets <- const_set(full_table, id = "id", member = "METABOLISM")
    #physeq <- filter_taxa(physeq, function(x) sum(x > 0)/length(x) >= 0.05, TRUE)
    sets <- unify_sets(physeq, sets)
    return(list(
        physeq = physeq,
        sets = sets
    ))
}

#' @param results Output of enrichment_analysis  
#' @param physeq This is the original physeq 
gingival_evaluate <- function(physeq, results){
    sample <- as(sample_data(physeq), "data.frame")
    sample <- sample %>% rownames_to_column(var = "sample_id") %>% 
        mutate(label = if_else(HMP_BODY_SUBSITE == "Supragingival Plaque",1,0)) %>% 
        dplyr::select(sample_id, label) %>% as_tibble()
    outcomes <- left_join(results, sample, by = "sample_id")
    return(outcomes)
}



#' Generate scores
#' @param set BiocSet
#' @param method Either ssgsea, gsva, wilcoxon or auc_scores
#' @param metric Either auc or inference
#' @param ... Additional arguments passed on to CBEA
enrichment_analysis <- function(physeq, set, method, 
                                label = NULL, metric, ...) {
    if (method %in% c("ssgsea", "gsva")) {
        # for ssgsea and gsva only pseudocount
        physeq <- transform_sample_counts(physeq, function(x) x + 1)
        scores <- alt_scores_physeq(physeq, set, method = method, preprocess = FALSE)
        scores <- scores %>% rownames_to_column(var = "sample_id") %>% as_tibble()
    } else if (method == "wilcoxon") {
        # wilcoxon with raw counts but adding pseudocounts to avoid zeroes 
        physeq <- transform_sample_counts(physeq, function(x) x + 1)
        if (metric == "auc") {
            output <- "scores"
        } else {
            output <- "sig"
        }
        scores <- wc_test_physeq(
            physeq = physeq, set = set,
            thresh = 0.05, alt = "two.sided",
            preprocess = FALSE, ...
        )
    } else if (method == "cbea") {
        # cbea requries transformation to proportions
        physeq <- transform_sample_counts(physeq, function(x) x + 1)
        #physeq <- transform_sample_counts(physeq, function(x) x / sum(x))
        scores <- CBEA::cbea(obj = physeq, set = set, thresh = 0.05, ...)
    }
    return(scores)
}

#' Old enrichment analysis function that does not use the CBEA package
# enrichment_analysis <- function(X, A, method, label, metric, ...){
#   if (method %in% c("ssgsea", "gsva")){
#     scores <- generate_alt_scores(X = X, A = A, method = method, preprocess = T, pcount = 1)
#   } else if (method %in% c("wilcoxon")){
#     if (metric == "auc"){
#       output <- "scores"
#     } else {
#       output <- "sig"
#     }
#     scores <- wc_test(X = X, A = A, thresh = 0.05, preprocess = T, pcount = 1, output = output)
#   } else {
#     scores <- cilr(X = X, A = A, resample = T, ..., maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
#   }
#   is.matrix(scores)
#   output <- enrichment_evaluate(scores = scores, results = label, metric = metric)
#   return(output)
# }

#' Getting random gene sets of different sizes
#' @param physeq Phyloseq object containing the data
#' @param size The size of the set
#' @param n_sets Number of sets of that size
get_rand_sets <- function(physeq, size, n_sets) {
    taxa <- taxa_names(physeq)
    set_list <- purrr::map(seq_len(n_sets), ~ {
        sets <- sample(taxa, size = size)
    })
    names(set_list) <- paste0("Set", seq_along(set_list))
    sets <- BiocSet::BiocSet(set_list)
    return(sets)
}
