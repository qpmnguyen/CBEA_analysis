library(tidyverse)
library(mia)
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
    sample <- as(colData(physeq), "data.frame")
    sample <- sample %>% rownames_to_column(var = "sample_id") %>% 
        mutate(label = if_else(HMP_BODY_SUBSITE == "Supragingival Plaque",1,0)) %>% 
        dplyr::select(sample_id, label) %>% as_tibble()
    outcomes <- left_join(results, sample, by = "sample_id") %>% as_tibble()
    return(outcomes)
}


#' Generate scores
#' @param set BiocSet
#' @param method Either ssgsea, gsva, wilcoxon or auc_scores
#' @param metric Either auc or inference
#' @param ... Additional arguments passed on to CBEA
enrichment_analysis <- function(physeq, set, method, 
                                label = NULL, metric, abund_values = "16SrRNA", ...) {
    if (method %in% c("ssgsea", "gsva")) {
        # for ssgsea and gsva only pseudocount
        assay(physeq, abund_values) <- assay(physeq, abund_values) + 1
        scores <- alt_scores_physeq(physeq, set, method = method, preprocess = FALSE, 
                                    abund_values = abund_values)
        scores <- scores %>% rownames_to_column(var = "sample_id") %>% as_tibble()
    } else if (method == "wilcoxon") {
        # wilcoxon with raw counts but adding pseudocounts to avoid zeroes 
        assay(physeq, abund_values) <- assay(physeq, abund_values) + 1
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
        args <- list(...)
        if ("control" %in% names(args)){
            if (args$distr %in% c("norm", "lst")){
                args$control <- NULL
            }
        } 
        # cbea requries transformation to proportions
        assay(physeq, abund_values) <- assay(physeq, abund_values) + 1
        physeq <- transformCounts(physeq, abund_values = abund_values, method = "relabundance", 
                                  name = "main_input")
        
        args <- c(args, list(
            obj = physeq, 
            set = set, 
            thresh = 0.05, 
            abund_values = "main_input"
        ))
        model <- do.call(cbea, args)
        print("Fitted model")
        scores <- tidy(model)
        if (!"sample_id" %in% colnames(scores)){
            scores <- scores %>% tibble::add_column(sample_id = colnames(assay(physeq, "main_input")), 
                                                    .before = 1)
        }
    }
    return(scores)
}

#' Getting random gene sets of different sizes
#' @param physeq Phyloseq object containing the data
#' @param size The size of the set
#' @param n_sets Number of sets of that size
get_rand_sets <- function(physeq, size, n_sets) {
    taxa <- rownames(rowData(physeq))
    set_list <- purrr::map(seq_len(n_sets), ~ {
        sets <- sample(taxa, size = size)
    })
    names(set_list) <- paste0("Set", seq_along(set_list))
    sets <- BiocSet::BiocSet(set_list)
    return(sets)
}
