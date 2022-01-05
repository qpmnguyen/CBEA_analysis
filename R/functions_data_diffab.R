library(DESeq2)
library(corncob)
library(glue)
library(binom)

#' Generate scores
#' @param set BiocSet
#' @param method Either ssgsea, gsva, wilcoxon or auc_scores
#' @param metric Either auc or inference
#' @param ... Additional arguments passed on to CBEA
enrichment_analysis <- function(physeq, set, method, 
                                label = NULL, metric, preprocess = TRUE, ...) {
    if (method %in% c("gsva", "wilcoxon")){
        preprocess <- FALSE
        message("Forcing preprocess to be false when evaluating GSVA")
    }
    if (preprocess == TRUE){
        physeq <- transform_sample_counts(physeq, function(x) x + 1)
        physeq <- transform_sample_counts(physeq, function(x) x / sum(x))
    }
    if (method %in% c("ssgsea", "gsva")) {
        scores <- alt_scores_physeq(physeq, set, method = method, preprocess = TRUE)
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
            preprocess = TRUE, ...
        )
    } else if (method == "cbea") {
        scores <- CBEA::cbea(obj = physeq, set = set, thresh = 0.05, ...)
    }
    return(scores)
}