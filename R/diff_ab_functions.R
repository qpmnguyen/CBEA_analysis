# This file holds functions to test for differential abundance 
# Quang Nguyen 

library(furrr)
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(corncob)

#' Interface function to perform differential abundance test 
#' @param method The method used to test differential abundance
#' @param thresh Threshold for adjustment 
#' @param adj Whether to perform p-value adjustment 
#' @param output What type of output "pvalue" or "significance" 
#' @param ... Additional arguments passed to the differential abundance functions
diff_ab <- function(physeq, method = c("cilr_welch", "cilr_wilcox", 
                                       "corncob", "deseq2"),
                    thresh, agg_level,
                    padj = FALSE, return = c("pvalue", "sig"), ...){
    method <- match.arg(method)
    return <- match.arg(return)
    if (method %in% c("corncob", "deseq2")){
        physeq <- tax_glom(physeq, taxrank = "GENUS")
        physeq <- transform_sample_counts(physeq, function(x) ifelse(x == 0, 1, x))
    } 
    results <- model_interface(physeq, method, agg_level, ...)
    if (padj == TRUE){
        results <- p.adjust(p = results, method = "BH")
    }
    if (return == "pvalue"){
        return(results)
    } else if (return == "sig"){
        return((results <= thresh)*1)
    }
}

#' Function to run models, with the defining feature labelled "group"
#' TODO: Add functionality to have different covariates being the discerning variable 
#' @param physeq Phyloseq type object containing taxtable and data points 
#' @param agg_level Level to aggregate variables too 
#' @param ... Additional arguments to assed to differentialTest, DESeq and cilr
model_interface <- function(physeq, method = c("cilr_welch", "cilr_wilcox",
                                       "corncob", "deseq2"), agg_level, ...){
    method <- match.arg(method)
    message(glue("Running {model}", model = method))
    if (method == "deseq2"){
        def <- list(test = "LRT", fitType = "local", reduced = ~1)
        if (missing(...)){
            args <- def
        } else {
            sup <- list(...)
            args <- merge_lists(defaults = def, supplied = sup)
        }
        deseq <- phyloseq_to_deseq2(physeq, ~factor(group))
        args$object <- deseq
        mod <- do.call(DESeq2::DESeq, args)
        res <- DESeq2::results(mod)
        sig <- res$pvalue
        names(sig) <- as.vector(tax_table(physeq)[rownames(res), agg_level])
    } else if (method == "corncob"){
        def <- list(formula = ~ group, phi.formula = ~ group, formula_null = ~1, 
                    phi.formula_null = ~ group, test = "LRT", boot = FALSE, 
                    fdr_cutoff = 0.05)
        if (missing(...)){
            args <- def
        } else {
            sup <- list(...)
            args <- merge_lists(defaults = def, supplied = sup)
        }
        args$data <- physeq
        mod <- do.call(corncob::differentialTest, args)
        sig <- mod$p
        names(sig) <- as.vector(tax_table(physeq)[names(sig), agg_level])
    } else if (stringr::str_detect(method, 'cilr')){
        def <- list(resample = T, distr = "norm", nperm = 5, adj = T, output = "zscore",
                    pcount = 1, transform = "prop",
                    maxrestarts=1000, epsilon = 1e-06, maxit= 1e5 )
        if (missing(...)){
            args <- def
        } else {
            sup <- list(...)
            args <- merge_lists(defaults = def, supplied = sup)
        }
        A <- taxtab2A(physeq@tax_table, agg_level = agg_level)
        X <- as(physeq@otu_table, "matrix")
        args$A <- A 
        args$X <- X 
        label <- physeq@sam_data$group
        idx <- which(label == 1)
        scores <- do.call(cilr, args)
        if (stringr::str_detect(method, "wilcox")){
            sig <- map_dbl(seq(ncol(scores)), .f = ~{
                wilcox.test(scores[-idx, .x], scores[idx, .x])$p.value
            })
        } else if (stringr::str_detect(method, "welch")) {
            sig <- map_dbl(seq(ncol(scores)), .f = ~{
                t.test(scores[-idx, .x], scores[idx, .x])$p.value
            })
        }
        names(sig) <- colnames(A)
    } 
    return(sig)
}

