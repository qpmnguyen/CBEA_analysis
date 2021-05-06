# This file holds functions to test for differential abundance 
# Quang Nguyen 

library(tidyverse)
library(phyloseq)
library(DESeq2)
library(corncob)
library(glue)
library(binom)

source("../../R/utils.R")
source("../../R/cilr.R")

#' Interface function to perform differential abundance test 
#' @param method The method used to test differential abundance
#' @param thresh Threshold for adjustment 
#' @param adj Whether to perform p-value adjustment 
#' @param data_type Whether data is wgs or 16S. Wgs require additional processing 
#' @param output What type of output "pvalue" or "significance" 
#' @param ... Additional arguments passed to the differential abundance functions
diff_ab <- function(physeq, 
                    method = c("cilr_welch", "cilr_wilcox", "corncob", "deseq2"),
                    agg_level,
                    data_type,
                    thresh = 0.05, 
                    padj = FALSE, 
                    return = c("pvalue", "sig"), prune = TRUE, 
                    ...){
    method <- match.arg(method)
    return <- match.arg(return)
    if (data_type == "wgs"){
        tax_counts <- otu_table(physeq) %>% as("matrix") %>% as.data.frame() %>% 
            rownames_to_column(var = "tax") %>%  
            filter(str_starts(tax, "s_")) %>% column_to_rownames(var = "tax")
        otu_table(physeq) <- otu_table(tax_counts, taxa_are_rows = T)
    }
    
    # filtering out singletons 
    if (prune == TRUE){
        physeq <- taxtab_prune(physeq, agg_level = agg_level)
    }
    
    if (method %in% c("corncob", "deseq2")){
        # aggregating  
        physeq <- tax_glom(physeq, taxrank = agg_level)
        # adding a pseudocount 
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
#' REMEMBER, GROUP IS THE DISCERNING VARIABLE  
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
            # deseq2 doesn't handle extra arguments well 
            idx <- which(names(sup) %in% c("distr", "adj", "output"))
            sup <- sup[-idx]
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
    # Model interfacing for cilr 
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
        # Creating A matrix 
        A <- taxtab2A(physeq@tax_table, agg_level = agg_level, full = FALSE)
        # Creating X matrix 
        X <- as(physeq@otu_table, "matrix")
        if (taxa_are_rows(physeq)){
            X <- t(X)
        }
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
    } else if (stringr::str_detect(method, "gsva")){
        def <- list(pcount = 1, method = "gsva")
        if (missing(...)){
            args <- def
        } else {
            sup <- list(...)
            args <- merge_lists(defaults = def, supplied = sup)
        }
        A <- taxtab2A(physeq@tax_table, agg_level = agg_level, full = FALSE)
        # Creating X matrix 
        X <- as(physeq@otu_table, "matrix")
        if (taxa_are_rows(physeq)){
            X <- t(X)
        }
        args$A <- A 
        args$X <- X 
        label <- physeq@sam_data$group
        idx <- which(label == 1)
        scores <- do.call(generate_alt_scores, args)
        if (stringr::str_detect(method, "wilcox")){
            sig <- map_dbl(seq(ncol(scores)), .f = ~{
                wilcox.test(scores[-idx, .x], scores[idx, .x])$p.value
            })
        } else if (stringr::str_detect(method, "welch")) {
            sig <- map_dbl(seq(ncol(scores)), .f = ~{
                t.test(scores[-idx, .x], scores[idx, .x])$p.value
            })
        }
        #TODO: GSVA AUTOMATICALLY DROPS SETS WITH SIZE == 1. THIS MEANS THAT THE LENGTH OF SIG
        # WILL NOT BE EQUIVALENT TO 
        names(sig) <- colnames(A)
    }
    return(sig)
}

#' @param nvec Named vector of 1 and 0s for "sig" output from diff_ab
#' @return prop is the proportion to be zero 
eval_function <- function(nvec, ci = FALSE){
    prop <- sum(nvec == 1)/length(nvec)
    if (ci == TRUE){
        conf <- binom.confint(prop, length(nvec), conf.level = 0.95, methods = "ac")
        return(tibble(est = prop, upper = conf$upper, lower = conf$lower))
    } else {
        return(prop)
    }
}



