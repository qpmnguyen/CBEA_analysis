library(furrr)
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(corncob)
source("R/cilr.R")
source("R/utils.R")
source("R/simulations.R")

dat <- zinb_simulation(200, spar = 0.2, eff_size = 2, s_rho = 0, b_rho = 0, 
                            n_tax = 300, n_inflate = 50, n_sets = 6, 
                            prop_set_inflate = 1, prop_inflate = 0.5, samp_prop = 0.5, 
                            method = "compensation", vary_params=TRUE, parameters=NULL)
physeq <- sim2phylo(dat)

physeq <- tax_glom(physeq, taxrank = "GENUS")
physeq <- transform_sample_counts(physeq, function(x) ifelse(x == 0, 1, x))

model_interface <- function(physeq, method, thresh, agg_level, ...){
    message(glue("Running {model}", model = "method"))
    if (method == "deseq2"){
        def <- list(test = "LRT", fitType = "local", reduced = ~1)
        if (missing(...)){
            args <- def
        } else {
            sup <- list(...)
            args <- merge_lists(defaults = def, supplied = sup)
        }
        deseq <- phyloseq_to_deseq2(physeq, ~as.factor(group))
        args$object <- deseq
        mod <- do.call(DESeq2::DESeq, args)
        #mod <- DESeq2::DESeq(deseq, test = "LRT", fitType = "parametric", reduced = ~1)
        res <- DESeq2::results(mod)
        sig <- (res$pvalue < thresh) * 1 
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
        sig <- (mod$p < thresh) * 1
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
        scores <- do.call(cilr, def)
        if (stringr::str_detect(method, "wilcox")){
            f <- wilcox.test
        } else if (stringr::str_detect(method, "welch")) {
            f <- t.test
        }
        
    } 
    return(sig)
}

model_interface(physeq, "corncob", 0.05)

agg_level <- "GENUS"
thresh <- 0.05
diff_ab_eval <- function(physeq, method, thresh, ...){
    message(glue("Running {model}", model = "method"))

        
}
