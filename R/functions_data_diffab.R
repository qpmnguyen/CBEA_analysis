library(DESeq2)
library(corncob)
library(glue)
library(binom)
library(CBEA)
requireNamespace("speedyseq", quietly = TRUE)

# retrieving 
source("R/functions_data_ss.R")
source("R/utils.R")


diff_ab <- function(physeq, eval, sets = NULL, method, thresh, return, ...){
    method <- match.arg(method, c("cbea", "corncob", "deseq2"))
    return <- match.arg(return, c("pval", "sig"))
    eval <- match.arg(eval, c("fdr", "rset"))
    if (!eval %in% c("fdr")){
        if (is.null(sets)){
            stop("Sets are required for this analysis")
        }
    }
    if (method %in% c("corncob", "deseq2")){
        physeq <- speedyseq::tax_glom(physeq, taxrank = "GENUS")
        if (method == "corncob"){
            
        }
            
    } else if (method == "cbea"){
        sets <- const_set_taxtable(tax_table(physeq), rank = "GENUS")
        mat <- otu_table(physeq) %>% as("matrix") %>% t()
        scores <- cbea(obj = mat, set = sets, ...)
        
    }
}


const_set_taxtable <- function(table, rank) {
    set <- NULL # R CMD NOTE
    id <- which(colnames(table) == rank)
    table <- as.data.frame(as(table, "matrix")[, seq_len(id)])
    # Get all names -- have to do this because sometimes
    # different phyla might have the same
    # genus or family.
    all_names <- apply(table, 1, function(i) {
        paste(i, sep = ";_;", collapse = ";_;")
    })
    # Getting all unique names by removing NAs
    unq_names <- stats::na.omit(unique(all_names))
    unq_names <- unq_names[!stringr::str_ends(unq_names, "NA")]
    # Get all set ids
    sets <- purrr::map(unq_names, ~ {
        names(all_names)[which(all_names == .x)]
    })
    names(sets) <- unq_names
    # Constructing sets
    sets <- BiocSet::BiocSet(sets)
    set_sizes <- BiocSet::es_elementset(sets) %>%
        dplyr::count(set, name = "size")
    # filter set by set sizes
    sets <- BiocSet::left_join_set(sets, set_sizes)
    return(sets)
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