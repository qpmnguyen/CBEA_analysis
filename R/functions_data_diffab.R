library(DESeq2)
library(corncob)
library(glue)
library(binom)
library(CBEA)
requireNamespace("speedyseq", quietly = TRUE)

# retrieving 
source("R/functions_data_ss.R")
source("R/utils.R")


diff_ab <- function(obj, eval, abund_values = "16SrRNA", sets = NULL, method, thresh, return, ...){
    method <- match.arg(method, c("cbea", "corncob", "deseq2"))
    return <- match.arg(return, c("pval", "sig"))
    eval <- match.arg(eval, c("fdr", "rset"))
    
    physeq <- mia::makePhyloseqFromTreeSummarizedExperiment(obj, abund_values = abund_values)
    
    if (!eval %in% c("fdr")){
        if (is.null(sets)){
            stop("Sets are required for this analysis")
        }
    }
    if (method %in% c("corncob", "deseq2")){
        physeq <- transform_sample_counts(physeq, function(x) x + 1)
        if (eval == "fdr"){
            physeq <- speedyseq::tax_glom(physeq, taxrank = "GENUS")
        } else if (eval == "rset"){
            grp_vec <- rep(NA, ntaxa(physeq))
            names(grp_vec) <- taxa_names(physeq)
            set_names <- sets %>% es_set() %>% pull(set)
            for (i in seq_along(set_names)){
                el_names <- sets %>% es_elementset() %>% filter(set == set_names[i]) %>% 
                    pull(element)
                idx <- which(names(grp_vec) %in% el_names)
                grp_vec[idx] <- set_names[i]
            }
            grp_vec <- as.vector(grp_vec)
            physeq <- speedyseq::merge_taxa_vec(physeq, grp_vec)
        }
        if (method == "corncob"){
            args <- list(formula = ~ group, phi.formula = ~ group, formula_null = ~1, 
                        phi.formula_null = ~ group, test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
            args$data <- physeq
            mod <- do.call(corncob::differentialTest, args)
            sig <- mod$p
            names(sig) <- as.vector(tax_table(physeq)[names(sig), "GENUS"])
        } else if (method == "deseq2"){
            args <- list(test = "LRT", fitType = "local", reduced = ~1)
            deseq <- phyloseq_to_deseq2(physeq, ~factor(group))
            args$object <- deseq
            mod <- do.call(DESeq2::DESeq, args)
            res <- DESeq2::results(mod)
            sig <- res$pvalue
            names(sig) <- as.vector(tax_table(physeq)[rownames(res), "GENUS"])
        }
    } else if (method == "cbea"){
        if (eval == "fdr"){
            sets <- const_set_taxtable(tax_table(physeq), rank = "GENUS")
        } 
        assay(obj, abund_values) <- assay(obj, abund_values) + 1
        obj <- mia::transformCounts(obj, 
                                       abund_values = abund_values, 
                                       method = "relabundance", name = "main_input")
        args <- list(...)
        args <- c(args, list(
            obj = obj,
            abund_values = "main_input", 
            set = sets
        ))
        # ensure that the large composition will always be used. 
        if ("distr" %in% names(args)){
            if (args$distr == "mnorm"){
                args <- c(args, list(
                    control = list(fix_comp = "large")
                ))
            }
        }
        mod <- do.call(cbea, args)
        scores <- tidy(mod) %>% add_column(sample_ids = colnames(assay(obj, "main_input")), .before = 1)
        group_var <- colData(obj)$group
        sig <- map_dbl(2:ncol(scores), ~{
            values <- scores %>% pull(.x) 
            t.test(x = values[which(group_var == 0)], y = values[which(group_var == 1)])$p.value
        })
        genera_names <- colnames(scores)[-1] %>% strsplit(";_;") %>% map_chr(., ~.x[[length(.x)]])
        names(sig) <- genera_names
    } else if (method == "gsva") {
        sig <- NULL
    }
    
    if (return == "sig"){
        sig <- ifelse(sig <= thresh, 1, 0)
    }
    return(sig)
}

eval_results <- function(obj, set, sig_vec){
    # first, let's extract names that should be differentially abundant 
    physeq <- mia::makePhyloseqFromTreeSummarizedExperiment(obj, "16SrRNA")
    
    target_elements <- set %>% es_elementset() %>% 
        filter(set %in% c("Aerobic", "Anaerobic")) %>%
        pull(element)
    
    g_sets <- const_set_taxtable(table = tax_table(physeq), "GENUS")
    
    target_diffab <- g_sets %>% es_elementset() %>% 
        filter(element %in% target_elements) %>% pull(set) %>% 
        unique() %>% strsplit(";_;") %>% map_chr(., ~.x[length(.x)])
    
    # second, we extract names that are actually differentially abundant 
    real_diffab <- names(sig_vec)[p.adjust(sig_vec) <= 0.1]
    
    # third, we examine the intersection and 
    inter <- length(intersect(target_diffab, real_diffab))
    
    # fourth, construct the binomial test 
    bintest <- binom.confint(x = inter, n = length(sig_vec), method = "ac")
    return(
        list(
            estimate = bintest$mean, 
            lower = bintest$lower,
            upper = bintest$upper
        )
    )
    
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