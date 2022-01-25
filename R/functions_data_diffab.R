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
        physeq <- speedyseq::tax_glom(physeq, taxrank = "GENUS")
        if (method == "corncob"){
            args <- list(formula = ~ group, phi.formula = ~ group, formula_null = ~1, 
                        phi.formula_null = ~ group, test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
            args$data <- physeq
            mod <- do.call(corncob::differentialTest, args)
            sig <- mod$p
            names(sig) <- as.vector(tax_table(physeq)[names(sig), "GENUS"])
        } else if (method == "DESeq2"){
            args <- list(test = "LRT", fitType = "local", reduced = ~1)
            deseq <- phyloseq_to_deseq2(physeq, ~factor(group))
            args$object <- deseq
            mod <- do.call(DESeq2::DESeq, args)
            res <- DESeq2::results(mod)
            sig <- res$pvalue
            names(sig) <- as.vector(tax_table(physeq)[rownames(res), "GENUS"])
        }
    } else if (method == "cbea"){
        sets <- const_set_taxtable(tax_table(physeq), rank = "GENUS")
        assay(obj, abund_values) <- assay(obj, abund_values) + 1
        obj <- mia::transformCounts(obj, 
                                       abund_values = abund_values, 
                                       method = "relabundance", name = "main_input")
        mia::tra
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
                    control = list(fit_comp = "large")
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