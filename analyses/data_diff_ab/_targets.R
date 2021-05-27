library(targets)
library(tarchetypes)
library(tidyverse)
library(phyloseq)
library(stringr)
library(future)

plan(multisession)
source("functions/diff_ab_functions.R")
set.seed(1020)

tar_option_set(error = "workspace")

# First, define the file paths and the file list 
fdr_files <- tibble(
    path = c("../../data/hmp_stool_16S.rds",
             "../../data/hmp_stool_wgs.rds"),
    dset = c("stool_16S", "stool_wgs"),
    type = c("16S", "WGS"),
    agg = c("GENUS", "Genus")
)

# Same thing for pwr files  
pwr_files <- tibble(
    path = c("../../data/hmp_supergingival_supragingival_16S.rds"),
    dset = c("gingival_16S"),
    type = c("16S"),
    agg = c("GENUS")
)

# Define evaluation grid which includes all the cILR variants  
eval_grid <- tar_target(eval_grid,{ 
    df <- cross_df(list(
        methods = c("cilr_wilcox", "cilr_welch"),
        distr = c("norm", "mnorm"),
        output = c("cdf", "zscore"),
        adj = c(TRUE, FALSE)))
    add_methods <- tibble(methods = c("deseq2", "corncob"))
    df <- dplyr::bind_rows(df, add_methods)
    df
})

# Define the analysis for fdr, which is mapped over files but repeated 1000 times  
# Since tar_rep does not support pattern = map(target_names), having to perform a loop within each data set
# Which actually doesn't take that much time 
fdr_analysis <- tar_map(unlist = FALSE, values = fdr_files, names = "dset",
    tar_rep(fdr_rep, {
        print(eval_grid)
        data <- readRDS(path)
        group <- rbinom(n = nsamples(data), size = 1, prob = 0.5)
        sample_data(data)$group <- as.factor(group)
        eval <- vector(length = nrow(eval_grid))
        for (i in seq_len(nrow(eval_grid))){
            result <- diff_ab(data, agg_level= agg, data_type = type,
                              method = eval_grid$methods[i], 
                              padj = FALSE, 
                              return = "sig",
                              output = eval_grid$output[i],
                              distr = eval_grid$distr[i], adj = eval_grid$adj[i])
            eval[i] <- eval_function(result, ci = FALSE)
        }
        eval_grid %>% mutate(eval = eval)
    }, batches = 50, reps = 10),
    tar_rds(save_file_fdr, {
        saveRDS(fdr_rep, file = glue("output/{dset}_fdr.rds", dset = dset))
    })
)

# Power analysis maps across eval_grid but there is only one data set and no repetition is required 
pwr_analysis <- tar_target(pwr, {
    # load_files 
    print("Load files")
    data <- readRDS(pwr_files$path)
    physeq <- data$physeq
    annotation <- data$annotation
    # define group 
    print("Define Group")
    group <- ifelse(sample_data(physeq)$HMP_BODY_SUBSITE == "Supragingival Plaque", 1, 0)
    sample_data(physeq)$group <- group %>% as.factor()
    
    # perform analysis 
    print("Perform Analysis")
    print(paste0("Method: ", eval_grid$methods))
    print(paste0("Distribution: ", eval_grid$distr))
    print(paste0("Adjustment: ", eval_grid$adj))
    print(paste0("Output: ", eval_grid$output))
    
    result <- diff_ab(physeq, agg_level = "GENUS", data_type = "16S", return = "sig", 
                      method = eval_grid$methods, 
                      distr = eval_grid$distr, 
                      output = eval_grid$output, 
                      adj = eval_grid$adj)
    # restrict annotation to Aerobic or Anaerobic only and consider all of them as labelled 
    print("Retrieving annotations")
    annotation <- annotation %>% filter(Meth %in% c("Aerobic", "Anaerobic"))
    print("Length of result object...")
    print(length(result))
    result_df <- tibble(Genera = names(result), values = unname(result))
    annotated_result <- inner_join(annotation, result_df, by = "Genera")
    eval <- eval_function(annotated_result$values, ci = TRUE)
    print("Finializing evaluation by binding evaluation with grid")
    tibble(eval_grid, eval)
}, pattern = map(eval_grid))

# saving files  
pwr_save_file <- tar_rds(pwr_save_file, 
                         saveRDS(pwr, file = glue("output/{dset}_pwr.rds", dset = pwr_files$dset)))

list(eval_grid, fdr_analysis, pwr_analysis, pwr_save_file)
