library(tidymodels)
library(tidyverse)
library(phyloseq)

source("../../R/cilr.R")
source("../../R/utils.R")
source("functions/step_functions.R")
physeq <- readRDS(file = "../../data/nielsen_ibd_wgs.rds")

process_wgs <- function(physeq, lab_col, case_label){
    label <- sample_data(physeq) %>% as("data.frame") %>% dplyr::select(!!lab_col)
    label <- ifelse(label == case_label, 1, 0) %>% as.factor()
    data <- otu_table(physeq) %>% as("matrix") %>% as.data.frame() %>% 
        dplyr::select(starts_with("s__")) 
    data <- dplyr::bind_rows(label, data)
    table <- tax_table(physeq) %>% unclass() %>% as.data.frame() %>% 
        rownames_to_column(var = "colnames") %>% filter(str_detect(colnames, "s__")) %>%
        column_to_rownames(var = "colnames")
    return(list(data = data, table = table))
}


