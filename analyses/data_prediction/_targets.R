library(targets)
library(tarchetypes)
library(tidyverse)
source("functions/pred_functions.R")

files <- tar_files(files, c("../../data/nielsen_ibd_wgs.rds", 
                            "../../data/hmp_stool_tongue_16S.rds",
                            "../../data/hmp_stool_tongue_wgs.rds"))
data <- tar_target(data, {
    dat <- readRDS(files)
    if (str_detect(files, "wgs")){
        if (str_detect(files, "ibd")){
            dat <- process_pred(dat, lab_col = "disease", "IBD", "wgs") 
        } else {
            dat <- process_pred(dat, lab_col = "body_site", "stool", "wgs")
        }
    } else if (str_detect(files,"16S")){
        dat <- process_pred(dat, lab_col = "HMP_BODY_SITE", "Gastrointestinal Tract", "16S")
    }
    dset_names <- str_split(files, "/")[[1]] %>% 
        tail(n = 1) %>% 
        strsplit(".", fixed = T) %>% .[[1]] %>% 
        head(n = 1)
    dat <- list(dat = dat, title = dset_names)
    dat
}, pattern = map(files))





data_ibd <- tar_rds(data_ibd, {
    readRDS(file = "../../data/nielsen_ibd_wgs.rds")
})

ibd_proc <- tar_target(ibd_proc, {
    process_wgs(data_ibd, lab_col = "disease", "IBD")
})


add_methods <- tibble(methods = c("gsva", "clr", "ssgsea"))
cilr_settings <- cross_df(list(
    method = c("cilr"),
    output = c("zscore", "cdf"),
    distr = c("mnorm", "norm"),
    adj = c(TRUE, FALSE)
))


cilr_agg <- tar_map(unlist = FALSE, values = cilr_settings,
    tar_target(model, {
        if ("Genus" %in% colnames(data$dat$table)){
            lvl <- "Genus"
        } else {
            lvl <- "GENUS"
        }
        agg <- generate_aggregation(data$dat, level = lvl, method = method, output = output, 
                                         distr = distr, adj = adj)
        list(agg = agg, title = data$title)
    }, pattern = map(data)),
    tar_target(eval, {
        metrics <- fit_and_eval(model$agg, nfolds = 10, task = "classification")[2,]
        tibble(dset = model$title, method = method, output = output, distr = distr, adj = adj, 
               auc = metrics$mean, upper = metrics$mean + metrics$std_err, 
               lower = metrics$mean - metrics$std_err)
    }, pattern = map(model))
)

other_agg <- tar_map(unlist = FALSE, values = add_methods, 
    tar_target(model, {
        if ("Genus" %in% colnames(data$dat$table)){
            lvl <- "Genus"
        } else {
            lvl <- "GENUS"
        }
        agg <- generate_aggregation(data$dat, level = lvl, method = methods)
        list(agg = agg, title = data$title)
    }, pattern = map(data)),
    tar_target(eval, {
        metrics <- fit_and_eval(model$agg, nfolds = 10, task = "classification")[2,]
        tibble(dset = model$title, method = methods, auc = metrics$mean, upper = metrics$mean + metrics$std_err, 
               lower = metrics$mean - metrics$std_err)
    }, pattern = map(model))
)


combine <- tar_combine(combine, other_agg[[2]], cilr_agg[[2]], 
                       command = dplyr::bind_rows(!!!.x))


list(files, data, cilr_agg, other_agg, combine)

