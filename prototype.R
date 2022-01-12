library(phyloseq)
library(tidyverse)
library(CBEA)
library(BiocSet)
source("R/functions_data_ss.R")
source("R/cilr.R")

# Supporting functions #####
physeq2X <- function(otu_table){
    X <- otu_table(otu_table)
    X <- as(X, "matrix")
    X <- t(X)
    return(X)
}

set2A <- function(sets, otu_table){
    s_lab <- sets %>% es_set() %>% pull(set)
    e_lab <- taxa_names(otu_table)
    A <- matrix(0,nrow = length(e_lab), ncol = length(s_lab))
    rownames(A) <- e_lab
    colnames(A) <- s_lab
    for (i in s_lab){
        el_names <- sets %>% es_elementset() %>% filter(set == i) %>% 
            pull(element)
        el_index <- which(rownames(A) %in% el_names)
        A[,i][el_index] <- 1
    }
    return(A)
}


set_sizes <- c(20,50,100,150,200)
eval_settings <- cross_df(
    list(
        models = c("cbea"),
        distr = c("norm", "mnorm"),
        adj = c(TRUE, FALSE),
        iter = 1:10
    )
)

do_old_code_randset <- function(otu_table, eval_settings){
    X <- physeq2X(otu_table)
    new_anal_old <- eval_settings %>% 
        mutate(rand_sets = list(get_rand_sets(otu_table, size = 100, n_sets = 1))) %>% 
        mutate(res = pmap(., function(models, distr, adj, iter, rand_sets){
            A <- set2A(rand_sets, otu_table)
            cilr(X = X, A = A, resample = TRUE, output = "sig", distr = distr, 
                 adj = adj, preprocess = FALSE)
        }))
    new_anal_old <- new_anal_old %>% mutate(eval = map_dbl(res, function(x){
        pred <- x %>% pull(Set1) %>% factor(levels = c(0,1))
        label <- factor(rep(0, length(pred)), levels = c(0,1))
        1 - yardstick::specificity_vec(truth = label, estimate = pred, event_level = "second")
    }))
    out <- new_anal_old %>% group_by(distr, adj) %>% summarize(mean = mean(eval)) 
    return(out)
}


do_new_code_randset <- function(otu_table, settings){
    new_anal <- settings %>% 
        mutate(sets = list(get_rand_sets(otu_table, size = 100, n_sets = 1))) %>%
        mutate(res = pmap(., function(models, distr, adj, iter, sets){
            cbea(obj = otu_table, set = sets, output = "sig", distr = distr, adj = adj, n_boot = 1)
        })) %>% 
        mutate(eval = map_dbl(res, function(x){
            pred <- x %>% pull(Set1) %>% factor(levels = c(0,1))
            label <- factor(rep(0, length(pred)), levels = c(0,1))
            # sum(aer == 1 & label == 0)/sum(label == 0)
            #1 - yardstick::specificity_vec(truth = label, 
            #                               estimate = pred, 
            #                               event_level = "second")
            sum(pred == 1 & label == 0)/sum(label == 0)
        }))
    #new_anal <- new_anal %>% mutate(eval_2 = map_dbl(res, function(x){
    #    pred <- x %>% pull(Set1) %>% factor(levels = c(0,1))
    #    label <- factor(rep(0, length(pred)), levels = c(0,1))
    #    return(1 - yardstick::spec_vec(truth = label, estimate = pred, event_level = "second"))
    #}))
    
    out <- new_anal %>% group_by(distr, adj) %>% summarize(mean = mean(eval))
    return(out)
}

# DATA PROC ####

df <- readRDS(file = "data/hmp_supergingival_supragingival_16S.rds")
df <- gingival_processing(data = df)
otu <- df$physeq
sets <- df$sets
otu <- transform_sample_counts(otu, function(x) x + 1)
otu <- transform_sample_counts(otu, function(x) x/sum(x))

df_hard <- readRDS(file = "data/ackerman_ibd_16S.rds")
df_hard <- transform_sample_counts(df_hard, function(x) x + 1)
df_hard <- transform_sample_counts(df_hard, function(x) x/sum(x))

join_df <- sample_data(otu) %>% as(.,"data.frame") %>% 
    rownames_to_column("sample_id") %>% 
    dplyr::select(sample_id, HMP_BODY_SUBSITE) %>% 
    mutate(label = if_else(HMP_BODY_SUBSITE == "Supragingival Plaque", 1, 0))

# NEW ANALYSIS WITH OLD AND NEW CODE ####
set.seed(1020)
do_new_code_randset(otu_table = otu, settings = eval_settings)






# REDO NEW CODE WITH NEW ANALYSIS ####  
# OLD CODE WITH NEW ANALYSIS ####
eval_settings %>% 
    mutate(rand_set = list(get_rand_sets(otu, size = 100, n_sets = 1))) %>%
    mutate(model = pmap(., function(models, distr, adj, iter, rand_set){
        cilr(X = , A = , resample = TRUE, output = "sig", 
             adj = adj, distr = distr)
    }))

# NEW CODE WITH OLD ANALYSIS ####  
# eval_settings <- bind_rows(eval_settings, 
#                            tibble(models = "wilcoxon", 
#                                   distr = NA_character_, adj = NA))
eval_settings <- cross_df(
    list(
        models = c("cbea"),
        distr = c("norm", "mnorm"),
        adj = c(TRUE, FALSE)
    )
)



eval_settings %>% mutate(res = pmap(., function(models, distr, adj){
    print(c(models, distr, adj))
    scores <- cbea(obj = otu, set = sets, output = "sig", 
                   distr = distr, adj = adj)
    return(scores)
})) %>% mutate(eval = map_dbl(res, function(x){
    match_df <- left_join(x, join_df) 
    aer <- match_df %>% pull(Aerobic)
    label <- match_df %>% pull(label)
    # sum(aer == 1 & label == 0)/sum(label == 0)
    1 - yardstick::specificity_vec(truth = as.factor(label), 
                           estimate = as.factor(aer), 
                           event_level = "second")
}))


# NEW CODE WITH DISTRIBUTION #### 
rand_set <- get_rand_sets(otu, size = 100, n_sets = 1)
true_scores <- cbea(obj = otu, set = rand_set, output = "raw", distr=NULL, adj=FALSE)
true_scores <- true_scores %>% rename(truth = Set1)

test_set <- get_rand_sets(otu, size = 100, n_sets = 100)
rand_scores <- cbea(obj = otu, set = test_set, output = "raw", distr = NULL, adj = FALSE)
plot_scores <- rand_scores %>% pivot_longer(cols = -c(sample_id))
distr <- fitdistrplus::fitdist(data = plot_scores$value, distr = "norm")
distr_mnorm <- mixtools::normalmixEM(x  = plot_scores$value)

plot_scores <- plot_scores %>% rename(permuted_values = value) %>% left_join(true_scores) %>% 
    dplyr::select(-name) %>% 
    mutate(fitted_norm = rnorm(nrow(.), 
                        mean = distr$estimate["mean"], 
                        sd = distr$estimate["sd"]), 
           fitted_mnorm = rnormmix(n = nrow(.), 
                            lambda = distr_mnorm$lambda, 
                            mu = distr_mnorm$mu, 
                            sigma = distr_mnorm$sigma)) %>% 
    pivot_longer(cols = -c(sample_id, truth))


ggplot(plot_scores) +
    geom_histogram(mapping = aes(x = truth), 
                   alpha = 0.7, fill = "gray40", bins = 40) + 
    geom_histogram(mapping = aes(x = value, 
                    color = name, fill = name), 
                 alpha = 0.4, bins = 40) + 
    scale_fill_npg() + 
    scale_color_npg() + 
    guides(color = "none", fill = "none") + 
    facet_wrap(~name, nrow = 3) + theme_bw()














