library(drake)
library(here)
library(glue)
library(tidyverse)

here::i_am("analyses/single_sample_data/plans.R")

source(here("analyses/single_sample_data/functions/enrichment.R"))


set.seed(1234)

plan_enrichment <- drake_plan(
    # get the data
    data = readRDS(file = here("data/hmp_supergingival_supragingival_16S.rds")) %>% enrichment_processing(),
    # extract the label 
    label = data$label,
    # get score and evaluate  
    results = 
    models = target({
        cross(list(
            resample = TRUE,
            output = c("cdf", "zscore"),
            distr = c("norm", "mnorm"),
            adj = c(TRUE, FALSE),
            maxrestarts=1000, 
            epsilon = 1e-06, 
            maxit= 1e5))
    }),
    score = target({
        arg <- rlist::list.append(models, X = data$X, A = data$A)
        pmap(arg, .f = cilr)
    }, dynamic = map(models))
  # cilr_score = target({
  #     cilr(data$X, data$A, output = output, distr = distr, adj = adj, resample = T,  
  #          maxrestarts=1000, epsilon = 1e-06, maxit= 1e5)
  #   },
  #   dynamic = map(models)
  # ),
  # other_score = target({
  #       generate_alt_scores(data$X, data$A, method = method, preprocess = T, pcount = 1) 
  #   },
  #   dynamic = cross(
  #     method = c("ssgsea", "gsva")
  #   )
  # ),
  # evaluate_cilr = target(
  #   criterion(scores= cilr_score, results = label),
  #   dynamic = map(cilr_score)
  # ), 
  # evaluate_other = target(
  #   criterion(scores = other_score, results = label),
  #   dynamic = map(other_score)
  # ),
  # summary = target(
  #   rbind(evaluate_cilr, evaluate_other),
  #   dynamic = combine(evaluate_cilr, evaluate_other)
  # ),
  # save = target({
  #     saveRDS(file = here("analyses/single_sample_data/results/enrichment.rds"))
  # })
)
drake::vis_drake_graph(plan_enrichment)
drake::make(plan_enrichment)
