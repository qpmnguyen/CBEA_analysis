library(drake)
library(here)
library(glue)

here::i_am("analyses/single_sample_data/plans.R")

source(here("analyses/single_sample_data/R/functions.R"))

plan <- drake_plan(
  raw_data = readRDS(file = here("data/hmp_supergingival_supragingival_16S.rds")),
  proc_data = process_data(raw_data),
  cilr_score = target(
    cilr(proc_data$X, proc_data$A, output = "cdf", distr = distr, adj = adj, resample = T,  
         maxrestarts=1000, epsilon = 1e-06, maxit= 1e5),
    transform = cross(
      adj = c(TRUE, FALSE),
      distr = c("norm", "mnorm")
    )
  ),
  other_score = target(
    generate_alt_scores(proc_data$X, proc_data$A, method = method, preprocess = T, pcount = 1),
    transform = cross(
      method = c("ssgsea", "gsva")
    )
  ),
  evaluate_cilr = target(
    eval(label = proc_data$label, cilr_score$Aerobic),
    transform = map(cilr_score)
  ), 
  evaluate_other = target(
    eval(label = proc_data$label, other_score$Aerobic),
    transform = map(other_score)
  ),
  summary = target(
    rbind(evaluate_cilr, evaluate_other),
    transform = combine(evaluate_cilr, evaluate_other)
  )
)
drake::vis_drake_graph(plan)
