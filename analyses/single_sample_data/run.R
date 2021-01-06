library(targets)
library(tarchetypes)
library(tidyverse)

here::i_am("run.R")
setwd(file.path("D:/research/teailr/analyses/single_sample_data/")) 

targets::tar_make()
# targets::tar_make_future()
