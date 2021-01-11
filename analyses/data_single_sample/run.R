library(targets)
library(tarchetypes)
library(tidyverse)

here::i_am("run.R")
setwd(file.path("D:/research/teailr/analyses/data_single_sample/")) 

targets::tar_make()
# targets::tar_visnetwork()
# targets::tar_make_future()
