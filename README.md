# Taxonomic enrichment analysis using competitive isometric log-ratios (cILR)

Method to perform enrichment analysis using isometric log ratio transformation for microbiome relative abundance data.  

Paper in preparation by Quang Nguyen.  

There are files with different use cases  
`*.pbs` files are submission files for the HPC Torque system at Dartmouth to run all analysis plans.  
`run.R` is the main interfacing R script for all `*.pbs` scripts  
`renv.lock` is the lockfile for `renv`  

There are various directories with different use cases  

`figures`: All the generated figures  
`R`: Contains all primary global cILR functions, data processing and plotting scripts  
`data`: Contains all the data sets available for usage in the manuscript  


The special folder is the `analyses` folder. Here we use the `targets` package to coordinate workflows for each piece of the analysis of the manuscript. 

* `data_` are workflows for real data analyses for the respective tasks.    
* `simulations_` are workflows for simulation analyses for the respective tasks.  

There are total of three tasks: `diff_ab` (Differential abundance), `single_sample` (single sample enrichment testing), and `pred`/`prediction` (predictive tasks).  For `single_sample` and `prediction`, there are mini-workflows that breaks up the big workflow into more managable chunks for running. For `single_sample`, this corresponds to the evaluation metrics. For `prediction`, this corresponds to the two predictive task type (classification and regression). As such, they have more folders and a master `*_functions` folders that hold functionality shared.   

In each of the sub-folders, there is a `_targets.R` file that contains all the rules. There is a `functions/` folder with any relevant functions, `output` with workflow outputs, and `data` for any relevant data sets used (copied from the master `data` folder). To run each workflow, navigate to the folder of interest (e.g. `setwd("analyses/data_single_sample")`) and execute the command `targets::tar_make()`/`targest::tar_make_future()`.  

For more information on the `targets` workflow, please refer to the [documentation](https://books.ropensci.org/targets/).   
The easiest way to run is to loop through all the directories using the positional arguments in `run.R`.  

