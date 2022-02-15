# CBEA: Competitive Balances for Taxonomic Enrichment Analysis. 

[![DOI](https://zenodo.org/badge/253318333.svg)](https://zenodo.org/badge/latestdoi/253318333)


Method (formerly known as cILR) to perform enrichment analysis using competitive balances for microbiome relative abundance data.     

Manuscript under review at PLoS Computational Biology. Pre-print [available on bioRxiv](https://www.biorxiv.org/content/10.1101/2021.09.07.459294v1.full). 

There are files with different use cases  
`*.pbs` files are submission files for the HPC Torque system at Dartmouth to run all analysis plans.  
`run.R` is the main interfacing R script for all `*.pbs` scripts  
`renv.lock` is the lockfile for `renv`  

There are various directories with different use cases  
`figures`: All the generated figures  
`R`: Contains all primary global cILR functions, data processing and plotting scripts  
`data`: Contains all the data sets available for usage in the manuscript  


All analyses are ran as projects using the `targets` package. Each sub-analysis is a mini project nested using the projects feature. For more documentation for using projects with `targets`, please refer to the corresponding documentation [online](https://books.ropensci.org/targets/). More specifically, each project has a name, a workflow instruction file (prefix `script_*.R`), and a data store. The file `targets.yaml` contains all information on sub-analyses available. To run each sub analysis, use the `run_new.R` file in script mode (e.g. `RScript run_new.R --project "sim_diffab"`).   
