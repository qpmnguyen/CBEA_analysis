#!/usr/bin/bash 
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate cbea

Rscript R/runtime.R
