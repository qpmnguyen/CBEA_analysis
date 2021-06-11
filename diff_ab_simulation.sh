#!/bin/bash -l

# declare the name for this job 

#SBATCH --job-name=DIFF_AB_SIMULATION

# Distributing jobs across 5 nodes with 20 cores each node 
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --hostlist=k[40-58]

# Requesting 4000 MB of RAM
#SBATCH --mem-per-cpu=4000

# request 60 hours of wall time
#SBATCH --time=60:00:00

# Dispatching job to standard partitions on discovery
#SBATCH --partition=standard

# specify your email address and when things are emailed

#SBATCH --mail-user=Quang.P.Nguyen.GR@dartmouth.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# By default, SLURM scripts execute in your home directory, not the
# directory from which they were submitted.

cd $SLURM_SUBMIT_DIR

# Run run.R as a script to start the targets pipeline
conda activate teailr
Rscript run.R --ncores 100 --dir "analyses/simulations_diff_ab/" --remove FALSE --parallel TRUE

