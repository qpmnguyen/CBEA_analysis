#SBATCH --job-name=DIFF_AB_SIMULATION

# Distributing jobs across 5 nodes with 20 cores each node 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

# Requesting RAM in MBs
#SBATCH --mem-per-cpu=20000

# requesting walltime 
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
Rscript run.R --ncores 10 --dir "analyses/simulations_diff_ab/" --remove FALSE --parallel TRUE --cluster TRUE