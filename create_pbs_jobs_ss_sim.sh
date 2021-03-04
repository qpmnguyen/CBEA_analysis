#!/bin/bash

# number of nodes per job (1 node = 16 cores = 64/128GB RAM)
NODES=1

# number of cores per job (1 core = 4/8GB RAM)
PPN=20

# email where to send job notifications ***MUST ESCAPE THE '@'*** e.g. fist.last\@dartmouth.edu
EMAIL="Quang.P.Nguyen.GR\@dartmouth.edu"

# desired wall time HH:MM:SS for each job
WALLTIME=20:00:00

# create user directory in scratch if it doesn't exist
if test -d /dartfs-hpc/scratch/$USER
then
	:
else
	mkdir /dartfs-hpc/scratch/$USER
fi

# reset jobs counter
count=0

# reads the last arguent passed on the command line -- must me '-x' if you want to submit jobs to the selecte queue, otherwise, it's just a test, no submission
lastarg=${@: -1}

declare -a dir=("analyses/simulations_single_sample_auc/" "analyses/simulations_single_sample_fdr/" "analyses/simulations_single_sample_pwr")
# == FOR LOOP
for i in "${dir[@]}" 
do
	# increase and prints the current job number
	((count+=1))
	echo "job $count"
	
	### creates a temporary PBS script in user's scratch space
	# PBS script file name
	jobname="/scratch/$USER/job_$count_${RANDOM}"
	
	jobname=${jobname//\ /\_}
	jobfilename="${jobname}.pbs"
	
	#echo "#PBS -j oe" >> ${jobfilename}
	
	# sets the requested number of nodes and cores
	echo "#PBS -l nodes=${NODES}:ppn=${PPN}" >> ${jobfilename}
	
	# sets requested walltime
	echo "#PBS -l walltime=${WALLTIME}" >> ${jobfilename}
	
	# run only on Intel nodes, only cellJ and other optional features
    #echo "#PBS -l feature=intel,cellJ,etc..." >> ${jobfilename}
	
	# sets which job events trigger e-mail to pre-defined e-mail address: (b)egin, (e)nd, (a)bort
	echo "#PBS -m bea" >> ${jobfilename}
	echo "#PBS -M $EMAIL" >> ${jobfilename}
	
	# sets the work directory
	echo "cd \$PBS_O_WORKDIR" >> ${jobfilename}
    echo "conda activate teailr" >> ${jobfilename}
	# build the command to be added to PBS script
	command="Rscript run.R --ncores 20 --dir $i"	
	
	# displays the command and adds it to the PBS script
	echo $command	
	echo $command >> ${jobfilename}

	# if the last argument on the command line was '-x', the PBS script is added to the queue
	if [ -n "$lastarg" ] && [ ${lastarg} == "-x" ]
	then
		echo "inserting in queue"
		qsub -N $jobname $jobfilename
	fi
	
	# displays PBS script path and content							
	#echo $jobname
	#cat $jobfilename
	
	# deletes PBS script (cleanup)
	rm -f $jobfilename
done
