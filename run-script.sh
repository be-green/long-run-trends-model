#!/bin/bash 
#SBATCH -n 24 #Request 4 tasks (cores)
#SBATCH -t 0-12:30 #Request runtime of 30 minutes
#SBATCH -p sched_mit_sloan_batch #Run on sched_engaging_default partition
#SBATCH --mem-per-cpu=8GB #Request 4G of memory per CPU
#SBATCH -o output_%j.txt #redirect output to output_JOBID.txt
#SBATCH -e error_%j.txt #redirect errors to error_JOBID.txt
#SBATCH --mail-type=BEGIN,END #Mail when job starts and ends
#SBATCH --mail-user=brice.c.green@gmail.com #email recipient
module load knitro/12.0.0
module load mit/matlab/2020a
cd ~/documents/long-run-trends-model
matlab -nodesktop -nosplash -r "optimize_model_multistart; exit;"