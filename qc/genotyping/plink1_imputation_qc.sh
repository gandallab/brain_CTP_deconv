#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.imputation_qc
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=6:00:00,h_data=60G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 1
# Email address to notify
#$ -M $USER
# Notify when
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
#module load gcc/4.9.3
module load plink

#==============================================================================

## substitute the command to run your code
## in the two lines below:

# Imputation QC filters
plink \
--bfile ${bfile}  \
--chr 1-22 \
--hwe 0.000001 \
--maf 0.01 \
--exclude ${exclude}  \
--make-bed \
--threads 1 \
--memory 60000 \
--out ${out}

#==============================================================================

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "



