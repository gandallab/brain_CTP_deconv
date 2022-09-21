#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=1:00:00,h_data=8G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 1
# Email address to notify
# Notify when
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc/4.9.3
module load R/4.0.2

## substitute the command to run your code
## in the two lines below:
echo '/usr/bin/time -v hostname'
/usr/bin/time -v hostname

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "

#==============================================================================

Rscript ${scripts}/convert_cgn_gr.R \
${cgn_dir} \
${filen} \
${out_dir}