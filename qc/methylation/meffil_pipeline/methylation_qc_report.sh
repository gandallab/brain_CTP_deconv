#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.meffil_qc
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=4:00:00,h_data=8G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 8
# Email address to notify
#$ -M chloe.yap@uq.net.au
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
module load R/4.0.2

#==============================================================================

# Arguments
# - Rscript directory for methylation_qc_report.R
# - qcdir: directory for QC outputs
# - samplesheet directory (assumed to be within qcdir); filen: study identifier (eg. "seed")
# - geno: directory to genotyping bfiles corresponding to this sample
# - number of cores for R to use (N.B. ensure ncpus argument above matches)

## substitute the command to run your code
## in the two lines below:
Rscript ${scripts}/methylation_qc_report.R \
${qcdir} \
"samplesheet_"${filen}".csv" \
${filen} \
${geno} \
8

#==============================================================================

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "


