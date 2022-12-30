#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.meffil_qc_normalise
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=10:00:00,h_data=40G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 1
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
# - Rscript directory for methylation_qc_normalise.R
# - qcdir: directory for QC outputs
# - filen: study identifier (eg. "seed")
# - pc: number of PCs
# - number of cores for R to use (N.B. ensure ncpus argument above matches)

Rscript ${scripts}/methylation_qc_normalise.R \
${qcdir} \
${filen} \
${pc} \
1
#>${qcdir}/${filen}_normalise_1core.o 2>${qcdir}/${filen}_normalise_1core.e

#==============================================================================

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "

