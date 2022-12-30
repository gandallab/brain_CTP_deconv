#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.gwas_linear
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=5:00:00,h_data=10G
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

${software}/gcta_1.93.2beta/gcta64 \
--fastGWA-mlm \
--bfile ${bfile} \
--pheno ${pheno} \
--covar ${covar} \
--qcovar ${qcovar} \
--keep ${keep} \
--reml-alg 1 \
--out ${out} \
--thread-num 1

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####
