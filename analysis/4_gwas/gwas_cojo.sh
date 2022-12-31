#### submit_job.sh START ####
#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.cojo
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=6:00:00,h_data=10G
## Modify the parallel environment
## and the number of cores as needed:
#$ -pe shared 1

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load the job environment:
. /u/local/Modules/default/init/modules.sh
## Edit the line below as needed:
module load gcc

${software}/gcta_1.93.2beta/gcta64 \
--bfile ${bfile} \
--cojo-file ${gwas} \
--cojo-slct \
--cojo-actual-geno \
--cojo-p ${pthresh} \
--out ${out}

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
#### submit_job.sh STOP ####
