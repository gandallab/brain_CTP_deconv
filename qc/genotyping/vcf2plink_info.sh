#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.vcf2plink_info
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=8:00:00,h_data=16G
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
module load bcftools
module load plink/1.90b3.45

#==============================================================================

# Arguments
# - impvcf (data) directory
# - imp (output) directory
# - name
# - reference directory for rsid liftover script

# Convert VCF to plink
plink \
--vcf ${impvcf}/${impvcfn}.vcf.gz \
--double-id \
--keep-allele-order \
--make-bed \
--out ${impdir}/${impvcfn} > ${impdir}/${impvcfn}".log"

# Create info file for the variants
bcftools \
query -f '%CHROM\t%ID\t%FILTER\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/MAF\t%INFO/ER2\t%INFO/R2\n' \
${impvcf}/${impvcfn}.vcf.gz > ${impdir}/${impvcfn}.info

#==============================================================================

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
