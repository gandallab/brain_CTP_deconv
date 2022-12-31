#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.vcf2plink_rsidliftover
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=10:00:00,h_data=10G
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
module load bcftools
module load plink
module load R/4.0.2

#==============================================================================

# Arguments
# - impvcf (data) directory
# - imp (output) directory
# - chromosome
# - reference directory for rsid liftover script

plink \
--vcf ${impvcf}/"chr"${chr}".dose.vcf.gz" \
--id-delim \
--keep-allele-order \
--make-bed \
--out ${impdir}/imputed_chr${chr} > ${impdir}/"imputed_chr"${chr}".log"

bcftools \
query -f '%CHROM\t%ID\t%FILTER\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/RefPanelAF\t%INFO/INFO\n' ${impvcf}/"chr"${chr}".dose.vcf.gz" > ${impdir}/"imputed_chr"${chr}".info"

# .bim files
Rscript ${ref}/rsid_liftover_aut_x.R \
${impdir}/"imputed_chr"${chr}".bim"

# .info files
cut -d ' ' -f2 ${impdir}/"imputed_chr"${chr}".bim" > ${impdir}/tmp_${chr}
paste -d ' ' ${impdir}/tmp_${chr} ${impdir}/"imputed_chr"${chr}".info" > ${impdir}/foo_${chr}
awk '{print $2, $1, $4, $5, $6, $7, $8, $9, $10, $11}' ${impdir}/foo_${chr} > ${impdir}/"imputed_chr"${chr}".info"
rm ${impdir}/tmp_${chr}
rm ${impdir}/foo_${chr}


#==============================================================================

# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "
