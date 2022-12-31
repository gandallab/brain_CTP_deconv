#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID.sbayesr
#$ -j y
## Edit the line below as needed:
#$ -l h_rt=16:00:00,h_data=50G
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
module load R/4.0.2
module load gcc

#==============================================================================
#
# PGS analysis (SBayesR, GCTBv2.0)
# http://cnsgenomics.com/software/gctb/#Bayesianalphabet
#
#==============================================================================

# Tutorial: http://cnsgenomics.com/software/gctb/#SummaryBayesianAlphabet

#==============================================================================
# Variables
#==============================================================================

echo Reading ldm from: ${mldm}
echo Reading sumstats: ${dir}
echo Taking the bfile ${bfile}
echo h2 = ${h2}
echo pi = "${pi}"
echo gamma = "${gamma}"
echo Will write output to: ${out}
echo Output name prefix: ${filen}
echo Target name: ${target}

#------------------------------------------------------------------------------
# 3. Run SBayesR
#------------------------------------------------------------------------------

${software}/gctb_2.03beta_Linux/gctb \
--sbayes R \
--mldm ${mldm} \
--gwas-summary ${dir} \
--hsq ${h2} \
--chain-length 50000 \
--burn-in 10000 \
--seed 12345 \
--no-mcmc-bin \
--pi "${pi}" \
--gamma "${gamma}" \
--exclude-mhc \
--out-freq 10 \
--out ${out}/${filen}

#------------------------------------------------------------------------------
# 4. Check SBayesR outputs
#------------------------------------------------------------------------------

Rscript ${scripts}/pgs_sbayesr_check.R \
${dir} \
${out}/${filen} \
${out} \
${filen}

#------------------------------------------------------------------------------
# 5. Generate PGS
#------------------------------------------------------------------------------

awk '{print $2, $5, $8 }' ${out}/${filen}.snpRes > ${out}/${filen}"_SBayesR_predictor.txt"

# PRS for the target individuals
${software}/plink2 \
--bfile ${bfile} \
--score ${out}/${filen}"_SBayesR_predictor.txt" \
--out ${out}/${target}_pgs.${filen}.sbayesr



