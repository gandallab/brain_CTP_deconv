software=~/project-gandalm/software
ref=~/shared-gandalm/brain_CTP/Data/genotyping/ref
impvcf=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD_synapse
impdir=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD_synapse
dbsnp_hg19=~/shared-gandalm/refGenomes/hg19/dbSNP
grm=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD_synapse/grm
qcdir=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD_synapse/qc_check

module load plink/1.90b624
module load bcftools/1.11
module load htslib/1.12
module load vcftools/0.1.16


#==============================================================================
#
# IMPUTED DATA
#
#==============================================================================

# Concatenate chromsomes into one VCF file
for chromosome in $impvcf/UCLA_ASD.chr*.dose.vcf.gz; do
    bcftools index -t $impvcf/UCLA_ASD.vcf.gz
done
bcftools concat $impvcf/UCLA_ASD.chr*.dose.vcf.gz -o $impvcf/UCLA_ASD.vcf.gz -Oz
bcftools index -t $impvcf/UCLA_ASD.vcf.gz

# Convert VCF to plink and generate an info file on the variants
qsub -v ref=${ref},\
impvcf=${impvcf},\
impdir=${impdir},\
impvcfn=UCLA_ASD \
${ref}/vcf2plink_info.sh

# Find SNPs with INFO<0.3
awk '$10<0.3' ${impdir}/UCLA_ASD.info | awk '{print $2}' > ${impdir}/SNPs_aut_info_less0.3.txt

# Imputation QC
# HWE 1-e6, MAF 0.01, INFO 0.3
qsub -v software=${software},\
bfile=${impdir}/UCLA_ASD,\
exclude=${impdir}/SNPs_aut_info_less0.3.txt,\
out=${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6 \
${ref}/plink1_imputation_qc.sh

# Update rsID
${software}/plink2 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6 \
--make-bed \
--set-all-var-ids @:# \
--output-chr 26 \
--rm-dup force-first \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_tmp

# Before continuing, run hg19_chrpos_rsid.sh to generate ${dbsnp_hg19}/00-All_coord_rsid_format_unique4_changeorder_out.txt
${software}/plink2 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_tmp \
--make-bed \
--update-name ${dbsnp_hg19}/00-All_coord_rsid_format_unique4_changeorder_out.txt \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid

#==============================================================================
#
# QC METRICS AND CHECKS
#
#==============================================================================

#-------------------------------------
# Align to 1KG for pop comparison
#-------------------------------------

# Check RAF + find allele outliers + change ref
${software}/plink2 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid \
--freq \
--keep-allele-order \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid

module load R/4.0.2
Rscript \
${ref}/RAF_outlier_filter_1KG.R \
${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid.afreq \
${ref}/1000GP_Phase3_combined_rsid_clean.legend \
${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid

${software}/plink2 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid \
--ref-allele ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid.afreq2recode 2 1 \
--exclude ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid.afreqoutlier \
--set-hh-missing \
--make-bed \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG

#--------------------------------------
# PCA
#--------------------------------------
# 1. Pre-PCA
# a) Filter out the rare variants
${software}/plink2 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG \
--maf 0.05 \
--make-bed \
--out ${grm}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.05

# b) sort the order of the SNPs for picking the shared set between reference and target data
awk '{print $2}' ${grm}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.05.bim | sort > \
${grm}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.05.snps.txt
# sort ${REF}.05.SNPs.txt > ${REF}.05.snps_sort.txt # did this in R
comm -12 ${ref}/1000GP_Phase3_combined_rsid_clean.legend.05.snps_sort.txt ${grm}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.05.snps.txt > \
${grm}/common.1KG.UCLA_ASD.SNPs.txt

# c) get the common SNPs for the reference and target data
${software}/plink2 \
--bfile ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05 \
--extract ${grm}/common.1KG.UCLA_ASD.SNPs.txt \
--make-bed \
--out ${grm}/1000G_phase3_20130502_combined_snpsonly.05.common_UCLA_ASD

${software}/plink2 \
--bfile ${grm}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.05 \
--maf 0.05 \
--extract ${grm}/common.1KG.UCLA_ASD.SNPs.txt \
--make-bed \
--out ${grm}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.05.05.common

#--------------------------------------
# 2. PC projection

# a) generate GRM for 1000G
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${grm}/1000G_phase3_20130502_combined_snpsonly.05.common_UCLA_ASD \
--make-grm \
--out ${grm}/1000G_phase3_20130502_combined_snpsonly.05.common_UCLA_ASD

# b) Take first 2 principal components of reference (1000G)
${software}/gcta_1.93.2beta/gcta64 \
--grm ${grm}/1000G_phase3_20130502_combined_snpsonly.05.common_UCLA_ASD  \
--pca 2 \
--out ${grm}/1000G_phase3_20130502_combined_snpsonly.05.common_UCLA_ASD_pca2

# c) Generate loadings of each SNP in the reference on the above PCs
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${grm}/1000G_phase3_20130502_combined_snpsonly.05.common_UCLA_ASD \
--pc-loading ${grm}/1000G_phase3_20130502_combined_snpsonly.05.common_UCLA_ASD_pca2 \
--out ${grm}/1000G_phase3_20130502_combined_snpsonly.05.common_UCLA_ASD_pca2_snp_loading

# d) Project onto reference sample
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${grm}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.05.05.common \
--project-loading ${grm}/1000G_phase3_20130502_combined_snpsonly.05.common_UCLA_ASD_pca2_snp_loading 2 \
--out ${grm}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.05.05.common_pca2

# Generate plot and population assignments
Rscript ${ref}/pca_plots_ancestry.R \
${ref}/1000GP_Phase3.sample \
${grm}/1000G_phase3_20130502_combined_snpsonly.05.common_UCLA_ASD_pca2.eigenvec \
${grm}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.05.05.common_pca2.proj.eigenvec \
UCLA_ASD \
${grm}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.05.05.common_pca2 \
${grm}

#--------------------------------------
# Create heterozygosity + missingness plots in R
#--------------------------------------
${software}/plink2 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG \
--missing \
--out ${qcdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG

# LD prune SNPs to generate het file on plink1.9 only
${software}/plink1.9 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG \
--indep 50 5 2 \
--out ${qcdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG

# --het file on plink1.9 only
${software}/plink1.9 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG \
--extract ${qcdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.prune.in \
--het \
--out ${qcdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG

Rscript ${ref}/heterozygosity_inbreeding_check.R \
${qcdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.het \
${qcdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_RAF_1KG.smiss \
${grm}/pc_and_population.txt \
${grm}/ \
"all"

#==============================================================================
#
# Extract EUR
#
#==============================================================================

grep EUR ${grm}/pc_and_population.txt | cut -f1,2 > ${grm}/UCLA_ASD_EUR.id

${software}/plink2 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid \
--make-bed \
--keep ${grm}/UCLA_ASD_EUR.id \
--hwe 0.000001 \
--maf 0.01 \
--geno 0.05 \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR

# Get allele frequency
${software}/plink2 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR \
--freq \
--keep-allele-order \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR

# Subset to brainCTP
${software}/plink2 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR \
--make-bed \
--keep ~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAP_LIBD_ASDbrain_genoid_IID_brainCTP.id \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP

#==============================================================================
#
# Calculate PCs for this dataset
#
#==============================================================================

# LD prune SNPs to generate het file on plink1.9 only
${software}/plink \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR \
--indep 50 5 2 \
--keep ~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAP_LIBD_ASDbrain_genoid_IID_brainCTP.id \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP

# Make GRM based on LD-pruned SNPs
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR \
--keep ~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAP_LIBD_ASDbrain_genoid_IID_brainCTP.id \
--extract ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP.prune.in \
--make-grm \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_ldprune

# Generate genotyping PCs
${software}/gcta_1.93.2beta/gcta64 \
--grm ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_ldprune \
--pca 20 \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_ldprune_pca20

# Identify outliers on PCA and remove
# Removing samples that are 3+ stds away from the mean on PC1, PC2 or PC3 (AN01093, AN00764)
${software}/plink \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR \
--indep 50 5 2 \
--keep ~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAP_LIBD_ASDbrain_genoid_IID_brainCTP.id \
--remove ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_ldprune_pca20.outlier \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_rmoutlier

# Make GRM based on LD-pruned SNPs
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR \
--keep ~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAP_LIBD_ASDbrain_genoid_IID_brainCTP.id \
--remove ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_ldprune_pca20.outlier \
--extract ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_rmoutlier.prune.in \
--make-grm \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_ldprune_rmoutlier

# Generate genotyping PCs
${software}/gcta_1.93.2beta/gcta64 \
--grm ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_ldprune_rmoutlier \
--pca 20 \
--out ${impdir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR_brainCTP_ldprune_rmoutlier_pca20
