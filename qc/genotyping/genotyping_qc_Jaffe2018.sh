genoid=Illumina_1M
genoid=Illumina_h650

orig_dir=~/shared-gandalm/GenomicDatasets/LIBD/Genotypes/LIBD-QCd
geno_dir=~/shared-gandalm/brain_CTP/Data/genotyping/Jaffe2018
software=~/project-gandalm/software
ref=~/shared-gandalm/brain_CTP/Data/genotyping/ref
impvcf=${geno_dir}/${genoid}/imputed
impdir=${geno_dir}/${genoid}/imputed
grm_dir=~/shared-gandalm/brain_CTP/Data/genotyping/Jaffe2018/grm

#==============================================================================
#
# Genotyping QC for LIBD samples
#
#==============================================================================

module load plink

# Extract data
${software}/plink2 \
--bfile ${orig_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid} \
--make-bed \
--keep ${geno_dir}/Jaffe2018_brnum.id \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018

# Genotyping QC filters
# - use same filters as Cindy to be concordant
# - use plink2 to deal with HWE sex differences automatically(?)
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018 \
--make-bed \
--mind 0.1 \
--geno 0.05 \
--hwe 0.000001 \
--maf 0.01 \
--output-chr 26 \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC

# Sex check
plink \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC \
--check-sex \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC
awk '{print $1, $2, $4}' ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC.sexcheck  | awk 'NR>1' > ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC.to.update.sex.txt

# Create frequency file
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC \
--freq \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC

# DO MANUALLY - REMOVE!
module load bcftools
module load vcftools
strand=${ref}/Human1M-Duov3_B-b37.Source.strand
#--------------------------------------
# 1. Flip minus strand 
# [https://gist.github.com/snewhouse/a02386e4facce2c95c8e0ef2740877c3]
cd ${geno_dir}/${genoid}
cp ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC.bim ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC.bim_ori

# Reassign IDs to match the strand file
Rscript ${ref}/format_bimfile_to_strand.R \
${strand} \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC.bim

# Remove duplicated SNPs so the update_build.sh doesn't crash
# Add merge-x flag here, as GSAv2 has split-x vs GSAv1 is merge-x
# Keep merge-x to make liftover script work
# NOT --sort-vars, or messes up the SNP order!
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC \
--rm-dup force-first \
--make-bed \
--merge-x \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC_forflip

# Run flip script
${ref}/update_build.sh \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC_forflip \
${strand} \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC_flipped

# Convert to rsid
cp ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC_flipped.bim ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC_flipped.bim_ori
Rscript ${ref}/rsid_liftover_aut_x.R  \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC_flipped.bim

# Sort variants, as have preserved the PAR region order in rsid_liftover_aut_x.R
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC_flipped \
--sort-vars \
--make-pgen \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC_flipped

${software}/plink2 \
--pfile ${DATA}/${prefix}_${cohort}_cleaned_aut_x_flipped \
--make-bfile \
--out ${DATA}/${prefix}_${cohort}_cleaned_aut_x_flipped


# Execute Will Rayner's pre-imputation checks
${software}/HRC-1000G-check-bim.pl \
-b ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC.bim \
-f ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC.afreq \
-r ${ref}/1000GP_Phase3_combined.legend \
-g \
-p ALL \
-o ${geno_dir}/${genoid}

sh ${geno_dir}/${genoid}/Run-plink.sh

# Merge together updated genotyping files (autosomes only)
ls ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated-chr*.bim | \
sed 's/.\{4\}$//' | grep -v 23 > ${geno_dir}/${genoid}/pmerge.list
plink \
--merge-list ${geno_dir}/${genoid}/pmerge.list \
--make-bed \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE

# Check RAF + find allele outliers + change ref
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE \
--freq \
--keep-allele-order \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF

Rscript \
${ref}/RAF_outlier_filter_1KG.R \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF.afreq \
${ref}/1000GP_Phase3_combined_rsid_clean.legend \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG

${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE \
--ref-allele ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.afreq2recode 2 1 \
--exclude ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.afreqoutlier \
--set-hh-missing \
--make-bed \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG

for chr in {1..23}; do \
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG \
--chr ${chr} \
--make-bed \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG_${chr}; \
done

# Convert to VCF (autosomes only for the moment)
for chr in {1..22}; do \
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG \
--chr ${chr} \
--recode vcf \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG_chr${chr}
done

module load vcftools
module load htslib
# http://www.htslib.org/download/
# PATH=~/project-gandalm/software/bin/:$PATH
for chr in {1..22}; do \
vcf-sort ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG_chr${chr}.vcf | bgzip -c > ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG_chr${chr}.vcf.gz
done

#==============================================================================
#
# QC METRICS AND CHECKS
#
#==============================================================================

#--------------------------------------
# PCA
#--------------------------------------
# 1. Pre-PCA
# a) Filter out the rare variants
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG \
--maf 0.05 \
--make-bed \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.05

# b) sort the order of the SNPs for picking the shared set between reference and target data
awk '{print $2}' ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.05.bim | sort > \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.05.snps.txt
# sort ${REF}.05.SNPs.txt > ${REF}.05.snps_sort.txt # did this in R
comm -12 ${ref}/1000GP_Phase3_combined_rsid_clean.legend.05.snps_sort.txt ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.05.snps.txt > \
${ref}/common.1KG.${genoid}.SNPs.txt

# c) get the common SNPs for the reference and target data
${software}/plink2 \
--bfile ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05 \
--extract ${ref}/common.1KG.${genoid}.SNPs.txt \
--make-bed \
--out ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_${genoid}

${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.05 \
--maf 0.05 \
--extract ${ref}/common.1KG.${genoid}.SNPs.txt \
--make-bed \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.05.common

#--------------------------------------
# 2. PC projection

# a) generate GRM for 1000G
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_${genoid} \
--make-grm \
--out ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_${genoid}

# b) Take first 2 principal components of reference (1000G)
${software}/gcta_1.93.2beta/gcta64 \
--grm ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_${genoid}  \
--pca 2 \
--out ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_${genoid}_pca2

# c) Generate loadings of each SNP in the reference on the above PCs
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_${genoid} \
--pc-loading ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_${genoid}_pca2 \
--out ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_${genoid}_pca2_snp_loading

# d) Project onto reference sample
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.05.common \
--project-loading ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_${genoid}_pca2_snp_loading 2 \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.05.common_pca2

# Generate plot and population assignments
Rscript ${ref}/pca_plots_ancestry.R \
${ref}/1000GP_Phase3.sample \
${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_${genoid}_pca2.eigenvec \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.05.common_pca2.proj.eigenvec \
LIBD_${genoid} \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.05.common \
${geno_dir}/${genoid}

#--------------------------------------
# Create heterozygosity + missingness plots in R
#--------------------------------------
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG \
--missing \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG

# LD prune SNPs to generate het file on plink1.9 only
${software}/plink1.9 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG \
--indep 50 5 2 \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG

# --het file on plink1.9 only
${software}/plink1.9 \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG \
--extract ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.prune.in \
--het \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG

Rscript ${ref}/heterozygosity_inbreeding_check.R \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.het \
${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG.smiss \
${geno_dir}/${genoid}/pc_and_population.txt \
${geno_dir}/${genoid}/ \
"all"

#--------------------------------------
# IBD
#--------------------------------------
plink \
--bfile ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG \
--genome \
--out ${geno_dir}/${genoid}/LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE_RAF_1KG

#==============================================================================
#
# IMPUTED DATA
#
#==============================================================================

for chr in {1..22}; do \
unzip ${geno_dir}/${genoid}/imputed/chr_${chr}.zip
done

for chr in {1..22}; do \
qsub -v ref=${ref},\
impvcf=${impvcf},\
impdir=${impdir},\
chr=${chr} \
${ref}/vcf2plink_rsidliftover.sh
done

for chr in {1..22}; do \
echo ${chr}
zcat ${impdir}/chr${chr}.info.gz | awk -F"\t" 'NR==1{print;next}$7>0.3' > ${impdir}/imputed_chr${chr}.info #to keep the header
done

ls ${impdir}/imputed_chr*.bim | grep -v chr23 | cut -d "." -f1 > ${impdir}/imputed_files_aut.txt

# Merge imputed, add indels, PAR
qsub -v impdir=${impdir},\
indeldir=${geno_dir}/${genoid},\
indelfilen=LIBD_szControl_Genotype_QCd_${genoid}_Jaffe2018_genoQC-updated_MERGE,\
filen=LIBD_${genoid} \
${ref}/merge_imputed_add_indels_PAR.sh

# b) info files
# Autosomes
ls ${impdir}/imputed_chr*.info | grep -v chr23 > ${impdir}/imputed_infofiles_aut.txt
rm ${impdir}/imputed_info_aut.txt # in case there is a pre-existing version
touch ${impdir}/imputed_info_aut.txt
while read file; do
  cat "$file" >> ${impdir}/imputed_info_aut.txt
done <${impdir}/imputed_infofiles_aut.txt

# Identify SNPs with INFO<0.3
awk '$10<0.3' ${impdir}/imputed_info_aut.txt | awk '{print $2}' > ${impdir}/SNPs_aut_info_less0.3.txt

# Imputation QC filters
qsub -v software=${software},\
bfile=${impdir}/imputed_MERGED_1_22_LIBD_${genoid},\
exclude=${impdir}/SNPs_aut_info_less0.3.txt,\
out=${impdir}/imputed_MERGED_1_22_LIBD_${genoid}_info0.3_maf01_hwe1e6 \
${ref}/plink1_imputation_qc.sh

#==============================================================================
#
# Merge panels
#
#==============================================================================

# Merge together genotype files from different genotyping arrays
# vi ${geno_dir}/merge_arrays/merge_arrays.list
# /shared-gandalm/brain_CTP/Data/genotyping/Jaffe2018/Illumina_1M/imputed/imputed_MERGED_1_22_LIBD_Illumina_1M
# /shared-gandalm/brain_CTP/Data/genotyping/Jaffe2018/Illumina_1M/imputed/imputed_MERGED_1_22_LIBD_Illumina_h650

# Run bmerge in plink1
touch ${geno_dir}/merge_arrays/merge_arrays.list
for i in Illumina_1M Illumina_h650; do \
ls ${geno_dir}/${i}/imputed/imputed_MERGED_1_22_LIBD_${i}_info0.3_maf01_hwe1e6.bim | sed 's/....$//' >> ${geno_dir}/merge_arrays/merge_arrays.list
done

qsub -v software=${software},\
mergelist=${geno_dir}/merge_arrays/merge_arrays.list,\
out=${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6 \
${ref}/plink1_merge.sh

# Update IDs to match methylation data n=462 left
# - decide not to use this - get around the IID problem later if needed for PGS
${software}/plink2 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6 \
--make-bed \
--update-ids ${geno_dir}/merge_arrays/brnum_methylsample.id \
--out ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_methylID

# Update rsID
${software}/plink2 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6 \
--make-bed \
--set-all-var-ids @:# \
--output-chr chr26 \
--rm-dup force-first \
--out ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_tmp

${software}/plink2 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_tmp \
--make-bed \
--update-name /u/project/gandalm/shared/GenomicDatasets/ABCD_r201_r1/impute/imputeABCD_July2020/results/TOPMED_postimputation-master/snp_dir/AllChr_Sorted_Tabdelim.txt \
--out ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid

# Extract EUR AFR
for anc in EUR AFR; do
${software}/plink2 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid \
--make-bed \
--keep ${geno_dir}/merge_arrays/pc_and_population_Illumina_1M_Illumina_h650_${anc}.txt \
--hwe 0.000001 \
--maf 0.01 \
--geno 0.05 \
--out ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}
done

# Frequency file for EUR AFR
for anc in EUR AFR; do
${software}/plink2 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc} \
--freq \
--out ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}
done

# Extract brainCTP from EUR AFR
for anc in EUR AFR; do
${software}/plink2 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc} \
--make-bed \
--keep ~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAP_LIBD_ASDbrain_genoid_IID_brainCTP.id \
--out ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP
done

# Frequency file for brainCTP subset of EUR AFR
for anc in EUR AFR; do
${software}/plink2 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP \
--freq \
--out ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP
done

#==============================================================================
#
# Make full GRM
#
#==============================================================================

# Illumina_1M: 152 AFR / 157 EUR / 20 other
# Illumina_h650: 64 AFR / 63 EUR / 5 other / 1 SAS

# cut -f2 ~/shared-gandalm/refGenomes/hapmap/hapmap3.bim > ~/shared-gandalm/brain_CTP/Data/genotyping/ref/hapmap3.snplist

#==============================================================================
# GRM based on all SNPs
#==============================================================================

for anc in EUR AFR; do
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6 \
--keep ${geno_dir}/merge_arrays/pc_and_population_Illumina_1M_Illumina_h650_${anc}.txt \
--make-grm \
--out ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_${anc}
done

#==============================================================================
# Generate PCs within-ancestry
#==============================================================================

impdir=~/shared-gandalm/brain_CTP/Data/genotyping/Jaffe2018/merge_arrays

# LD prune SNPs
# - otherwise, first few PCs might just be capturing LD rather than pop strat
for anc in EUR AFR; do
${software}/plink1.9 \
--bfile ${impdir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_${anc} \
--indep 50 5 2 \
--out ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_${anc}
done

# Generate genotyping PCs within-ancestry
for anc in EUR AFR; do
${software}/gcta_1.93.2beta/gcta64 \
--grm ${impdir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6 \
--extract ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_${anc}.prune.in \
--keep ${geno_dir}/merge_arrays/pc_and_population_Illumina_1M_Illumina_h650_${anc}.txt \
--pca 20 \
--out ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_${anc}_pca20
done

#---------------------
# brainCTP only
#---------------------
# Generate GRM
for anc in EUR AFR; do
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP \
--make-grm \
--out ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP
done

# LD prune SNPs
for anc in EUR AFR; do
${software}/plink1.9 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP \
--indep 50 5 2 \
--out ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP
done

# Generate genotyping PCs within-ancestry
for anc in EUR AFR; do
${software}/gcta_1.93.2beta/gcta64 \
--grm ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP \
--extract ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP.prune.in \
--pca 20 \
--out ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_pca20
done

# List outliers and plot PCs
for anc in EUR AFR; do
    julia ${ref}/pca_outlier.jl \
        ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_pca20.eigenvec \
        ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_ldprune_pca20
done

#---------------------

# Remove outliers (Sample137/Br1878, Sample153/Br1876, Sample664/Br1684)
# LD prune SNPs
for anc in EUR AFR; do
${software}/plink1.9 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP \
--remove ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_ldprune_pca20.outlier \
--indep 50 5 2 \
--out ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_rmoutliers
done

# Generate genotyping PCs within-ancestry
for anc in EUR AFR; do
${software}/gcta_1.93.2beta/gcta64 \
--grm ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP \
--extract ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_rmoutliers.prune.in \
--remove ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_ldprune_pca20.outlier \
--pca 20 \
--out ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_pca20_rmoutliers
done

# List outliers and plot PCs
for anc in EUR AFR; do
    julia ${ref}/pca_outlier.jl \
        ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_pca20_rmoutliers.eigenvec \
        ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_pca20_rmoutliers
done

# Run outlier removal twice on AFR
for anc in AFR; do
    cat ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_pca20_rmoutliers.outlier >> \
        ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_ldprune_pca20.outlier
    # LD prune SNPs and remove outliers
    ${software}/plink1.9 \
        --bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP \
        --remove ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_ldprune_pca20.outlier \
        --indep 50 5 2 \
        --out ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_rmoutliers
    # Generate genotyping PCs within-ancestry
    ${software}/gcta_1.93.2beta/gcta64 \
        --grm ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP \
        --extract ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_rmoutliers.prune.in \
        --remove ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_ldprune_pca20.outlier \
        --pca 20 \
        --out ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_pca20_rmoutliers
    # List outliers and plot PCs
    julia ${ref}/pca_outlier.jl \
        ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_pca20_rmoutliers.eigenvec \
        ${grm_dir}/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_${anc}_brainCTP_pca20_rmoutliers
done

