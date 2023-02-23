#==============================================================================
#
# Genotyping QC for ROSMAP samples
#
#==============================================================================

genoid=Illumina_HumanOmniExpress
filen=chop.rosmap.euam.vFinal.382only

genoid=Affymetrix_GeneChip6.0
filen=ROSMAP_arrayGenotype

orig_dir=~/shared-gandalm/GenomicDatasets/ROSMAP/ROSMAP_genotyping/genotyping
geno_dir=~/shared-gandalm/GenomicDatasets/ROSMAP/ROSMAP_genotyping/genotyping
software=~/project-gandalm/software
ref=~/shared-gandalm/brain_CTP/Data/genotyping/ref
gendir=${geno_dir}/${genoid}/genotyped
impvcf=${geno_dir}/${genoid}/imputed
impdir=${geno_dir}/${genoid}/imputed

module load plink

# Check missingness (with --mind 0.1, lose 592 people!)
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/download/${filen} \
--missing \
--update-ids ~/shared-gandalm/GenomicDatasets/ROSMAP/ROSMAP_metadata/update_genoid_FIDIID.txt \
--out ${gendir}/check_missing/ROSMAP_genotype_${genoid}

for chr in {1..22}; do \
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/download/${filen} \
--missing \
--chr ${chr} \
--update-ids ~/shared-gandalm/GenomicDatasets/ROSMAP/ROSMAP_metadata/update_genoid_FIDIID.txt \
--out ${gendir}/check_missing/ROSMAP_genotype_${genoid}_${chr}
done

# - check per-chromoome missingness

# Genotyping QC filters
# - use same filters as Cindy to be concordant
# - use plink2 to deal with HWE sex differences automatically(?)
${software}/plink2 \
--bfile ${geno_dir}/${genoid}/download/${filen} \
--missing \
--make-bed \
--mind 0.1 \
--geno 0.05 \
--hwe 0.000001 \
--maf 0.01 \
--output-chr 26 \
--update-ids ~/shared-gandalm/GenomicDatasets/ROSMAP/ROSMAP_metadata/update_genoid_FIDIID.txt \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC

# Sex check
plink \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC \
--check-sex \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC
awk '{print $1, $2, $4}' ${gendir}/ROSMAP_genotype_${genoid}_genoQC.sexcheck  | awk 'NR>1' > ${gendir}/ROSMAP_genotype_${genoid}_genoQC.to.update.sex.txt

# Create frequency file
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC \
--freq \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC

# Execute Will Rayner's pre-imputation checks (with TOPMed in HRC format)
# - N.B. hg38 format ... so would probably have to liftover beforehand
# ${software}/HRC-1000G-check-bim.pl \
# -b ${geno_dir}/${genoid}/ROSMAP_genotype_${genoid}_genoQC.bim \
# -f ${geno_dir}/${genoid}/ROSMAP_genotype_${genoid}_genoQC.afreq \
# -r ${ref}/PASS.Variantsbravo-dbsnp-all.tab.gz \
# -h \
# -o ${geno_dir}/${genoid}

# Execute Will Rayner's pre-imputation checks (with HRC, as ROSMAP is predominantly EUR)
${software}/HRC-1000G-check-bim.pl \
-b ${gendir}/ROSMAP_genotype_${genoid}_genoQC.bim \
-f ${gendir}/ROSMAP_genotype_${genoid}_genoQC.afreq \
-r ${ref}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab \
-h \
-o ${gendir}

# Execute Will Rayner's pre-imputation checks (with 1KG, ALL (maybe should be EUR only))
# ${software}/HRC-1000G-check-bim.pl \
# -b ${geno_dir}/${genoid}/ROSMAP_genotype_${genoid}_genoQC.bim \
# -f ${geno_dir}/${genoid}/ROSMAP_genotype_${genoid}_genoQC.afreq \
# -r ${ref}/1000GP_Phase3_combined.legend \
# -g \
# -p EUR \
# -o ${geno_dir}/${genoid}

sh ${gendir}/Run-plink.sh

# Merge together updated genotyping files (autosomes only)
#ls ${geno_dir}/${genoid}/ROSMAP_genotype_${genoid}_genoQC-updated-chr*.bim | sed 's/.\{4\}$//' | grep -v 23 > ${geno_dir}/${genoid}/pmerge.list
ls ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated-chr*.bim | sed 's/.\{4\}$//' > ${gendir}/pmerge.list
plink \
--merge-list ${gendir}/pmerge.list \
--make-bed \
--output-chr M \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE

# Remove the people with within-chromosome missingness ONLY for AFFY chip
genoid=Affymetrix_GeneChip6.0
filen=ROSMAP_arrayGenotype
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE \
--make-bed \
--remove ${gendir}/affy_chr8_chr11_chr12_chr19_miss0.3.rmind \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE

# Convert to VCF (autosomes only for the moment)
for chr in {1..22}; do \
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE \
--chr ${chr} \
--recode vcf \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}
done
# chrX, Michigan Imputation throws an error
# - may be due to heterozygous haploid issue for males https://www.biostars.org/p/247101/
# - try --set-hh-missing
# - separate males (run SHAPEIT2 with no phasing) and female (run EAGLE for phasing)
# - separate nonPAR vs PAR regions
chr=X
plink \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE \
--chr ${chr} \
--split-x b37 no-fail \
--make-bed \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR
plink \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR \
--chr ${chr} \
--not-chr XY \
--make-bed \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR
plink \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR \
--set-hh-missing \
--make-bed \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh \
--keep-males \
--make-bed \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh_males
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh \
--keep-females \
--make-bed \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh_females
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh_males \
--chr ${chr} \
--recode vcf \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh_males
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh_females \
--chr ${chr} \
--recode vcf \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh_females

# To get bgzip to work:
module load vcftools
module load htslib
# http://www.htslib.org/download/
# PATH=~/project-gandalm/software/bin/:$PATH
for chr in {1..22}; do \
vcf-sort ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}.vcf | bgzip -c > ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}.vcf.gz
done
chr=X
vcf-sort ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh_males.vcf | bgzip -c > ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh_males.vcf.gz
vcf-sort ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh_females.vcf | bgzip -c > ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_chr${chr}_noPAR_hh_females.vcf.gz

#==============================================================================
#
# QC METRICS AND CHECKS
#
#==============================================================================

#-------------------------------------
# Identify the samples with high chr12 missingness (among other regions) ...
# --> DECISION:
# - remove individuals with high missingness (>0.3) for chr8, 11, 12, 19 as these chunks have lots of SNPs in them
# - don't worry about the other chunks that fail as don't have many SNPs in them anyway (chr9, chr17)
# - unsure what to do about the chr19 chunk (2597 SNPs fail)
#-------------------------------------
for i in 12 11 8 9 17 19 20; do \
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_Affymetrix_GeneChip6.0_genoQC-updated-chr${i} \
--missing \
--out ${gendir}/imputation_error_chr${i}_check
done

# From here, read the .smiss files into R (eg. imputation_error_chr12_check.smiss)
# - manually check missingness
# - eg. ind 3 individuals with chr12 missingness >0.95!!
# MAP37865636
# MAP50109503
# MAP38606123
# Put these IDs in imputation_error_chr12_11_8_9_17_19_20_remove.id

#-------------------------------------
# Align to 1KG for pop comparison
#-------------------------------------

# Check RAF + find allele outliers + change ref
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE \
--freq \
--keep-allele-order \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE

module load R/4.0.2
Rscript \
${ref}/RAF_outlier_filter_1KG.R \
${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE.afreq \
${ref}/1000GP_Phase3_combined_rsid_clean.legend \
${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE

${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE \
--ref-allele ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE.afreq2recode 2 1 \
--exclude ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE.afreqoutlier \
--set-hh-missing \
--make-bed \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG

#--------------------------------------
# PCA
#--------------------------------------
# 1. Pre-PCA
# a) Filter out the rare variants
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG \
--maf 0.05 \
--make-bed \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.05

# b) sort the order of the SNPs for picking the shared set between reference and target data
awk '{print $2}' ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.05.bim | sort > ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.05.snps.txt
# sort ${REF}.05.SNPs.txt > ${REF}.05.snps_sort.txt # did this in R
comm -12 ${ref}/1000GP_Phase3_combined_rsid_clean.legend.05.snps_sort.txt ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.05.snps.txt > ${ref}/common.1KG.ROSMAP_${genoid}.SNPs.txt

# c) get the common SNPs for the reference and target data
${software}/plink2 \
--bfile ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05 \
--extract ${ref}/common.1KG.ROSMAP_${genoid}.SNPs.txt \
--make-bed \
--out ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_ROSMAP_${genoid}

${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.05 \
--maf 0.05 \
--extract ${ref}/common.1KG.ROSMAP_${genoid}.SNPs.txt \
--make-bed \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.05.common

#--------------------------------------
# 2. PC projection

# a) generate GRM for 1000G
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_ROSMAP_${genoid} \
--make-grm \
--out ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_ROSMAP_${genoid}

# b) Take first 2 principal components of reference (1000G)
${software}/gcta_1.93.2beta/gcta64 \
--grm ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_ROSMAP_${genoid}  \
--pca 2 \
--out ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_ROSMAP_${genoid}_pca2

# c) Generate loadings of each SNP in the reference on the above PCs
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_ROSMAP_${genoid} \
--pc-loading ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_ROSMAP_${genoid}_pca2 \
--out ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_ROSMAP_${genoid}_pca2_snp_loading

# d) Project onto reference sample
${software}/gcta_1.93.2beta/gcta64 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.05.common \
--project-loading ${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_ROSMAP_${genoid}_pca2_snp_loading 2 \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.05.common_pca2

# Generate plot
Rscript ${ref}/pca_plots_ancestry.R \
${ref}/1000GP_Phase3.sample \
${ref}/1000G/1000G_phase3_20130502_combined_snpsonly.05.common_ROSMAP_${genoid}_pca2.eigenvec \
${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.05.common_pca2.proj.eigenvec \
ROSMAP_${genoid} \
${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.05.common \
${gendir}

#--------------------------------------
# Create heterozygosity + missingness plots in R
#--------------------------------------
${software}/plink2 \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG \
--missing \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG

# LD prune SNPs to generate het file on plink1.9 only
plink \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG \
--indep 50 5 2 \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG

# --het file on plink1.9 only
plink \
--bfile ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG \
--extract ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.prune.in \
--het \
--out ${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG

Rscript ${ref}/heterozygosity_inbreeding_check.R \
${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.het \
${gendir}/ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE_RAF_1KG.smiss \
${geno_dir}/pc_and_population.txt \
${geno_dir}/ \
"all"

#==============================================================================
#
# IMPUTED DATA
#
#==============================================================================

for chr in {1..22}; do \
unzip -P "eDmAInpRG80svR" ${geno_dir}/${genoid}/imputed/chr_${chr}.zip
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

ls ${impdir}/imputed_chr*.bim | grep -v chr23 | sed 's/....$//' > ${impdir}/imputed_files_aut.txt

# Merge imputed, add indels, PAR
qsub -v impdir=${impdir},\
indeldir=${gendir},\
indelfilen=ROSMAP_genotype_${genoid}_genoQC-updated_HRC_MERGE,\
filen=ROSMAP_${genoid} \
${ref}/merge_imputed_add_indels_PAR.sh

# b) info files
# # All
# ls ${impdir}/imputed_chr*.info > ${impdir}/imputed_infofiles.txt
# rm ${impdir}/imputed_info.txt # in case there is a pre-existing version
# touch ${impdir}/imputed_info.txt
# while read file; do
#   cat "$file" >> ${impdir}/imputed_info.txt
# done <${impdir}/imputed_infofiles.txt

# Autosomes
ls ${impdir}/imputed_chr*.info | grep -v chr23 > ${impdir}/imputed_infofiles_aut.txt
rm ${impdir}/imputed_info_aut.txt # in case there is a pre-existing version
touch ${impdir}/imputed_info_aut.txt
while read file; do
  cat "$file" >> ${impdir}/imputed_info_aut.txt
done <${impdir}/imputed_infofiles_aut.txt

# Identify SNPs with INFO<0.3
awk '$10<0.3' ${impdir}/imputed_info_aut.txt | awk '{print $2}' > ${impdir}/SNPs_aut_info_less0.3.txt
#awk '$10<0.8' ${impdir}/imputed_info_aut.txt | awk '{print $2}' > ${impdir}/SNPs_aut_info_less0.8.txt

# Imputation QC filters (put as job as otherwise run out of memory)
qsub -v software=${software},\
bfile=${impdir}/imputed_MERGED_1_22_ROSMAP_${genoid},\
exclude=${impdir}/SNPs_aut_info_less0.3.txt,\
out=${impdir}/imputed_MERGED_1_22_ROSMAP_${genoid}_info0.3_maf01_hwe1e6 \
${ref}/plink1_imputation_qc.sh

#==============================================================================
#
# Merge panels
#
#==============================================================================

# Run bmerge in plink1 (put as job as otherwise run out of memory)
touch ${geno_dir}/merge_arrays/merge_arrays.list
for i in Illumina_HumanOmniExpress Affymetrix_GeneChip6.0; do \
ls ${geno_dir}/${i}/imputed/imputed_MERGED_1_22_ROSMAP_${i}_info0.3_maf01_hwe1e6.bim | sed 's/....$//' >> ${geno_dir}/merge_arrays/merge_arrays.list
done

qsub -v software=${software},\
mergelist=${geno_dir}/merge_arrays/merge_arrays.list,\
out=${geno_dir}/merge_arrays/imputed_MERGED_1_22_ROSMAP_info0.3_maf01_hwe1e6 \
${ref}/plink1_merge.sh

#==============================================================================
#
# Make full GRM
#
#==============================================================================

${software}/gcta_1.93.2beta/gcta64 \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_ROSMAP_info0.3_maf01_hwe1e6 \
--make-grm \
--out ${geno_dir}/grm/imputed_MERGED_1_22_ROSMAP_info0.3_maf01_hwe1e6

# LD prune SNPs to generate het file on plink1.9 only
plink \
--bfile ${geno_dir}/merge_arrays/imputed_MERGED_1_22_ROSMAP_info0.3_maf01_hwe1e6 \
--indep 50 5 2 \
--out ${geno_dir}/grm/imputed_MERGED_1_22_ROSMAP_info0.3_maf01_hwe1e6

# Generate genotyping PCs
${software}/gcta_1.93.2beta/gcta64 \
--grm ${geno_dir}/grm/imputed_MERGED_1_22_ROSMAP_info0.3_maf01_hwe1e6 \
--extract ${geno_dir}/grm/imputed_MERGED_1_22_ROSMAP_info0.3_maf01_hwe1e6.prune.in \
--pca 20 \
--out ${geno_dir}/grm/imputed_MERGED_1_22_ROSMAP_info0.3_maf01_hwe1e6_pca20
