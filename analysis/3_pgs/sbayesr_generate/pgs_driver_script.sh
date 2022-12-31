#==============================================================================
# 
# PGS
#
#==============================================================================

software=~/project-gandalm/software
sumstats_orig=~/shared-gandalm/GWAS
ref=~/shared-gandalm/brain_CTP/Data/genotyping/ref
sumstats=~/shared-gandalm/brain_CTP/Data/integration/pgs/sumstats
scripts=~/shared-gandalm/brain_CTP/Scripts/integration/pgs
mldm=~/shared-gandalm/brain_CTP/Data/genotyping/ref/ldm/band_ukb_10k_hm3/ukb10k.mldm
pgs=~/shared-gandalm/brain_CTP/Data/integration/pgs/scores
bfile_target=~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAP_ASDbrain_LIBD_imputed_MERGE_1_22_QC

module load R/4.0.2
module load gcc

#==============================================================================
# Munge sumstats
#==============================================================================

# N.B. the Jansen et al. 2019 Alzheimer's disease sumstats ran into convergence issues
azd_file=ALZ.Marioni.2018/AD_UKB_parents_IGAP_07Feb2019.txt
Rscript ${scripts}/prs_format_gwas.R \
${sumstats_orig}/${azd_file} \
"SNP" "A1" "A2" "freq" "b" "se" "p" "N" \
"NA" \
${sumstats} \
"ALZ_Marioni2018.sumstats"

# N.B. the SCZ PGC3 sumstats ran into convergence issues
scz_file=SCZ.Pardinas.PGC.2018/SCZ.CLOZUK.2018.COJO.rsid.removedoublehead.ma
#awk {'print $0" "35802'} ${sumstats_orig}/${scz_file} > ${sumstats}/SCZ.CLOZUK_rsFiltered_oneMHCsnp.txt_N
Rscript ${scripts}/prs_format_gwas.R \
${sumstats_orig}/${scz_file} \
"SNP" "A1" "A2" "freq" "b" "se" "p" "n" \
"NA" \
${sumstats} \
"SCZ_Pardinas2018.sumstats"

height_file=Height.Yengo.2018/Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz
Rscript ${scripts}/prs_format_gwas.R \
${sumstats_orig}/${height_file} \
"SNP" "Tested_Allele" "Other_Allele" "Freq_Tested_Allele_in_HRS" "BETA" "SE" "P" "N" \
"NA" \
${sumstats} \
"Height_Yengo2018.sumstats"

ea_file=EduYears.Lee.SSGAC.2018/GWAS_EA_excl23andMe.txt
#awk {'print $0"\t"766345'} ${sumstats_orig}/${ea_file} > ${sumstats_orig}/${ea_file}_N
Rscript ${scripts}/prs_format_gwas.R \
${sumstats_orig}/${ea_file}_N \
"SNP" "A1" "A2" "EAF" "Beta" "SE" "Pval" "766345" \
"NA" \
${sumstats} \
"EduYears_Lee2018.sumstats"

mdd_file=MDD.Howard.PGC.2019/MDD_PGC_UKB_depression_genome-wide.txt
#awk {'print $0" "500199'} ${sumstats_orig}/${mdd_file} > ${sumstats_orig}/MDD_PGC_UKB_depression_genome-wide.txt_N
Rscript ${scripts}/prs_format_gwas.R \
${sumstats_orig}/${mdd_file}_N \
"MarkerName" "A1" "A2" "Freq" "LogOR" "StdErrLogOR" "P" "500199" \
"NA" \
${sumstats} \
"MDD_Howard2019.sumstats"

# Already have pre-processed
# - Chronotypes (Jones 2019)
# - IQ (Savage 2018) 

# Get ASD from previous work (this has logOR and also has AF)
# - rename as ASD.Grove.2018_logOR.sumstats

# Intersect HM3 SNPs
filen=SCZ_Pardinas2018_v2
for filen in AZD_Marioni2018 SCZ_Pardinas2018 ASD_Grove2018_logOR Height_Yengo2018 Chronotype_Jones2019 IQ_Savage2018 EduYears_Lee2018 MDD_Howard2019 
do
qsub -v scripts=${scripts},\
f1=${sumstats}/${filen}.sumstats,\
c1=1,\
h1="header1",\
f2=${ref}/hapmap3_autosome.snplist,\
c2=1,\
h2="noheader2",\
mc="SNP",\
ot="separate",\
od=${sumstats},\
n1=${filen}"_HM3_int1.sumstats" \
${scripts}/intersect_2files.sh
done

# Intersect GWAS and target SNPs
for filen in AZD_Marioni2018 SCZ_Pardinas2018 ASD_Grove2018_logOR Height_Yengo2018 Chronotype_Jones2019 IQ_Savage2018 EduYears_Lee2018 MDD_Howard2019 
do
snplist=~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAPb37_ASDbrainb37_LIBDb38_sharedrsid_snplist.bim
qsub -v scripts=${scripts},\
f1=${sumstats}/${filen}"_HM3_int1.sumstats",\
c1=1,\
h1="header1",\
f2=${snplist},\
c2=2,\
h2="noheader2",\
mc="SNP",\
ot="separate",\
od=${sumstats},\
n1=${filen}"_HM3_target_int2.sumstats" \
${scripts}/intersect_2files.sh
done

#filen=AZD_Jansen2019
filen=AZD_Marioni2018
snplist=~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAPb37_ASDbrainb37_LIBDb38_sharedrsid_snplist.bim
qsub -v scripts=${scripts},\
f1=${sumstats}/${filen}".sumstats",\
c1=1,\
h1="header1",\
f2=${snplist},\
c2=2,\
h2="noheader2",\
mc="SNP",\
ot="separate",\
od=${sumstats},\
n1=${filen}"_target_int2.sumstats" \
${scripts}/intersect_2files.sh

filen=SCZ_Pardinas2018
snplist=~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAPb37_ASDbrainb37_LIBDb38_sharedrsid_snplist.bim
qsub -v scripts=${scripts},\
f1=${sumstats}/${filen}"_HM3_int1.sumstats",\
c1=1,\
h1="header1",\
f2=${snplist},\
c2=2,\
h2="noheader2",\
mc="SNP",\
ot="separate",\
od=${sumstats},\
n1=${filen}"_HM3_target_int2.sumstats" \
${scripts}/intersect_2files.sh

filen=ASD_Grove2018_logOR
snplist=~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAPb37_ASDbrainb37_LIBDb38_sharedrsid_snplist.bim
qsub -v scripts=${scripts},\
f1=${sumstats}/${filen}"_HM3_int1.sumstats",\
c1=1,\
h1="header1",\
f2=${snplist},\
c2=2,\
h2="noheader2",\
mc="SNP",\
ot="separate",\
od=${sumstats},\
n1=${filen}"_HM3_target_int2.sumstats" \
${scripts}/intersect_2files.sh

# Plot allele frequency
# Need to flip *all* of the below because *.afreq* FILE GOES BY ALT_FREQS
# Keep separated as some GWAS need allele flipping
# Remove outliers
#imp_dir=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD
#freq=${imp_dir}/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR.afreq
for filen in AZD_Marioni2018 SCZ_Pardinas2018 ASD_Grove2018_logOR Height_Yengo2018 Chronotype_Jones2019 IQ_Savage2018 EduYears_Lee2018 MDD_Howard2019 
do
freq=~/shared-gandalm/brain_CTP/Data/genotyping/megaanalysis/ROSMAP_ASDbrain_LIBD_imputed_MERGE_1_22_QC.afreq
Rscript ${scripts}/compare_raf.R \
${freq} \
"ID" "ALT" "REF" "ALT_FREQS" \
${sumstats}/${filen}"_HM3_target_int2.sumstats" \
"SNP" "A1" "A2" "FREQ" \
"NA" \
${sumstats} \
"ROSMAP_ASDbrain_LIBD_"${filen}"_afreq.png" \
${filen}"_HM3_target_AF_int3.sumstats"
done

#==============================================================================
# Generate PGS (SBayesR)
#==============================================================================

# SBayesR
filen=AZD_Marioni2018
target=ROSMAP_ASDbrain_LIBD
#target=ROSMAP
#bfile=~/shared-gandalm/GenomicDatasets/ROSMAP/ROSMAP_WGS/imputed/HRC/ROSMAP_b37_JointAnalysis_wgsQC_rmdup_rsid_HRC_imputedQC_MERGE_1_22_ROSMAPonly_EUR
qsub -v software=${software},\
scripts=${scripts},\
mldm=${mldm},\
dir=${sumstats}/${filen}"_HM3_target_AF_int3.sumstats",\
bfile=${bfile_target},\
h2=0.7,\
pi=0.95\\,0.03\\,0.01\\,0.01,\
gamma=0\\,0.01\\,0.1\\,1,\
out=${pgs},\
filen=${filen},\
target=${target} \
${scripts}/pgs_sbayesr_mldm.sh

filen=SCZ_Pardinas2018
target=ROSMAP_ASDbrain_LIBD
#target=LIBD
#bfile=~/shared-gandalm/brain_CTP/Data/genotyping/Jaffe2018/merge_arrays/imputed_MERGED_1_22_LIBD_info0.3_maf01_hwe1e6_rsid_EUR
qsub -v software=${software},\
scripts=${scripts},\
mldm=${mldm},\
dir=${sumstats}/${filen}"_HM3_target_AF_int3.sumstats",\
bfile=${bfile_target},\
h2=0.7,\
pi=0.95\\,0.03\\,0.01\\,0.01,\
gamma=0\\,0.01\\,0.1\\,1,\
out=${pgs},\
filen=${filen},\
target=${target} \
${scripts}/pgs_sbayesr_mldm.sh

filen=ASD_Grove2018_logOR
target=ROSMAP_ASDbrain_LIBD
#target=ASDbrain
#bfile=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR
qsub -v software=${software},\
scripts=${scripts},\
mldm=${mldm},\
dir=${sumstats}/${filen}"_HM3_target_AF_int3.sumstats",\
bfile=${bfile_target},\
h2=0.5,\
pi=0.95\\,0.03\\,0.01\\,0.01,\
gamma=0\\,0.01\\,0.1\\,1,\
out=${pgs},\
filen=${filen},\
target=${target} \
${scripts}/pgs_sbayesr_mldm.sh

filen=Height_Yengo2018
target=ROSMAP_ASDbrain_LIBD
#target=ASDbrain
#bfile=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR
qsub -v software=${software},\
scripts=${scripts},\
mldm=${mldm},\
dir=${sumstats}/${filen}"_HM3_target_AF_int3.sumstats",\
bfile=${bfile_target},\
h2=0.8,\
pi=0.95\\,0.03\\,0.01\\,0.01,\
gamma=0\\,0.01\\,0.1\\,1,\
out=${pgs},\
filen=${filen},\
target=${target} \
${scripts}/pgs_sbayesr_mldm.sh

filen=IQ_Savage2018
target=ROSMAP_ASDbrain_LIBD
#target=ASDbrain
#bfile=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR
qsub -v software=${software},\
scripts=${scripts},\
mldm=${mldm},\
dir=${sumstats}/${filen}"_HM3_target_AF_int3.sumstats",\
bfile=${bfile_target},\
h2=0.3,\
pi=0.95\\,0.03\\,0.01\\,0.01,\
gamma=0\\,0.01\\,0.1\\,1,\
out=${pgs},\
filen=${filen},\
target=${target} \
${scripts}/pgs_sbayesr_mldm.sh

filen=Chronotype_Jones2019
target=ROSMAP_ASDbrain_LIBD
#target=ASDbrain
#bfile=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR
qsub -v software=${software},\
scripts=${scripts},\
mldm=${mldm},\
dir=${sumstats}/${filen}"_HM3_target_AF_int3.sumstats",\
bfile=${bfile_target},\
h2=0.8,\
pi=0.95\\,0.03\\,0.01\\,0.01,\
gamma=0\\,0.01\\,0.1\\,1,\
out=${pgs},\
filen=${filen},\
target=${target} \
${scripts}/pgs_sbayesr_mldm.sh

filen=EduYears_Lee2018
target=ROSMAP_ASDbrain_LIBD
#target=ASDbrain
#bfile=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR
qsub -v software=${software},\
scripts=${scripts},\
mldm=${mldm},\
dir=${sumstats}/${filen}"_HM3_target_AF_int3.sumstats",\
bfile=${bfile_target},\
h2=0.4,\
pi=0.95\\,0.03\\,0.01\\,0.01,\
gamma=0\\,0.01\\,0.1\\,1,\
out=${pgs},\
filen=${filen},\
target=${target} \
${scripts}/pgs_sbayesr_mldm.sh

filen=MDD_Howard2019
target=ROSMAP_ASDbrain_LIBD
#target=ASDbrain
#bfile=~/shared-gandalm/GenomicDatasets/PsychENCODE/Genotypes/IndividualStudies/UCLA_ASD/UCLA_ASD_info0.3_maf01_hwe1e6_rsid_EUR
qsub -v software=${software},\
scripts=${scripts},\
mldm=${mldm},\
dir=${sumstats}/${filen}"_HM3_target_AF_int3.sumstats",\
bfile=${bfile_target},\
h2=0.4,\
pi=0.95\\,0.03\\,0.01\\,0.01,\
gamma=0\\,0.01\\,0.1\\,1,\
out=${pgs},\
filen=${filen},\
target=${target} \
${scripts}/pgs_sbayesr_mldm.sh
