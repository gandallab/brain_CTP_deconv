#==============================================================================
#
# Data cleaning: Luo et al. base-resolution single-cell methylation profiles
#
#==============================================================================

#==============================================================================
# Directories
#==============================================================================

data=~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020
scripts=~/shared-gandalm/brain_CTP/Scripts/methylation/qc
out=~/shared-gandalm/brain_CTP/Data/reference_cell_profile/Luo2020

#==============================================================================
# Download pseudobulk single-cell methylation data from Luo et al. 2022
#==============================================================================

# List of the following files
# allc_Exc_L1-3_CUX2.tsv.gz
# allc_Exc_L2-4_RORB.tsv.gz
# allc_Exc_L4_PLCH1.tsv.gz
# allc_Exc_L4-5_FOXP2.tsv.gz
# allc_Exc_L4-5_TOX.tsv.gz
# allc_Exc_L4-6_LRRK1.tsv.gz
# allc_Exc_L5-6_PDZRN4.tsv.gz
# allc_Exc_L6_TLE4.tsv.gz
# allc_Exc_L6_TSHZ2.tsv.gz
# allc_Inh_CGE_LAMP5.tsv.gz
# allc_Inh_CGE_NDNF.tsv.gz
# allc_Inh_CGE_VIP.tsv.gz
# allc_Inh_CGE-MGE_CHST9.tsv.gz
# allc_Inh_MGE_B3GAT2.tsv.gz
# allc_Inh_MGE_CALB1.tsv.gz
# allc_Inh_MGE_PVALB.tsv.gz
# allc_Inh_MGE_UNC5B.tsv.gz
# allc_NonN_Astro_FGF3R.tsv.gz
# allc_NonN_Endo.tsv.gz
# allc_NonN_Micro.tsv.gz
# allc_NonN_Oligo_MBP.tsv.gz
# allc_NonN_OPC.tsv.gz
# allc_hs_fc_UMB_412.tsv.gz
# allc_Outlier.tsv.gz

#==============================================================================
# extract_CGN.sh: Extract CpGs (CG*) (+ optional coverage threhsold)
#==============================================================================

ls ${data}/allc* | cut -f1 -d '.' > ${data}/filen.tmp

cd ${data}

# extract_CGN.sh is a wrapper for:
# zcat ${filen}.tsv.gz | grep CG[A-Z] | awk -v x="$thresh" '$6 >= x' > ${filen}_CGN_c${thresh}.tsv

while read p; do
qsub -v filen="$p",\
thresh=0 \
${scripts}/extract_CGN.sh
done <${data}/filen.tmp

#==============================================================================
# convert_cgn_gr.sh: wrapper to convert to GenomicRanges format
#==============================================================================

ls ${data}/allc* | cut -f11 -d '/' | grep .gz | cut -f1 -d '.' | sed 's/allc_//g' > ${data}/filen_all

while read p; do
echo "$p"
qsub -v scripts=${scripts},\
cgn_dir=${data}/allc_${p}_CGN_c0.tsv,\
filen=allc_${p}_CGN_c0_aut.rds,\
out_dir=${out} \
${scripts}/convert_cgn_gr.sh
done <${data}/filen_all

#==============================================================================
# Overlap with CpG sites
#==============================================================================

# - read in 450K annotation file in hg19 (reduce to small #)
# - overlap on positional coordinates (check how good this is)
# - calculate methylation beta (M / (M + U + pseudocount))

# SUM NEAREST PROBES TO ILMN450K PROBE
# - download 450K and EPIC manifest files from https://zwdzwd.github.io/InfiniumAnnotation
# - overlap CpGs that are in both of the 2 manifest files
ilmn450k_dir=~/shared-gandalm/brain_CTP/Data/methylation/reference/HM450.hg19.manifest_aut_mask_EPICoverlap.rds

while read p; do
echo "$p"
qsub -v scripts=${scripts},\
cgn_dmr_dir=${data}/allc_${p}_CGN_c0_aut.rds,\
ilmn450k_dir=${ilmn450k_dir},\
bp_thresh="100",\
coverage="10",\
filen=allc_${p}_CGN_aut_ilmn450kepic_aggto100bp_c10.rds,\
out_dir=${out} \
${scripts}/overlap_cgn_dmr_ilmn450kepic.sh
done <${data}/filen_all

# Take the directories for the full files
#ls ${data}/*_CGN_c10_aut.rds | grep -v Outlier | grep -v hs_fc_UMB_412 > ${data}/cgn.dir
ls ${data}/*_CGN_c10_aut_ilmn450k_0bp.rds | grep -v Outlier | grep -v hs_fc_UMB_412 > ${data}/cgn_ilmn450k.dir
ls ${data}/*_CGN_c10_aut_ilmn450k_0bp.rds | grep -v Outlier > ${data}/cgn_ilmn450k_all.dir

# Take the directories of the aggto100 files
ls ${data}/*_CGN_aut_ilmn450kepic_aggto100bp_c10.rds | grep -v Outlier | grep -v hs_fc_UMB_412 > ${data}/cgn_ilmn450kepic_aggto100bp_c10.dir
ls ${data}/allc* | cut -f11 -d '/' | grep .gz | cut -f1 -d '.' | sed 's/allc_//g' | grep -v Outlier | grep -v hs_fc_UMB_412 > ${data}/celltype_n_v2

#==============================================================================
# Munge into matrix of cell type (rows) x probe
#==============================================================================

# See .Rmd document

#==============================================================================
# Get histograms of coverage per celltype
#==============================================================================
