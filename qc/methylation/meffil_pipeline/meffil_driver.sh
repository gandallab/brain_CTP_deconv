#==============================================================================
#
# brain_CTP: methylation QC pipeline
#
#==============================================================================

#==============================================================================
# Jaffe et al.
#==============================================================================

#filen=Jaffe2018_age18
filen=Jaffe2018_age18_plate
#filen=Jaffe2018_fetal
filen=Jaffe2018_fetal_plate
filen=Jaffe2018_age0

scripts=~/shared-gandalm/brain_CTP/Scripts/methylation/qc/meffil_pipeline
qc=~/shared-gandalm/brain_CTP/Data/methylation/Jaffe2018/processed
analysis=~/shared-gandalm/brain_CTP/Data/methylation/Jaffe2018/analysis

geno_dir="no_geno" ##### UPDATE LATER
ref=~/shared-gandalm/brain_CTP/Data/methylation/reference

#==============================================================================
# ROSMAP
#==============================================================================

filen=ROSMAP
scripts=~/shared-gandalm/brain_CTP/Scripts/methylation/qc/meffil_pipeline
qc=~/shared-gandalm/brain_CTP/Data/methylation/ROSMAP/processed
analysis=~/shared-gandalm/brain_CTP/Data/methylation/ROSMAP/analysis

geno_dir="no_geno" ##### UPDATE LATER
ref=~/shared-gandalm/brain_CTP/Data/methylation/reference

download=~/shared-gandalm/GenomicDatasets/ROSMAP/methylation/download
raw=~/shared-gandalm/GenomicDatasets/ROSMAP/methylation/idat

# touch ${raw}/SYNAPSE_METADATA_MANIFEST_aggregated.tsv
# for i in 5772325072 5772325089 5815381002 5815381005 5815381006 5815381007 5815381008 5815381009 5815381012 5815381014 5815381015 5815381017 5815381018 5815381019 5815381023 5815381025 5815381026 5815381027 5815381028 5815381029 5815381030 5821203010 5821203011 5821203016 5822020001 5822020002 5822020003 5822020004 5822020005 5822020006 5822020007 5822020008 5822020009 5822020010 5822038001 5822038002 5822038003 5822038004 5822038005
# do
# cat ${download}/${i}/SYNAPSE_METADATA_MANIFEST.tsv >> ${raw}/SYNAPSE_METADATA_MANIFEST_aggregated.tsv
# done
# 
# for i in 5772325072 5772325089 5815381002 5815381005 5815381006 5815381007 5815381008 5815381009 5815381012 5815381014 5815381015 5815381017 5815381018 5815381019 5815381023 5815381025 5815381026 5815381027 5815381028 5815381029 5815381030 5821203010 5821203011 5821203016 5822020001 5822020002 5822020003 5822020004 5822020005 5822020006 5822020007 5822020008 5822020009 5822020010 5822038001 5822038002 5822038003 5822038004
# do
# cat ${download}/${i}/SYNAPSE_METADATA_MANIFEST.tsv >> ${raw}/SYNAPSE_METADATA_MANIFEST_aggregated.tsv
# done
#==============================================================================
# Hannon/Spiers et al.
#==============================================================================

filen=HannonSpiers2015
scripts=~/shared-gandalm/brain_CTP/Scripts/methylation/qc/meffil_pipeline
qc=~/shared-gandalm/brain_CTP/Data/methylation/HannonSpiers2015/processed_Hannon
analysis=~/shared-gandalm/brain_CTP/Data/methylation/HannonSpiers2015/analysis
ref=~/shared-gandalm/brain_CTP/Data/methylation/reference



# QC report
qsub -v scripts=${scripts},\
qcdir=${qc},\
filen=${filen},\
geno=${geno_dir} \
${scripts}/methylation_qc_report.sh

# methylation_qc_rmoutliers.R script (run interactively)

# QC normalisation
# - age18: 15 PCs? or 10 (the last value written below)?
# - fetal: 5 PCs
# - ROSMAP: 10 PCs
qsub -v scripts=${scripts},\
qcdir=${qc},\
filen=${filen},\
pc=10 \
${scripts}/methylation_qc_normalise.sh

# qsub -v scripts=${scripts},\
# qcdir=${qc},\
# filen=${filen},\
# pc=5 \
# ${scripts}/methylation_qc_normalise.sh

# Normalised beta
Rscript ${scripts}/methylation_osca_prep.R \
${qc}/norm.beta_${filen}.Robj \
${qc}/norm.beta_${filen}.txt

# FOR ROSMAP with big batch effects
# - read into R and run ComBat for batch and plate
filen=ROSMAP
filen=ROSMAP_ComBatbatchplate_neg0
filen=ROSMAP_ComBatplate_neg0
filen=ROSMAP_SNMplate_neg0

filen=Jaffe2018_age18_ComBatplate_neg0
filen=Jaffe2018_age0 ##### TO DO

# Output from meffil is in transposed format
~/osca \
--tefile ${qc}/norm.beta_${filen}.txt \
--methylation-beta \
--no-fid \
--make-bod \
--out ${qc}/${filen}

# Update the annotation file
~/osca \
--befile ${qc}/${filen} \
--update-opi ${ref}/450K.opi

# Apply filters
~/osca \
--befile ${qc}/${filen} \
--sd-min 0.02 \
--extract-probe ${ref}/450K_aut.probe \
--exclude-probe ${ref}/450K_mask.probe \
--make-bod \
--out ${analysis}/${filen}_sd02_aut_mask

# Output correlation matrix #####
~/osca \
--befile ${analysis}/${filen}_sd02_aut_mask \
--save-r2 \
--out ${analysis}/${filen}_sd02_aut_mask_rsq

# Output matrix
~/osca \
--befile ${analysis}/${filen}_sd02_aut_mask \
--make-efile \
--out ${analysis}/${filen}_sd02_aut_mask.txt 

# Output matrix
~/osca \
--befile ${analysis}/${filen}_sd02_aut_mask \
--make-tefile \
--out ${analysis}/${filen}_sd02_aut_mask_t.txt 

#cut -f2- ${analysis}/${filen}_sd02_aut_mask_t.txt > ${analysis}/${filen}_sd02_aut_mask_t.glint.txt 
sed -e 's/\t/ /g' ${analysis}/${filen}_sd02_aut_mask_t.txt > ${analysis}/${filen}_sd02_aut_mask_t_spacedelim.txt 

# Run PCA analysis on filtered probes
~/osca \
--befile ${analysis}/${filen}_sd02_aut_mask \
--pca 20 \
--out ${analysis}/${filen}_sd02_aut_mask

#---------------------
# ASDbrain
qc=~/shared-gandalm/brain_CTP/Data/methylation/ASD_methylation_brain/analysis
filen=KCL_R01MH094714_ASD_Illumina450K_PFC
analysis=~/shared-gandalm/brain_CTP/Data/methylation/ASD_methylation_brain/analysis
ref=~/shared-gandalm/brain_CTP/Data/methylation/reference

# Apply filters
~/osca \
--befile ${qc}/${filen} \
--sd-min 0.02 \
--extract-probe ${ref}/450K_aut.probe \
--exclude-probe ${ref}/450K_mask.probe \
--make-efile \
--out ${analysis}/${filen}_sd02_aut_mask.txt

# ROSMAP
qc=~/shared-gandalm/brain_CTP/Data/methylation/ROSMAP/processed
filen=ROSMAP
analysis=~/shared-gandalm/brain_CTP/Data/methylation/ROSMAP/analysis
ref=~/shared-gandalm/brain_CTP/Data/methylation/reference

# Apply filters
~/osca \
--befile ${qc}/${filen} \
--sd-min 0.02 \
--extract-probe ${ref}/450K_aut.probe \
--exclude-probe ${ref}/450K_mask.probe \
--make-efile \
--out ${analysis}/${filen}_sd02_aut_mask.txt

#---------------------
# Remove sd02 filter
~/osca \
--befile ${qc}/${filen} \
--extract-probe ${ref}/450K_aut.probe \
--exclude-probe ${ref}/450K_mask.probe \
--make-bod \
--out ${analysis}/${filen}_aut_mask

# Output matrix
~/osca \
--befile ${analysis}/${filen}_aut_mask \
--make-efile \
--out ${analysis}/${filen}_aut_mask.txt 

# Output matrix
~/osca \
--befile ${analysis}/${filen}_aut_mask \
--make-tefile \
--out ${analysis}/${filen}_aut_mask_t.txt 

# Run PCA analysis on filtered probes
~/osca \
--befile ${analysis}/${filen}_aut_mask \
--pca 10 \
--out ${analysis}/${filen}_aut_mask

#---------------------
# REFERENCE for PCA
# Convert to .bod file
~/osca \
--tefile ~/shared-gandalm/brain_CTP/Data/reference_cell_profile/funnorm_full_BrBlDZMWL_beta.txt \
--methylation-beta \
--no-fid \
--make-bod \
--out ~/shared-gandalm/brain_CTP/Data/reference_cell_profile/brain_ref_450k

~/osca \
--befile ~/shared-gandalm/brain_CTP/Data/reference_cell_profile/brain_ref_450k \
--extract-probe ${ref}/450K_aut_rmmask_rmsnp.probe \
--make-bod \
--out ~/shared-gandalm/brain_CTP/Data/reference_cell_profile/brain_ref_450k_aut_mask_rmsnp

#--extract-probe ${ref}/450K_aut.probe \
#--exclude-probe ${ref}/450K_mask.probe \

# Run PCA analysis on filtered probes
~/osca \
--befile ~/shared-gandalm/brain_CTP/Data/reference_cell_profile/brain_ref_450k_aut_mask_rmsnp \
--pca 20 \
--out ~/shared-gandalm/brain_CTP/Data/reference_cell_profile/brain_ref_450k_aut_mask_rmsnp

