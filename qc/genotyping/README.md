# README: genotyping QC

# Scripts description

- driver scripts
  - `genotyping_qc_ASDbrain.sh`: driver script for UCLA_ASD genotyping QC [@Daniel: please add?]
  - `genotyping_qc_Jaffe2018.sh`: driver script for LIBD genotyping QC [@Daniel: please add?]
  - `genotyping_qc_ROSMAP.sh`: driver script for ROSMAP genotyping QC [@Daniel: please add?]
- scripts called within driver scripts
  - `RAF_outlier_filter_1KG.R`: check alleles against 1KG
  - `pca_plots_ancestry.R`: PCA plots for population stratification/ancestry
  - `heterozygosity_inbreeding_check.R`: generates plots to QC by heterozygosity and inbreeding
  - `vcf2plink_rsidliftover.sh`: convert imputation VCF output to PLINK bfile format. Also calls `rsid_liftover_aut_x.R`
  - `merge_imputed_add_indels_PRA.sh`: 
  - `plink1_imputation_qc.sh`: merge imputation files and imputation QC

Software used:

- PLINKv1.90b6.24
- GCTAv1.93.2beta
- bcftoolsv1.11
- vcftoolsv0.1.16

Please see the manuscript Methods for further detail on the genotyping QC pipeline.

