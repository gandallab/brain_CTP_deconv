# README: methylation QC

# Overview

## ROSMAP

- started with raw methylation .idat files
- data downloaded from Synapse accession [`syn7357283`](https://www.synapse.org/#!Synapse:syn7357283)

## LIBD

- started with raw methylation .idat files
- data downloaded from GEO accession [`GSE74193`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74193)

## UCLA_ASD

- data already QC-ed using an ~identical pipeline (Wong et al. 2018 Human Molecular Genetics).
- data downloaded from Synapse [`syn8263588`](https://www.synapse.org/#!Synapse:syn8263588)

# Scripts description: `meffil_pipeline/`

- `meffil_driver.sh`: driver script
- `methylation_qc_report.R` and `methylation_qc_report.sh`: scripts to generate the meffil QC reports (see `meffil_output/`)
- `methylation_qc_rmoutliers.R`: script to remove outliers identified by `methylation_qc_report.*`
- `methylation_qc_normalise.R` and `methylation_qc_normalise.sh`: scripts to generate the meffil normalisation report
- `methylation_osca_prep.R`: script to munge data for input into OSCA

# QC for ROSMAP and LIBD using meffil [Min et al. 2018 Bioinformatics](https://academic.oup.com/bioinformatics/article/34/23/3983/5042224): `meffil_output/`

meffil outputs .html files with QC results:

- `qc.report_ROSMAP.html`: meffil QC report, informs sample/probe exclusion
- `qc.report_Jaffe2018_age0.html`: meffil QC report, informs sample/probe exclusion
- `normalization_report_ROSMAP.html`: meffil normalisation report demonstrating control of batch effects after quantile normalisation for ROSMAP samples
- `normalization_report_Jaffe2018_fetal_platerandom.html.html`: meffil normalisation report demonstrating control of batch effects after quantile normalisation for LIBD samples under 18 years of age
- `normalization_report_Jaffe2018_age18_platefixed.html.html`: meffil normalisation report demonstrating control of batch effects after quantile normalisation for LIBD samples 18 years or older

## A note on batch effects ...

Note that there were prominent batch effects in the LIBD dataset after quantile normalisation, and that these were attributable to plate and slide effects (for more detail, see `methylation_unsupervised_eda_Jaffe2018.html` -> Section 1.3 "Plot PCs" which also shows that these effects existed in the original paper). 

There were also large batch effects in the ROSMAP dataset (related to thermocycler batch), which we additionally attempted to mitigate using ComBat and SNM. However, these batch correction methods induced negative beta values, which created troubles for CTP deconvolution, which were not easily dealt with ... Therefore, for the purposes of deconvolution, we elected not to use these additional batch correction methods, included batch as a fixed covariate in the analyses, and paid careful attention to its statistical effects in interpreting our results.
