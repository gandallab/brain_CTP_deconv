# Polygenic score (PGS) generation using SBayesR

PGS were generated for ROSMAP, LIBD and UCLA_ASD participants for various neuropsychiatric traits

`pgs_driver_script.sh` is the driver script for the pipeline, which performs the following steps:

- formatting: `intersect_2files.sh`
- Run SBayesR (GCTB v2.03 beta), check outputs and generate PGS using PLINK2: `pgs_sbayesr_mldm.sh`, `pgs_sbayesr_check.R`. SBayesR tutorial can be found [here](https://cnsgenomics.com/software/gctb/#Tutorial)