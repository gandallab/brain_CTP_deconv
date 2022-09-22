# README: brain_CTP deconvolution script

This directory contains scripts used to perform brain_CTP deconvolution as described in Yap et al. 2022.

## Scripts

`estimate_ctp_houseman.R` : script to deconvolve brain CTPs

**Start here** for "out-of-the-box" deconvolution as performed in Yap et al. 2022

**Inputs**

-   bulk methylation beta matrix: see `methylation_data.tar.gz` for examples (data used in the paper)

-   processed marker probes: `Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_beta.rds`

`marker_probe_selection/` : folder contains additional scripts demonstrating how the marker probes in `Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_beta.rds` were selected.

1.  `reference_munge.sh` : driver .sh script for initial pre-processing of raw sequencing count data (large files, best done on HPC)

    -   Calls:

        -   `extract_CGN.sh` -\> wrapper for `extract_CGN.R`

            -   Function: from raw sequencing count data, extract reads mapping to CG\*. These represent CpG methylation sites, of which a subset are assayed with the Illumina 450K and EPIC array chips

        -   `convert_cgn_gr.sh` -\> wrapper for `convert_cgn_gr.R`

            -   Function: convert text files to GenomicRanges format as these files are easier to manipulate and overlap

        -   `overlap_cgn_dmr_ilmn450kepic.sh` -\> wrapper for `overlap_cgn_dmr_ilmn450kepic.R`

            -   Function: overlap sequenced CGN sites with Illumina 450K/EPIC sites

2.  `reference_munge.Rmd` : .Rmd script munging output from `overlap_cgn_dmr_ilmn450kepic.sh` from long format to wide format

3.  `find_dmr_extremes.R` : .R script identifying marker probes

    1.  Aggregate multiple Excitatory cell sub-types into Exc, and Inhibitory cell sub-types into Inh, as this gives better CpG coverage (within some cell sub-types, the sequencing is somewhat sparse)

    2.  Calculate %methylation at each CpG site: M/(M+U)

    3.  Identify marker CpG probes on the basis of being

        1.  \>=60% methylated in only 1 cell-type and \<= 40% methylated in all other cell-types (ie. a marker probe that is up-methylated)

        2.  \<=40% methylated in only 1 cell-type and \>=60% methylated in all other cell-types (ie. a marker probe that is down-methylated)

    4.  Format and output as `Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_beta.rds`

## Data

Input bulk methylation data (ROSMAP, LIBD, UCLA_ASD):

-   `methylation_data.tar.gz`

Single-cell methylome sequencing data used to identify marker probes

-   raw data: `www.[Ask Chongyuan]`

-   sequencing read counts across methylation sites overlapping with Illumina 450K/EPIC array CpG probes (read counts summed within +/-50bp of the array CpG probe), then filtering for probes with coverage \>=10

    -   Methylated counts: `ilmn450kepic_aggto100bp_celltype_c10_M.txt`

    -   Total counts: `ilmn450kepic_aggto100bp_celltype_c10_MU.txt`

-   sequencing read counts at selected marker probes:

    -   Methylated counts: `ilmn450kepic_aggto100bp_c10_extremes_split6040all_celltype7_M.rds`

    -   Unmethylated counts: `ilmn450kepic_aggto100bp_c10_extremes_split6040all_celltype7_U.rds`

    -   Total counts: `ilmn450kepic_aggto100bp_c10_extremes_split6040all_celltype7_MU.rds`

-   processed marker probes for input into `estimate_ctp_houseman.R` : `Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_beta.rds`
