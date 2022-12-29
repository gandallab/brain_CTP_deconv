# README: brain_CTP deconvolution script

This directory contains scripts used to perform brain_CTP deconvolution as described in Yap et al. 2023

## Data

Available within sub-directory: `data/`

### Target dataset (bulk methylation array data to deconvolve): 

-   `methylation_data.tar.gz` contains compressed ROSMAP, LIBD, UCLA_ASD beta matrices (used in paper)

### Methylation single-cell sequencing reference data (Luo et al. 2022)

The below files are listed in order of how they are generated

1.  raw data: available from NCBI GEO/SRA with accession number GSE140493 (Luo et al. 2022 Cell Genomics) [https://www.sciencedirect.com/science/article/pii/S2666979X22000271]

2.   sequencing read counts across methylation sites overlapping with Illumina 450K/EPIC array CpG probes (read counts summed within +/-50bp of the array CpG probe), then filtering for probes with coverage \>=10

    -   Methylated counts: `ilmn450kepic_aggto100bp_celltype_c10_M.txt`

    -   Total counts: `ilmn450kepic_aggto100bp_celltype_c10_MU.txt`

3.   sequencing read counts at selected marker probes:

    -   Methylated counts: `ilmn450kepic_aggto100bp_c10_extremes_split6040all_celltype7_M.rds`

    -   Unmethylated counts: `ilmn450kepic_aggto100bp_c10_extremes_split6040all_celltype7_U.rds`

    -   Total counts: `ilmn450kepic_aggto100bp_c10_extremes_split6040all_celltype7_MU.rds`

4.   processed marker probes for input into `estimate_ctp_houseman.R` : `Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_beta.rds`

	- 	**Start here for out-of-the-box brain CTP deconvolution**


### Methylation array reference data (Guintivano et al. 2013)

-   processed methylation array cell-type marker probes: `dlpfc_450k_guintivano_rowftest_cell2.rds` (see Identification of cell-type specific methylome profiles -> Methylation array)

## Scripts

### Deconvolution of brain CTPs

**Start here** for "out-of-the-box" deconvolution, as performed in Yap et al. 2023

The relevant script is `estimate_ctp_houseman.R` which includes:

-   requisite functions: 
	
		- 	`projectCelType` (performs Houseman 2012 deconvolution implemented in minfi by Aryee et al. 2014 Bioinformatics) 
		
		-	`getErrorPerSample` (indicates quality of deconvolution, see Seiler-Vellame et al. 2022 biorXiv)

-   load libraries and read in inputs (see **Data**):

		-   Target dataset of bulk methylation beta matrix: see `methylation_data.tar.gz` for examples (data used in the paper)

		-   Reference dataset of methylation sequencing marker probes: `Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_beta.rds`
		
		- 	Reference dataset of NeuN+/- array marker probes as comparison: `dlpfc_450k_guintivano_rowftest_cell2.rds`

-   munge reference and target datasets

-   deconvolution using the 2 references 
		
		- 	methylation sequencing with 7 cell-types 
		
		- 	methylation array with 2 cell-types (for comparison)

-   plot deconvolved CTPs


### Identification of cell-type specific methylome profiles

These scripts are included within the sub-directory: `marker_probe_selection/`

#### Methylation sequencing (Luo et al. 2022, 7 cell-types)

The below scripts demonstrate how the marker probes in `Luo2020_extremes_dmr_ilmn450kepic_aggto100bp_c10_cell7_split6040all_beta.rds` were selected.

This is the reference of choice used in Yap et al. 2023

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

#### Methylation array (Guintivano et al. 2014, 2 cell-types)

This data is available via Bioconductor (Jaffe and Kaminsky et al. 2022) [http://bioconductor.org/packages/release/data/experiment/html/FlowSorted.DLPFC.450k.html]

1. `reference_munge_DLPFC_Guintivano.R`: identify marker probes using rowFtests - the default method in `minfi::estimateCellCounts`