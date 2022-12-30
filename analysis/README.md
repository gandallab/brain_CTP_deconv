# README: `analysis/`

## `1_deconvolution_method_comparison/`

These scripts compare various deconvolution methods, informing our final choice of algorithm.

Folder contains one .Rmd script per dataset, accompanied by a .html output file:

- `ctp_validation_Wong2018_ASDmethylation.Rmd`: UCLA_ASD dataset
- `ctp_validation_Jaffe2018_age0.Rmd`: LIBD dataset
- `ctp_validation_ROSMAP.Rmd`: ROSMAP dataset

Further detail is provided in the manuscript (Supplementary Note 1), but as an overview, the comparisons included:

| Prefix in .Rmd| Reference data   | Marker probe identification method   | Algorithm         |
|---------------|------------------|---------------------|-------------------|
| hseq          | sequencing Luo et al. 2022  | "extremes": %methylation 60/40  | Houseman          |
| mcc           | sequencing Luo et al. 2022  | "extremes": %methylation 60/40  | methylCC          |
| celfie        | sequencing Luo et al. 2022  | "extremes": %methylation 60/40  | CelFIE            |
| h             | array (aggregated datasets) | TOAST: eBayes          | Houseman          |
| sSV           | (reference free) | NA                    | smartSVA          |
| vae           | (reference free) | NA                    | VAE               |

We also experimented with a variety of marker probe identification methods.

- For the sequencing-based reference data, we also tried a chi-squared method but found the "extremes" method was superior.
- For the array-based reference data, we also tried the rowFtests method implemented in minfi. The results were extremely similar.

## `2_ctp_differences/`

Analyses for brain CTP differences associated with diagnosis, age and sex.

Folder contains paired .Rmd and .html files for:

- `ctp_differences`: includes exploratory data anlyses for diagnostic, age and sex differences for the ROSMAP, LIBD and UCLA_ASD datasets, experimenting with covariates, transformations and statistical methods.
- `ctp_differences_local_plots`: this is really a continuation of above but this contains more of the finalised results of this analyses. 
- `ctp_differences_replication_BDR`: equivalent analyses in the BDR replication dataset (for Alzheimer's disease)

## `3_pgs/`

Analyses for associations between brain CTPs and polygenic scores (PGS) for various neuropsychiatric traits.
We used SBayesR to generate polygenic scores.

Folder contains paired .Rmd and .html files for:

- `pgs_analysis`: exploratory, finalised, and sensitivity analyses for these analyses.

Scripts for PGS generation are provided in the `qc/` folder.

## `4_gwas/`

For the GWAS association method, we took a linear model approach within each dataset (ie. ROSMAP, LIBD, UCLA_ASD). We had also tried MLMA-LOCO but ran into convergence errors, probably because of our small datasets.

For the brain CTP phenotypes, we took 2 approaches:

1. The n=7 brain CTP phenotypes were clr-transformed (to deal with proportionality), with subsequent inverse normal transformation (to induce normality as clr-transformation can create outliers).
2. Performing compositionally-aware PCA on the raw brain CTPs, and taking these PCs as phenotypes (CTP_PCs). This avoids the need for clr-transformation, but the drawback is that drawing conclusions about brain CTPs requires the reader to map each CTP_PC to its respective PCA loadings. For example, CTP_PC1 below corresponds to low Exc/Inh/Astro proportions and high OPC/Micro proportions.

<img width="227" alt="image" src="https://user-images.githubusercontent.com/19381296/210031331-5935b756-42b5-43b6-8f10-824940d3e67c.png">

We then meta-analysed the 3 datasets using METAL.

We also performed COJO to identify conditionally independent loci.

Folder contains:

- 

