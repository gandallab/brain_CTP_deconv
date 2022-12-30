# README: `analysis/`

## `1_deconvolution_method_comparison/`

These scripts compare various deconvolution methods, informing our final choice of algorithm.

There is one .Rmd script per dataset, accompanied by a .html output file:

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

Contains paired .Rmd and .html files for:

- `ctp_differences`: includes exploratory data anlyses for diagnostic, age and sex differences for the ROSMAP, LIBD and UCLA_ASD datasets, experimenting with covariates, transformations and statistical methods.
- `ctp_differences_local_plots`: this is really a continuation of above but this contains more of the finalised results of this analyses. 
- `ctp_differences_replication_BDR`: equivalent analyses in the BDR replication dataset (for Alzheimer's disease)

## `3_pgs/`

Analyses for associations between brain CTPs and polygenic scores for various neuropsychiatric traits.

Contains paired .Rmd and .html files for:

- `pgs_analysis`: exploratory, finalised, and sensitivity analyses for these analyses.

## `4_gwas/`

