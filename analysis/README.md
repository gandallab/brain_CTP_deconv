# README: analysis/ folder

## 1_deconvolution_method_comparison/

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

## 2_ctp_differences/



## 3_pgs/

## 4_gwas/
