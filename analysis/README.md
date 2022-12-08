# README: analysis/ folder

## 1_deconvolution_method_comparison/

These scripts compare various deconvolution methods, informing our final choice of algorithm.

The comparisons include:

| Prefix        | Reference        | Differential test   | Algorithm         |
|---------------|------------------|---------------------|-------------------|
| hseq          | Luo et al. 2022  | %methylation 60/40  | Houseman          |
| mcc           | Luo et al. 2022  | ???                 | methylCC          |
| celfie        | Luo et al. 2022  | ???                 | CelFIE            |
| h             | aggregate array  | ??? TOAST?          | Houseman          |
| sSV           | (reference free) |                     | smartSVA          |
| vae           | (reference free) |                     | VAE               |

- Houseman with sequencing  %methylation differential probe ()

There is one .Rmd script per dataset, accompanied by a .html output file:

- `ctp_validation_Wong2018_ASDmethylation.Rmd`: UCLA_ASD dataset
- `ctp_validation_Jaffe2018_age0.Rmd`: LIBD dataset
- `ctp_validation_ROSMAP.Rmd`: ROSMAP dataset

## 2_ctp_differences/



## 3_pgs/

## 4_gwas/