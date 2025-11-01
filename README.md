# rrBLUP Genomic Prediction

This repository contains an R script for estimating **Genomic Estimated Breeding Values (GEBVs)** in chile pepper using **rrBLUP**.

## Author
**Ehtisham Khokhar**  
New Mexico State University  
📧 ehtishamshakeel@gmail.com  

## Description
This R script performs:
- SNP data processing from HapMap files  
- Genomic prediction using `rrBLUP::mixed.solve`  
- 5-fold cross-validation for model performance  
- Estimation of SNP effects and identification of top markers  

## Dependencies
```r
library(xlsx)
library(readxl)
library(rrBLUP)
library(ggplot2)

## Reference
Endelman, J. B. (2011). Ridge Regression and Other Kernels for Genomic Selection with R Package rrBLUP. The Plant Genome, 4(3), 250–255.
