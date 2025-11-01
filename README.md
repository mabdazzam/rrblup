# rrBLUP Genomic Prediction

This repository contains an R script for estimating **Genomic Estimated
Breeding Values (GEBVs)** using the **rrBLUP** package.  The pipeline applies
**ridge regression best linear unbiased prediction (RR-BLUP)** to predict
genomic values based on genome-wide SNP data.

## Author
**Ehtisham Khokhar**  New Mexico State University
ehtishamshakeel@gmail.com

## Description
This R script is designed to perform **genomic prediction** by
integrating genotype and phenotype data.  It allows users to assess the
accuracy of prediction models through **cross-validation** and extract **SNP
effects** for further interpretation.

### Workflow Overview
1. **Data Import and Preparation**
   - Reads SNP marker data from a HapMap file (`.hmp.txt`)
   - Imports phenotypic data for target traits (`.xlsx`)
   - Converts genotypes into numeric format suitable for rrBLUP
   - Handles missing data with imputation

2. **Genomic Prediction using rrBLUP**
   - Uses `rrBLUP::mixed.solve` to estimate marker effects and GEBVs
   - Implements **ridge regression** to account for multicollinearity among
     SNPs
   - Calculates predicted trait values for each individual

3. **Cross-Validation**
   - Employs **5-fold cross-validation** to evaluate predictive ability
   - Splits the dataset into training and testing sets
   - Calculates correlation (r) between predicted and observed phenotypes

## Output
Running this script will generate:

- **Predicted GEBVs** for each genotype
- **Cross-validation statistics:** fold-wise correlations and mean predictive
  ability (r)
- **SNP effect estimates:** average SNP effects across folds
- **Top SNPs** based on absolute effect size
- **Plots:**
  - Observed vs predicted trait scatter plots
  - Cross-validation fold barplots
- **CSV files (optional):**
  - `DSFG_CV_fold_correlations.csv`
  - `DSFG_CV_predictions.csv`
  - `DSFG_Top50_SNPs_avgCV.csv`
  - `Top50_SNPs_DSFG.csv`


## Reference
Endelman, J. B. (2011). **Ridge Regression and Other
Kernels for Genomic Selection with R Package rrBLUP.** Plant Genome*,
4(3), 250–255.
