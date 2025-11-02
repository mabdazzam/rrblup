# rrBLUP Genomic Prediction

This repository contains an R script for estimating **Genomic Estimated
Breeding Values (GEBVs)** using the **rrBLUP** package. It applies ridge
regression best linear unbiased prediction (RR-BLUP) to predict genomic values
from genome-wide SNP data by integrating genotype (HapMap) and phenotype
(trait) data, and it also lets users evaluate prediction accuracy through
cross-validation and extract SNP effects for interpretation.

## Workflow Overview
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

### Provided input
This repository is shipped with input files placed in the `inputs/` directory:

- `inputs/NMSU150_KNNimp_BeagleImp.hmp.txt`
- `inputs/mydata_means.xlsx`

You do **not** need to rename or download these to run the script. The script
will also create an `outputs/` directory (if it does not exist) and write all
results there.

### Windows (PowerShell)

1. Create workspace
```powershell
cd $HOME
New-Item -ItemType Directory -Force -Path .\src\rrblup | Out-Null
cd .\src\rrblup
```
2. Clone and run
```powershell
git clone https://github.com/ehtishamsk/rrblup.git
cd .\rrblup

# 3. run (reads from .\inputs, writes to .\outputs)
Rscript .\prediction.R
```
### Linux / macOS

```sh
mkdir -p ~/usr/local/src && cd ~/usr/local/src
git clone git@github.com:ehtishamsk/rrblup.git
cd rrblup

# run (reads from ./inputs, writes to ./outputs)
Rscript prediction.R

# check outputs
ls outputs
```

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


## Author
**Ehtisham Khokhar**  New Mexico State University
ehtishamshakeel@gmail.com

## Reference
Endelman, J. B. (2011). **Ridge Regression and Other
Kernels for Genomic Selection with R Package rrBLUP.** *The  Plant Genome*,
4(3), 250–255.
