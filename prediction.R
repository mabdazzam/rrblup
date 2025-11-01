# Title       : Genomic Predictions  
# Author      : Ehtisham Khokhar
# University  : New Mexico State University
# Email       : ehtishamshakeel@gmail.com 
# Date        : 2025
# Purpose     : Estimate the Genomic Estimated Breeding Values (GEBVs)
# Description : This script uses the rrBLUP package to perform genomic predictions 
#               and estimate SNP effects and GEBVs based SNP markers 

## 1. Initializing
#  -  Import necessary libraries for data manipulation and visualization
library(xlsx)
library(readxl)
library(rrBLUP)
library(ggplot2)

## 2. Import Data
#  -  Set the working directory and import the data
#  -  This section groups and arranges the genotypes from these datasets
data_dir <- "C:/Users/ehtis/OneDrive - New Mexico State University/SUNNY/Research Projects/Mechanical Harvest Projects/genomic prediction/rrblup"

# show what's inside
if (dir.exists(data_dir)) {
  print(list.files(data_dir))
} else {
  stop("data_dir does not exist: ", data_dir)
}

getwd()
list.files()

## 3. Import Genotypic and Phenotypic Data
#  -  Prepare the data: Perform quality control on the HapMap file — filter markers by MAF (≥ 0.05), remove highly heterozygous markers (heterozygosity ≤ 0.2), and impute missing genotypes using LD-KNNi in TASSEL (with optional extra imputation using Beagle)
#  -  Finalize the dataset: Make sure the cleaned and imputed dataset is ready for use
#  -  Import the file: Read the processed HapMap file for downstream analysis

hapmap_file <- file.path(data_dir, "NMSU150_KNNimp_BeagleImp.hmp.txt")
pheno_file  <- file.path(data_dir, "mydata_means.xlsx")

# error handling for non existant files
if (!file.exists(hapmap_file)) {
  stop("genotype file not found: ", hapmap_file)
}
if (!file.exists(pheno_file)) {
  stop("phenotype file not found: ", pheno_file)
}

# read the file
hapmap_data <- read.table("NMSU150_KNNimp_BeagleImp.hmp.txt",
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE,
                          check.names = FALSE,   # avoids R altering column names
                          quote = "",            # prevents issues with quotes
                          comment.char = "")     # prevents skipping lines starting with #
# View first few rows
head(hapmap_data)
# import the pheno data 
pheno <- read_excel(pheno_file)
str(pheno)
length(unique(pheno$geno)) # number of Genotypes

## 4. Function to Convert HapMap Genotypes to Numeric Format
#  -  Purpose: Convert genotype data from character format (e.g., "AA", "AG", "GG") into numeric format (-1, 0, 1) for use in rrBLUP genomic prediction.
#  -  How it works:
#1.  Takes a vector of genotypes for a single marker (row from the HapMap matrix)
#2.  Converts "NA" strings to actual NA values
#3.  Identifies the unique alleles present in that marker (ignores NAs)
#4.  Assigns numeric codes according to rrBLUP scheme
#-   homozygous for the first allele = -1
#-   heterozygous (first allele + second allele) = 0
#-   homozygous for the second allele = 1
#5.  If the genotype does not match the expected pattern or marker has more than 2 alleles, it returns NA for that genotype
#-   Works marker by marker (row by row)
#-   Assumes diploid, bi-allelic SNPs
#-   Handles missing data (NA) safely
convert_geno <- function(geno_row) {
  # Replace "NA" strings with actual NA
  geno_row[geno_row == "NA"] <- NA
  
  # Get unique alleles ignoring NA
  alleles <- unique(unlist(strsplit(na.omit(geno_row), split = "")))
  
  # If not bi-allelic, return NAs
  if(length(alleles) != 2) {
    return(rep(NA, length(geno_row)))
  }
  
  # Convert to numeric codes
  sapply(geno_row, function(g) {
    if(is.na(g)) return(NA)
    if(g == paste0(alleles[1], alleles[1])) return(-1)
    if(g == paste0(alleles[2], alleles[2])) return(1)
    if(all(g %in% alleles)) return(0)
    return(NA)  # Unexpected genotype
  })
}
# Keep Marker IDs and GBS columns
geno_cols <- grep("^GBS", colnames(hapmap_data))  # Genotype columns only
marker_ids <- hapmap_data$`rs#`                   # Marker IDs

# Apply conversion row by row
geno_numeric <- t(apply(hapmap_data[, geno_cols], 1, convert_geno))

# Convert to data frame and restore column names
geno_numeric <- data.frame(rs. = marker_ids, geno_numeric, check.names = FALSE)
colnames(geno_numeric)[-1] <- colnames(hapmap_data)[geno_cols]

# View first few rows of numeric genotype data
head(geno_numeric)

# Remove the rs. column temporarily for transposition
geno_only <- geno_numeric[, -1]

# Transpose the matrix: now rows = individuals, columns = markers
geno_transposed <- t(geno_only)

# Convert to data frame and keep marker IDs as column names
geno_transposed <- data.frame(geno_transposed)

colnames(geno_transposed) <- geno_numeric$rs.  # marker IDs as column names
rownames(geno_transposed) <- colnames(geno_only)  # samples as row names

## 5. Train the model and summerize the results
# Use A.mat() to impute missing genotypes
Z_imputed <- A.mat(geno_transposed, return.imputed = TRUE)$imputed
str(Z_imputed)

# total number of markers to train the model 
cat("Total marker:", ncol(Z_imputed), "\n")

# Run mixed.solve with Z (marker effects)
trained_model_DSFG <- mixed.solve(y = pheno$DSFG, Z = Z_imputed)

trained_model_DSFG$beta   # overall mean

markers_effect_DSFG <- trained_model_DSFG$u     # SNP effects
BLUE_DSFG <- trained_model_DSFG$beta            # intercept (fixed effect)

# Predict genomic values
predict_train_DSFG <- as.matrix(Z_imputed) %*% markers_effect_DSFG

# Add BLUEs to get final predictions
BLUE_DSFG <- as.numeric(trained_model_DSFG$beta)   # force it into a scalar
predicted_train_result_DSFG <- as.vector(predict_train_DSFG) + BLUE_DSFG

# Summary of predictions
summary(predicted_train_result_DSFG)

# Correlation between observed and predicted DSFR
cor(as.vector(pheno$DSFG), predicted_train_result_DSFG, use = "complete")

# Create a data frame for plotting
DSFG_Obs_Pred_plot <- data.frame(
  Observed = pheno$DSFG,
  Predicted = predicted_train_result_DSFG
)

# Create scatter plot with 1:1 reference line
DSFG_plot <- ggplot(DSFG_Obs_Pred_plot, aes(x = Observed, y = Predicted)) +
  geom_point(color = "darkgreen", size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "red", size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Observed vs Predicted Destemming Force (DSFG) - Training Set",
    x = "Observed DSFG",
    y = "Predicted DSFG"
  ) +
  annotate(
    "text",
    x = min(DSFG_Obs_Pred_plot$Observed, na.rm = TRUE) + 5,
    y = max(DSFG_Obs_Pred_plot$Predicted, na.rm = TRUE) - 5,
    label = paste0("r = ",
                   round(cor(DSFG_Obs_Pred_plot$Observed,
                             DSFG_Obs_Pred_plot$Predicted, use = "complete.obs"), 2)),
    color = "blue",
    size = 5
  )

# Display the plot
DSFG_plot

# Save the plot
# ggsave("DSFG_Obs_Pred_plot.png", plot = DSFR_plot, width = 10, height = 6, dpi = 600, bg = "white")

## 6. Create 5 folds for cross-validation
set.seed(123)  # for reproducibility
n <- nrow(Z_imputed)
K <- 5  # number of folds
y <- pheno$DSFG


folds <- sample(rep(1:K, length.out = n))  # assign each individual to one fold

# Initialize vector to store prediction correlations
fold_cor <- numeric(K)

# Cross-validation loop
for (k in 1:K) {
  cat("\nRunning Fold", k, "...\n")
  # Training and validation indices
  train_index <- which(folds != k)
  valid_index <- which(folds == k)
  
  # Partition data
  Z_train <- Z_imputed[train_index, ]
  y_train <- y[train_index]
  Z_valid <- Z_imputed[valid_index, ]
  y_valid <- y[valid_index]
  
  # Train model on training set
  model <- mixed.solve(y = y_train, Z = Z_train)
  
  # Extract effects
  marker_effects <- model$u
  intercept <- as.numeric(model$beta)
  
  # Predict for validation individuals
  y_pred <- as.vector(Z_valid %*% marker_effects) + intercept
  
  # Compute correlation between observed and predicted
  fold_cor[k] <- cor(y_valid, y_pred, use = "complete.obs")
  
  cat("  Fold", k, "correlation (r):", round(fold_cor[k], 3), "\n")
}

# Compute overall predictive ability
mean_r <- mean(fold_cor, na.rm = TRUE)
sd_r <- sd(fold_cor, na.rm = TRUE)


cat("Mean predictive ability (r):", round(mean_r, 3), "\n")
cat("SD of predictive ability:", round(sd_r, 3), "\n")

# Visualize fold results
cv_results <- data.frame(Fold = 1:K, Correlation = fold_cor)


ggplot(cv_results, aes(x = factor(Fold), y = Correlation)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = round(Correlation, 2)), vjust = -0.5, size = 4) +
  theme_minimal(base_size = 14) +
  labs(
    title = "5-Fold Cross-Validation Predictive Ability (DSFG)",
    x = "Fold Number",
    y = "Predictive Ability (r)"
  ) +
  ylim(0, 1)

## 7. Function for k-fold cross-validation
#-    This function will be use to run CV for multiple traits

genomic_CV <- function(y_vec, Z_matrix, trait_name,
                       k_folds = 5, seed = 123,
                       save_results = TRUE,
                       top_n_snps = 50) {
  library(rrBLUP)
  
  # Input check
  if (length(y_vec) != nrow(Z_matrix)) {
    stop("Length of phenotype vector does not match number of individuals in genotype matrix.")
  }
  
  # Create folds
  set.seed(seed)
  n <- length(y_vec)
  folds <- sample(rep(1:k_folds, length.out = n))
  
  # Storage
  fold_cor <- numeric(k_folds)
  all_predictions <- data.frame(Observed = numeric(0),
                                Predicted = numeric(0),
                                Fold = integer(0))
  snp_effects_list <- list()  # store SNP effects per fold
  
  
  # Cross-validation loop
  
  for (k in 1:k_folds) {
    cat("\nRunning Fold", k, "of", k_folds, "...\n")
    
    train_idx <- which(folds != k)
    valid_idx <- which(folds == k)
    
    Z_train <- Z_matrix[train_idx, ]
    y_train <- y_vec[train_idx]
    Z_valid <- Z_matrix[valid_idx, ]
    y_valid <- y_vec[valid_idx]
    
    # Train rrBLUP model
    model <- mixed.solve(y = y_train, Z = Z_train)
    marker_effects <- model$u
    intercept <- as.numeric(model$beta)
    
    # Store SNP effects
    snp_effects_list[[k]] <- marker_effects
    
    # Predict on validation set
    y_pred <- as.vector(Z_valid %*% marker_effects) + intercept
    
    # Compute correlation
    fold_cor[k] <- cor(y_valid, y_pred, use = "complete.obs")
    
    # Store predictions
    all_predictions <- rbind(
      all_predictions,
      data.frame(Observed = y_valid, Predicted = y_pred, Fold = k)
    )
    
    cat("  Fold", k, "correlation (r):", round(fold_cor[k], 3), "\n")
  }
  
  # Compute cross-validated statistics
  
  mean_r <- mean(fold_cor, na.rm = TRUE)
  sd_r <- sd(fold_cor, na.rm = TRUE)
  
  cat("\nMean predictive ability (r):", round(mean_r, 3), "\n")
  cat("SD of predictive ability:", round(sd_r, 3), "\n")
  
  
  # Compute average SNP effects across folds
  
  cat("\nAveraging SNP effects across folds...\n")
  snp_effects_matrix <- do.call(cbind, snp_effects_list)
  avg_snp_effects <- rowMeans(snp_effects_matrix, na.rm = TRUE)
  
  snp_effects_df <- data.frame(
    SNP = colnames(Z_matrix),
    AvgEffect = avg_snp_effects,
    AbsEffect = abs(avg_snp_effects)
  )
  
  snp_effects_sorted <- snp_effects_df[order(-snp_effects_df$AbsEffect), ]
  top_snps <- head(snp_effects_sorted, top_n_snps)
  
  
  # Save results (optional)
  
  if (save_results) {
    fold_file <- paste0(trait_name, "_CV_fold_correlations.csv")
    pred_file <- paste0(trait_name, "_CV_predictions.csv")
    snp_file <- paste0(trait_name, "_Top", top_n_snps, "_SNPs_avgCV.csv")
    
    write.csv(data.frame(Fold = 1:k_folds, Correlation = fold_cor),
              fold_file, row.names = FALSE)
    write.csv(all_predictions, pred_file, row.names = FALSE)
    write.csv(top_snps, snp_file, row.names = FALSE)
    
    cat("\nCSV files saved:\n", fold_file, "\n", pred_file, "\n", snp_file, "\n")
  }
  
  
  # Return all outputs
  
  return(list(
    fold_cor = fold_cor,
    mean_r = mean_r,
    sd_r = sd_r,
    predictions = all_predictions,
    snp_effects = snp_effects_sorted,
    top_snps = top_snps
  ))
}

## 8. Call the function 'genomic_CV' for different traits
# Run 5-fold CV for DSFG and save CSV results
result_DSFG <- genomic_CV(
  y_vec = pheno$DSFG,
  Z_matrix = Z_imputed,
  trait_name = "DSFG",
  k_folds = 5
)

# Access results
result_DSFG$mean_r        # mean predictive ability
result_DSFG$fold_cor      # fold-wise correlation
result_DSFG$predictions   # observed vs predicted for all individuals

# Get SNP effects
snp_effects <- model$u

# Combine with SNP names
snp_df <- data.frame(
  SNP = colnames(Z_imputed),
  Effect = snp_effects
)

# Calculate absolute effect
snp_df$AbsEffect <- abs(snp_df$Effect)

# Sort by largest absolute effect
snp_df_sorted <- snp_df[order(-snp_df$AbsEffect), ]

# Select top 50 SNPs
top50_snps <- head(snp_df_sorted, 50)

# View
head(top50_snps)

# Save to CSV (optional)
write.csv(top50_snps, "Top50_SNPs_DSFG.csv", row.names = FALSE)



