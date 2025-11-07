# Title       : Genomic Predictions  
# Author      : Ehtisham Khokhar
# University  : New Mexico State University
# Email       : ehtishamshakeel@gmail.com 
# Date        : 2025
# Purpose     : Estimate the Genomic Estimated Breeding Values (GEBVs)
# Description : This script uses the rrBLUP package to perform genomic predictions 
#               and estimate SNP effects and GEBVs based SNP markers 

## 1. Initializing
req_pkgs <- c("readxl", "rrBLUP", "ggplot2", "xlsx")

for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

library(xlsx)
library(readxl)
library(rrBLUP)
library(ggplot2)

## 2. Configuration
data_dir <- "C:/Users/ehtis/OneDrive - New Mexico State University/SUNNY/Research Projects/Mechanical Harvest Projects/genomic prediction/rrblup"
out_dir  <- file.path(data_dir, "outputs")

hapmap_fname <- "NMSU150_KNNimp_BeagleImp.hmp.txt"
pheno_fname  <- "mydata_means.xlsx"
trait_name <- "RED"

if (!dir.exists(data_dir)) stop("data_dir does not exist: ", data_dir)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cv_fold_file_tpl <- file.path(out_dir, "%s_CV_fold_correlations.csv")
cv_pred_file_tpl <- file.path(out_dir, "%s_CV_predictions.csv")
cv_snps_file_tpl <- file.path(out_dir, "%s_Top%d_SNPs_avgCV.csv")
top50_snps_file  <- file.path(out_dir, "Top50_SNPs_DSFG.csv")

## 3. Import Genotypic and Phenotypic Data
hapmap_file <- file.path(data_dir, hapmap_fname)
pheno_file  <- file.path(data_dir, pheno_fname)

if (!file.exists(hapmap_file)) stop("genotype file not found: ", hapmap_file)
if (!file.exists(pheno_file)) stop("phenotype file not found: ", pheno_file)

hapmap_data <- read.table(hapmap_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
head(hapmap_data)

pheno <- read_excel(pheno_file)
str(pheno)
length(unique(pheno$geno))

## 4. Convert HapMap to Numeric
convert_geno <- function(geno_row) {
  geno_row[geno_row == "NA"] <- NA
  alleles <- unique(unlist(strsplit(na.omit(geno_row), split = "")))
  if(length(alleles) != 2) return(rep(NA, length(geno_row)))
  sapply(geno_row, function(g) {
    if(is.na(g)) return(NA)
    if(g == paste0(alleles[1], alleles[1])) return(-1)
    if(g == paste0(alleles[2], alleles[2])) return(1)
    if(all(g %in% alleles)) return(0)
    return(NA)
  })
}

geno_cols <- grep("^GBS", colnames(hapmap_data))
marker_ids <- hapmap_data$`rs#`
geno_numeric <- t(apply(hapmap_data[, geno_cols], 1, convert_geno))
geno_numeric <- data.frame(rs. = marker_ids, geno_numeric, check.names = FALSE)
colnames(geno_numeric)[-1] <- colnames(hapmap_data)[geno_cols]

geno_only <- geno_numeric[, -1]
geno_transposed <- t(geno_only)
geno_transposed <- data.frame(geno_transposed)
colnames(geno_transposed) <- geno_numeric$rs.
rownames(geno_transposed) <- colnames(geno_only)  # sample IDs

## 5. Impute missing genotypes and align phenotypes
Z_imputed <- A.mat(geno_transposed, return.imputed = TRUE)$imputed
cat("Total marker:", ncol(Z_imputed), "\n")

# Align phenotype vector with genotype rownames
if(!all(pheno$geno %in% rownames(Z_imputed))) stop("Some genotypes in phenotype file are missing in genotype matrix")
y_vec <- pheno[[trait_name]][match(rownames(Z_imputed), pheno$geno)]

## 6. Train model and summarize
trained_model <- mixed.solve(y = y_vec, Z = Z_imputed)
BLUE <- as.numeric(trained_model$beta)
BLUP <- as.vector(as.matrix(Z_imputed) %*% trained_model$u)
predicted_train_result <- BLUE + BLUP

obs_pred_df <- data.frame(
  Genotype  = rownames(Z_imputed),
  Observed  = y_vec,
  Predicted = predicted_train_result,
  BLUE      = rep(BLUE, length(y_vec)),
  BLUP      = BLUP
)

summary(predicted_train_result)
cor(y_vec, predicted_train_result, use = "complete.obs")

# Scatter plot
plot_title <- paste("Observed vs Predicted", trait_name, "- Training Set")
trait_plot <- ggplot(obs_pred_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "darkgreen", size = 3, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal(base_size = 14) +
  labs(title = plot_title, x = paste("Observed", trait_name), y = paste("Predicted", trait_name)) +
  annotate("text", x = min(obs_pred_df$Observed, na.rm = TRUE)+5, y = max(obs_pred_df$Predicted, na.rm = TRUE)-5,
           label = paste0("r = ", round(cor(obs_pred_df$Observed, obs_pred_df$Predicted, use = "complete.obs"),2)),
           color = "blue", size = 5)
trait_plot

plot_file <- file.path(out_dir, paste0(trait_name, "_Obs_Pred_plot.png"))
ggsave(filename = plot_file, plot = trait_plot, width = 10, height = 6, dpi = 600, bg = "white")

## 7. 5-fold CV function (with genotype alignment)
genomic_CV <- function(y_vec, Z_matrix, pheno_names, trait_name,
                       k_folds = 5, seed = 123, save_results = TRUE,
                       top_n_snps = 50) {
  library(rrBLUP)
  set.seed(seed)
  
  # Align phenotype vector with Z_matrix rows
  y_vec <- y_vec[match(rownames(Z_matrix), pheno_names)]
  
  n <- length(y_vec)
  folds <- sample(rep(1:k_folds, length.out = n))
  
  fold_cor <- numeric(k_folds)
  all_predictions <- data.frame(Genotype = character(0), Observed = numeric(0), Predicted = numeric(0),
                                BLUE = numeric(0), BLUP = numeric(0), Fold = integer(0))
  snp_effects_list <- list()
  
  for(k in 1:k_folds) {
    train_idx <- which(folds != k)
    valid_idx <- which(folds == k)
    
    Z_train <- Z_matrix[train_idx, ]
    Z_valid <- Z_matrix[valid_idx, ]
    y_train <- y_vec[train_idx]
    y_valid <- y_vec[valid_idx]
    
    model <- mixed.solve(y = y_train, Z = Z_train)
    BLUE <- as.numeric(model$beta)
    BLUP <- as.vector(Z_valid %*% model$u)
    y_pred <- BLUE + BLUP
    
    fold_cor[k] <- cor(y_valid, y_pred, use = "complete.obs")
    snp_effects_list[[k]] <- model$u
    
    all_predictions <- rbind(
      all_predictions,
      data.frame(
        Genotype = rownames(Z_valid),
        Observed = y_valid,
        Predicted = y_pred,
        BLUE = rep(BLUE, length(y_valid)),
        BLUP = BLUP,
        Fold = k
      )
    )
    
    cat("Fold", k, "r =", round(fold_cor[k],3), "\n")
  }
  
  mean_r <- mean(fold_cor, na.rm = TRUE)
  sd_r <- sd(fold_cor, na.rm = TRUE)
  
  # Average SNP effects
  snp_effects_matrix <- do.call(cbind, snp_effects_list)
  avg_snp_effects <- rowMeans(snp_effects_matrix, na.rm = TRUE)
  snp_effects_df <- data.frame(SNP = colnames(Z_matrix), AvgEffect = avg_snp_effects, AbsEffect = abs(avg_snp_effects))
  snp_effects_sorted <- snp_effects_df[order(-snp_effects_df$AbsEffect), ]
  top_snps <- head(snp_effects_sorted, top_n_snps)
  
  if(save_results) {
    write.csv(data.frame(Fold = 1:k_folds, Correlation = fold_cor),
              sprintf(cv_fold_file_tpl, trait_name), row.names = FALSE)
    write.csv(all_predictions, sprintf(cv_pred_file_tpl, trait_name), row.names = FALSE)
    write.csv(top_snps, sprintf(cv_snps_file_tpl, trait_name, top_n_snps), row.names = FALSE)
    cat("CSV files saved in:", out_dir, "\n")
  }
  
  return(list(fold_cor = fold_cor, mean_r = mean_r, sd_r = sd_r,
              predictions = all_predictions, snp_effects = snp_effects_sorted,
              top_snps = top_snps))
}



## 8. Run CV
result_trait <- genomic_CV(
  y_vec = pheno[[trait_name]],
  Z_matrix = Z_imputed,
  pheno_names = pheno$geno,
  trait_name = trait_name,
  k_folds = 5
)

# Access results
result_trait$predictions  # Genotype, BLUE, BLUP, Predicted, Fold
result_trait$fold_cor
result_trait$mean_r


# Check
all(obs_pred_df$Genotype == rownames(Z_imputed)) # Should return TRUE

cv_check <- result_trait$predictions
# Check each fold
for(f in unique(cv_check$Fold)) {
  fold_genos <- cv_check$Genotype[cv_check$Fold == f]
  if(length(fold_genos) != length(unique(fold_genos))) {
    warning(paste("Duplicate genotypes in fold", f))
  }
  cat("Fold", f, "genotypes check passed:", all(fold_genos %in% rownames(Z_imputed)), "\n")
}