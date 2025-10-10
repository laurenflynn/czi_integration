# step5 QC of all round
# We organze the QC and save the post QC data in this step.





library(dplyr)
#install.packages('Seurat')
library(Seurat)
library(patchwork)
library(umap)
library(hdf5r)
library(ggplot2)
library(Nebulosa)
library(sctransform)
library(monocle3)
library(magrittr)
library(R.utils)
#remotes::install_github('satijalab/seurat-wrappers')
library(Matrix)
library(SeuratWrappers)
library(SingleCellExperiment)
library(reshape2)





# 1. Import Data
# 1-1) Round 1
FL_207_091423 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_207_091423_s.rds")
FL_31_091423 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_31_091423_s.rds")
FL_296_091423 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_296_091423_s.rds")
FL_8_091423R <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_8_091423R_s.rds")
FL_171_091423R <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_171_091423R_s.rds")
# FL_1_091823 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_1_091823_s.rds") # removed due to abnormal percent.mt distribution
FL_12_091823 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_12_091823_s.rds")
FL_6_091823 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_6_091823_s.rds")
FL_625_091823 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_625_091823_s.rds")
FL_199_091823 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_199_091823_s.rds")
FL_9_091923 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_9_091923_s.rds")
FL_17_092023 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_17_092023_s.rds")
# FL_26_092023 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_26_092023_s.rds") # removed due to abnormal percent.mt distribution
# FL_9_092023 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_9_092023_s.rds") # removed due to abnormal percent.mt distribution
# FL_390_092023 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_390_092023_s.rds") # removed due to abnormal percent.mt distribution
# FL_171_092023 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_171_092023_s.rds") # removed due to abnormal percent.mt distribution
FL_45_092023 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_45_092023_s.rds")
FL_2_092123 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_2_092123_s.rds")
FL_5_092123 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_5_092123_s.rds")
FL_18_092123 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_18_092123_s.rds")
FL_296_092123 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_296_092123_s.rds")
FL_95_092123 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_95_092123_s.rds")
FL_62_092123 <- readRDS("output/round_2/output_4_step3_Round1_createSeuratObjects/FL_62_092123_s.rds")



# 1-2) Round 2
FL_22NA_72624 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_22NA_72624.rds")
FL_24_8724 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_24_8724.rds")
FL_24_72924 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_24_72924.rds")
FL_25NA_72624 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_25NA_72624.rds")
FL_28_72524 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_28_72524.rds")
FL_28_72624 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_28_72624.rds")
FL_33_72924 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_33_72924.rds")
FL_53_72924 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_53_72924.rds")
FL_59_72624 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_59_72624.rds")
FL_62_72524 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_62_72524.rds")
FL_62_72924 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_62_72924.rds")
FL_167_72624 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_167_72624.rds")
FL_283_72924 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_283_72924.rds")
FL_313_72924 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_313_72924.rds")
FL_314_8724 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_314_8724.rds")
FL_315_8724 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_315_8724.rds")
FL_316_72524 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_316_72524.rds")
FL_321_72624 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_321_72624.rds")
FL_323_8724 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_323_8724.rds")
FL_154_72524 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_154_72524.rds")
FL_316_8724 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_316_8724.rds")
FL_319_8724 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_319_8724.rds")
FL_320_72924 <- readRDS("output/round_2/output_4_step1_createSeuratObjects/FL_320_72924.rds")





# 2. Processing Data
# Adding in percentage of transcripts that map to mitochondrial genes
# 2-1) Round 1
FL_207_091423[["percent.mt"]] <- PercentageFeatureSet(FL_207_091423, pattern = "^MT-")
FL_31_091423[["percent.mt"]] <- PercentageFeatureSet(FL_31_091423, pattern = "^MT-")
FL_296_091423[["percent.mt"]] <- PercentageFeatureSet(FL_296_091423, pattern = "^MT-")
FL_8_091423R[["percent.mt"]] <- PercentageFeatureSet(FL_8_091423R, pattern = "^MT-")
FL_171_091423R[["percent.mt"]] <- PercentageFeatureSet(FL_171_091423R, pattern = "^MT-")
FL_12_091823[["percent.mt"]] <- PercentageFeatureSet(FL_12_091823, pattern = "^MT-")
FL_6_091823[["percent.mt"]] <- PercentageFeatureSet(FL_6_091823, pattern = "^MT-")
FL_625_091823[["percent.mt"]] <- PercentageFeatureSet(FL_625_091823, pattern = "^MT-")
FL_199_091823[["percent.mt"]] <- PercentageFeatureSet(FL_199_091823, pattern = "^MT-")
FL_9_091923[["percent.mt"]] <- PercentageFeatureSet(FL_9_091923, pattern = "^MT-")
FL_17_092023[["percent.mt"]] <- PercentageFeatureSet(FL_17_092023, pattern = "^MT-")
FL_45_092023[["percent.mt"]] <- PercentageFeatureSet(FL_45_092023, pattern = "^MT-")
FL_2_092123[["percent.mt"]] <- PercentageFeatureSet(FL_2_092123, pattern = "^MT-")
FL_5_092123[["percent.mt"]] <- PercentageFeatureSet(FL_5_092123, pattern = "^MT-")
FL_18_092123[["percent.mt"]] <- PercentageFeatureSet(FL_18_092123, pattern = "^MT-")
FL_296_092123[["percent.mt"]] <- PercentageFeatureSet(FL_296_092123, pattern = "^MT-")
FL_95_092123[["percent.mt"]] <- PercentageFeatureSet(FL_95_092123, pattern = "^MT-")
FL_62_092123[["percent.mt"]] <- PercentageFeatureSet(FL_62_092123, pattern = "^MT-")



# 2-2) Round 2
FL_22NA_72624[["percent.mt"]] <- PercentageFeatureSet(FL_22NA_72624, pattern = "^MT-")
FL_24_8724[["percent.mt"]] <- PercentageFeatureSet(FL_24_8724, pattern = "^MT-")
FL_24_72924[["percent.mt"]] <- PercentageFeatureSet(FL_24_72924, pattern = "^MT-")
FL_25NA_72624[["percent.mt"]] <- PercentageFeatureSet(FL_25NA_72624, pattern = "^MT-")
FL_28_72524[["percent.mt"]] <- PercentageFeatureSet(FL_28_72524, pattern = "^MT-")
FL_28_72624[["percent.mt"]] <- PercentageFeatureSet(FL_28_72624, pattern = "^MT-")
FL_33_72924[["percent.mt"]] <- PercentageFeatureSet(FL_33_72924, pattern = "^MT-")
FL_53_72924[["percent.mt"]] <- PercentageFeatureSet(FL_53_72924, pattern = "^MT-")
FL_59_72624[["percent.mt"]] <- PercentageFeatureSet(FL_59_72624, pattern = "^MT-")
FL_62_72524[["percent.mt"]] <- PercentageFeatureSet(FL_62_72524, pattern = "^MT-")
FL_62_72924[["percent.mt"]] <- PercentageFeatureSet(FL_62_72924, pattern = "^MT-")
FL_167_72624[["percent.mt"]] <- PercentageFeatureSet(FL_167_72624, pattern = "^MT-")
FL_283_72924[["percent.mt"]] <- PercentageFeatureSet(FL_283_72924, pattern = "^MT-")
FL_313_72924[["percent.mt"]] <- PercentageFeatureSet(FL_313_72924, pattern = "^MT-")
FL_314_8724[["percent.mt"]] <- PercentageFeatureSet(FL_314_8724, pattern = "^MT-")
FL_315_8724[["percent.mt"]] <- PercentageFeatureSet(FL_315_8724, pattern = "^MT-")
FL_316_72524[["percent.mt"]] <- PercentageFeatureSet(FL_316_72524, pattern = "^MT-")
FL_321_72624[["percent.mt"]] <- PercentageFeatureSet(FL_321_72624, pattern = "^MT-")
FL_323_8724[["percent.mt"]] <- PercentageFeatureSet(FL_323_8724, pattern = "^MT-")
FL_154_72524[["percent.mt"]] <- PercentageFeatureSet(FL_154_72524, pattern = "^MT-")
FL_316_8724[["percent.mt"]] <- PercentageFeatureSet(FL_316_8724, pattern = "^MT-")
FL_319_8724[["percent.mt"]] <- PercentageFeatureSet(FL_319_8724, pattern = "^MT-")
FL_320_72924[["percent.mt"]] <- PercentageFeatureSet(FL_320_72924, pattern = "^MT-")





# 3. QC percent.mt < 1%
x <- 1

# 3-1) Round1
FL_207_091423 <- subset(FL_207_091423, subset = percent.mt < x)
FL_31_091423 <- subset(FL_31_091423, subset = percent.mt < x)
FL_296_091423 <- subset(FL_296_091423, subset = percent.mt < x)
FL_8_091423R <- subset(FL_8_091423R, subset = percent.mt < x)
FL_171_091423R <- subset(FL_171_091423R, subset = percent.mt < x)
FL_12_091823 <- subset(FL_12_091823, subset = percent.mt < x)
FL_6_091823 <- subset(FL_6_091823, subset = percent.mt < x)
FL_625_091823 <- subset(FL_625_091823, subset = percent.mt < x)
FL_199_091823 <- subset(FL_199_091823, subset = percent.mt < x)
FL_9_091923 <- subset(FL_9_091923, subset = percent.mt < x)
FL_17_092023 <- subset(FL_17_092023, subset = percent.mt < x)
FL_45_092023 <- subset(FL_45_092023, subset = percent.mt < x)
FL_2_092123 <- subset(FL_2_092123, subset = percent.mt < x)
FL_5_092123 <- subset(FL_5_092123, subset = percent.mt < x)
FL_18_092123 <- subset(FL_18_092123, subset = percent.mt < x)
FL_296_092123 <- subset(FL_296_092123, subset = percent.mt < x)
FL_95_092123 <- subset(FL_95_092123, subset = percent.mt < x)
FL_62_092123 <- subset(FL_62_092123, subset = percent.mt < x)



# 3-2) Round2
FL_22NA_72624 <- subset(FL_22NA_72624, subset = percent.mt < x)
FL_24_8724 <- subset(FL_24_8724, subset = percent.mt < x)
FL_24_72924 <- subset(FL_24_72924, subset = percent.mt < x)
FL_25NA_72624 <- subset(FL_25NA_72624, subset = percent.mt < x)
FL_28_72524 <- subset(FL_28_72524, subset = percent.mt < x)
FL_28_72624 <- subset(FL_28_72624, subset = percent.mt < x)
FL_33_72924 <- subset(FL_33_72924, subset = percent.mt < x)
FL_53_72924 <- subset(FL_53_72924, subset = percent.mt < x)
FL_59_72624 <- subset(FL_59_72624, subset = percent.mt < x)
FL_62_72524 <- subset(FL_62_72524, subset = percent.mt < x)
FL_62_72924 <- subset(FL_62_72924, subset = percent.mt < x)
FL_167_72624 <- subset(FL_167_72624, subset = percent.mt < x)
FL_283_72924 <- subset(FL_283_72924, subset = percent.mt < x)
FL_313_72924 <- subset(FL_313_72924, subset = percent.mt < x)
FL_314_8724 <- subset(FL_314_8724, subset = percent.mt < x)
FL_315_8724 <- subset(FL_315_8724, subset = percent.mt < x)
FL_316_72524 <- subset(FL_316_72524, subset = percent.mt < x)
FL_321_72624 <- subset(FL_321_72624, subset = percent.mt < x)
FL_323_8724 <- subset(FL_323_8724, subset = percent.mt < x)
FL_154_72524 <- subset(FL_154_72524, subset = percent.mt < x)
FL_316_8724 <- subset(FL_316_8724, subset = percent.mt < x)
FL_319_8724 <- subset(FL_319_8724, subset = percent.mt < x)
FL_320_72924 <- subset(FL_320_72924, subset = percent.mt < x)





### Move to the 5 th section ##################################################
# This section is to test and plot
# 4. To get all the nCount_RNA (umi_counts) 
# 4-1) Extract UMI counts for each sample and sotre in DF
seurat_objects <- list(
  FL_207_091423 = FL_207_091423,
  FL_31_091423 = FL_31_091423,
  FL_296_091423 = FL_296_091423,
  FL_8_091423R = FL_8_091423R,
  FL_171_091423R = FL_171_091423R,
  FL_12_091823 = FL_12_091823,
  FL_6_091823 = FL_6_091823,
  FL_625_091823 = FL_625_091823,
  FL_199_091823 = FL_199_091823,
  FL_9_091923 = FL_9_091923,
  FL_17_092023 = FL_17_092023,
  FL_45_092023 = FL_45_092023,
  FL_2_092123 = FL_2_092123,
  FL_5_092123 = FL_5_092123,
  FL_18_092123 = FL_18_092123,
  FL_296_092123 = FL_296_092123,
  FL_95_092123 = FL_95_092123,
  FL_62_092123 = FL_62_092123,
  FL_22NA_72624 = FL_22NA_72624,
  FL_24_8724 = FL_24_8724,
  FL_24_72924 = FL_24_72924,
  FL_25NA_72624 = FL_25NA_72624,
  FL_28_72524 = FL_28_72524,
  FL_28_72624 = FL_28_72624,
  FL_33_72924 = FL_33_72924,
  FL_53_72924 = FL_53_72924,
  FL_59_72624 = FL_59_72624,
  FL_62_72524 = FL_62_72524,
  FL_62_72924 = FL_62_72924,
  FL_167_72624 = FL_167_72624,
  FL_283_72924 = FL_283_72924,
  FL_313_72924 = FL_313_72924,
  FL_314_8724 = FL_314_8724,
  FL_315_8724 = FL_315_8724,
  FL_316_72524 = FL_316_72524,
  FL_321_72624 = FL_321_72624,
  FL_323_8724 = FL_323_8724,
  FL_154_72524 = FL_154_72524,
  FL_316_8724 = FL_316_8724,
  FL_319_8724 = FL_319_8724,
  FL_320_72924 = FL_320_72924
)

# Create a list to store the UMI count DataFrames
umi_count_dfs <- list()

# Loop through each Seurat object, extract UMI counts, and store in a DataFrame
for (sample_name in names(seurat_objects)) {
  umi_counts <- seurat_objects[[sample_name]]@meta.data$nCount_RNA
  
  # Create a DataFrame with cell names as row names (index) and UMI counts as a column
  umi_df <- data.frame(
    Cell = rownames(seurat_objects[[sample_name]]@meta.data),
    UMI_Counts = umi_counts
  )
  
  # Store in the list
  umi_count_dfs[[sample_name]] <- umi_df
}

# Example: View UMI count DataFrame for one sample
head(umi_count_dfs$FL_207_091423)



# 4-2) Combine all UMI counts for Kruskal-Wallis Test
# Combine all UMI counts into a single DataFrame for Kruskal-Wallis test
all_umi_counts <- bind_rows(lapply(names(umi_count_dfs), function(sample_name) {
  umi_count_dfs[[sample_name]] %>%
    mutate(Sample = sample_name)  # Add sample name as a column
}))

# View the combined dataset
head(all_umi_counts)

# Save
write.csv(all_umi_counts, "output/round_2/output_4_step5_qc_of_round1AND2/UMI_counts_all_samples.csv")



# 4-3) Perform Kruskan-Wallis Test
kruskal_test <- kruskal.test(UMI_Counts ~ Sample, data = all_umi_counts)
# Print the test result
print(kruskal_test)
pairwise.wilcox.test(all_umi_counts$UMI_Counts, all_umi_counts$Sample, p.adjust.method = "bonferroni")

# Convert the matrix into a long format for ggplot
p_matrix <- pairwise.wilcox.test(all_umi_counts$UMI_Counts, all_umi_counts$Sample, p.adjust.method = "bonferroni")$p.value
p_matrix[is.na(p_matrix)] <- 1  # Replace NA values with 1 for visualization
p_long <- melt(as.matrix(p_matrix))

# Plot as heatmap
ggplot(p_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", na.value = "white") +
  labs(title = "Wilcoxon Test Pairwise P-values",
       x = "Sample", y = "Sample", fill = "p-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Benjamini-Yekutieli
pairwise.wilcox.test(all_umi_counts$UMI_Counts, all_umi_counts$Sample, p.adjust.method = "BY")

# Convert the matrix into a long format for ggplot
p_matrix <- pairwise.wilcox.test(all_umi_counts$UMI_Counts, all_umi_counts$Sample, p.adjust.method = "BY")$p.value
p_matrix[is.na(p_matrix)] <- 1  # Replace NA values with 1 for visualization
p_long <- melt(as.matrix(p_matrix))

# Plot as heatmap
ggplot(p_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", na.value = "white") +
  labs(title = "Wilcoxon Test Pairwise P-values",
       x = "Sample", y = "Sample", fill = "p-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Add log-transformed UMI counts to the dataframe
all_umi_counts$log_UMI_Counts <- log1p(all_umi_counts$UMI_Counts)

# View the updated dataframe
head(all_umi_counts)


# Run pairwise Wilcoxon test with Holm correction on log-transformed data
p_matrix_log_holm <- pairwise.wilcox.test(all_umi_counts$log_UMI_Counts, all_umi_counts$Sample, p.adjust.method = "holm")$p.value

# Replace NA values with 1 for visualization
p_matrix_log_holm[is.na(p_matrix_log_holm)] <- 1

# View the updated p-value matrix
print(p_matrix_log_holm)

# Convert the matrix into a long format for ggplot
p_matrix <- pairwise.wilcox.test(all_umi_counts$log_UMI_Counts, all_umi_counts$Sample, p.adjust.method = "BY")$p.value
p_matrix[is.na(p_matrix)] <- 1  # Replace NA values with 1 for visualization
p_long <- melt(as.matrix(p_matrix))

# Plot as heatmap
ggplot(p_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", na.value = "white") +
  labs(title = "Wilcoxon Test Pairwise P-values",
       x = "Sample", y = "Sample", fill = "p-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Modify values for visualization
p_long$significance <- ifelse(p_long$value > 0.05, "Not Significant", "Significant")
p_long
# Create heatmap plot
ggplot(p_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", limits = c(0, 0.05), na.value = "white") +
  labs(title = "Pairwise Wilcoxon Test - p-value Heatmap",
       x = "Sample",
       y = "Sample",
       fill = "p-value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# IQR x 2.5 qc

# Function to filter Seurat object based on IQR × 2.5 UMI threshold
iqr_qc_filter_umi <- function(seurat_object) {
  
  # Extract UMI counts
  umi_counts <- seurat_object$nCount_RNA
  
  # Compute IQR
  Q1 <- quantile(umi_counts, 0.25)
  Q3 <- quantile(umi_counts, 0.75)
  IQR <- Q3 - Q1
  high_umi_threshold <- Q3 + (2.5 * IQR)
  
  # Subset Seurat object to keep only valid cells
  filtered_seurat <- subset(seurat_object, subset = nCount_RNA <= high_umi_threshold)
  
  # Return the filtered Seurat object
  return(filtered_seurat)
}


FL_207_091423 <- iqr_qc_filter_umi(FL_207_091423)
FL_31_091423 <- iqr_qc_filter_umi(FL_31_091423)
FL_296_091423 <- iqr_qc_filter_umi(FL_296_091423)
FL_8_091423R <- iqr_qc_filter_umi(FL_8_091423R)
FL_171_091423R <- iqr_qc_filter_umi(FL_171_091423R)
FL_12_091823 <- iqr_qc_filter_umi(FL_12_091823)
FL_6_091823 <- iqr_qc_filter_umi(FL_6_091823)
FL_625_091823 <- iqr_qc_filter_umi(FL_625_091823)
FL_199_091823 <- iqr_qc_filter_umi(FL_199_091823)
FL_9_091923 <- iqr_qc_filter_umi(FL_9_091923)
FL_17_092023 <- iqr_qc_filter_umi(FL_17_092023)
FL_45_092023 <- iqr_qc_filter_umi(FL_45_092023)
FL_2_092123 <- iqr_qc_filter_umi(FL_2_092123)
FL_5_092123 <- iqr_qc_filter_umi(FL_5_092123)
FL_18_092123 <- iqr_qc_filter_umi(FL_18_092123)
FL_296_092123 <- iqr_qc_filter_umi(FL_296_092123)
FL_95_092123 <- iqr_qc_filter_umi(FL_95_092123)
FL_62_092123 <- iqr_qc_filter_umi(FL_62_092123)
FL_22NA_72624 <- iqr_qc_filter_umi(FL_22NA_72624)
FL_24_8724 <- iqr_qc_filter_umi(FL_24_8724)
FL_24_72924 <- iqr_qc_filter_umi(FL_24_72924)
FL_25NA_72624 <- iqr_qc_filter_umi(FL_25NA_72624)
FL_28_72524 <- iqr_qc_filter_umi(FL_28_72524)
FL_28_72624 <- iqr_qc_filter_umi(FL_28_72624)
FL_33_72924 <- iqr_qc_filter_umi(FL_33_72924)
FL_53_72924 <- iqr_qc_filter_umi(FL_53_72924)
FL_59_72624 <- iqr_qc_filter_umi(FL_59_72624)
FL_62_72524 <- iqr_qc_filter_umi(FL_62_72524)
FL_167_72624 <- iqr_qc_filter_umi(FL_167_72624)
FL_283_72924 <- iqr_qc_filter_umi(FL_283_72924)
FL_313_72924 <- iqr_qc_filter_umi(FL_313_72924)
FL_314_8724 <- iqr_qc_filter_umi(FL_314_8724)
FL_315_8724 <- iqr_qc_filter_umi(FL_315_8724)
FL_316_72524 <- iqr_qc_filter_umi(FL_316_72524)
FL_321_72624 <- iqr_qc_filter_umi(FL_321_72624)
FL_323_8724 <- iqr_qc_filter_umi(FL_323_8724)
FL_154_72524 <- iqr_qc_filter_umi(FL_154_72524)
FL_316_8724 <- iqr_qc_filter_umi(FL_316_8724)
FL_319_8724 <- iqr_qc_filter_umi(FL_319_8724)
FL_320_72924 <- iqr_qc_filter_umi(FL_320_72924)




seurat_objects <- list(
  FL_207_091423 = FL_207_091423,
  FL_31_091423 = FL_31_091423,
  FL_296_091423 = FL_296_091423,
  FL_8_091423R = FL_8_091423R,
  FL_171_091423R = FL_171_091423R,
  FL_12_091823 = FL_12_091823,
  FL_6_091823 = FL_6_091823,
  FL_625_091823 = FL_625_091823,
  FL_199_091823 = FL_199_091823,
  FL_9_091923 = FL_9_091923,
  FL_17_092023 = FL_17_092023,
  FL_45_092023 = FL_45_092023,
  FL_2_092123 = FL_2_092123,
  FL_5_092123 = FL_5_092123,
  FL_18_092123 = FL_18_092123,
  FL_296_092123 = FL_296_092123,
  FL_95_092123 = FL_95_092123,
  FL_62_092123 = FL_62_092123,
  FL_22NA_72624 = FL_22NA_72624,
  FL_24_8724 = FL_24_8724,
  FL_24_72924 = FL_24_72924,
  FL_25NA_72624 = FL_25NA_72624,
  FL_28_72524 = FL_28_72524,
  FL_28_72624 = FL_28_72624,
  FL_33_72924 = FL_33_72924,
  FL_53_72924 = FL_53_72924,
  FL_59_72624 = FL_59_72624,
  FL_62_72524 = FL_62_72524,
  FL_62_72924 = FL_62_72924,
  FL_167_72624 = FL_167_72624,
  FL_283_72924 = FL_283_72924,
  FL_313_72924 = FL_313_72924,
  FL_314_8724 = FL_314_8724,
  FL_315_8724 = FL_315_8724,
  FL_316_72524 = FL_316_72524,
  FL_321_72624 = FL_321_72624,
  FL_323_8724 = FL_323_8724,
  FL_154_72524 = FL_154_72524,
  FL_316_8724 = FL_316_8724,
  FL_319_8724 = FL_319_8724,
  FL_320_72924 = FL_320_72924
)

# Create a list to store the UMI count DataFrames
umi_count_dfs <- list()

# Loop through each Seurat object, extract UMI counts, and store in a DataFrame
for (sample_name in names(seurat_objects)) {
  umi_counts <- seurat_objects[[sample_name]]@meta.data$nCount_RNA
  
  # Create a DataFrame with cell names as row names (index) and UMI counts as a column
  umi_df <- data.frame(
    Cell = rownames(seurat_objects[[sample_name]]@meta.data),
    UMI_Counts = umi_counts
  )
  
  # Store in the list
  umi_count_dfs[[sample_name]] <- umi_df
}

# Example: View UMI count DataFrame for one sample
head(umi_count_dfs$FL_207_091423)



# 4-2) Combine all UMI counts for Kruskal-Wallis Test
# Combine all UMI counts into a single DataFrame for Kruskal-Wallis test
all_umi_counts <- bind_rows(lapply(names(umi_count_dfs), function(sample_name) {
  umi_count_dfs[[sample_name]] %>%
    mutate(Sample = sample_name)  # Add sample name as a column
}))

# View the combined dataset
head(all_umi_counts)

# Save
write.csv(all_umi_counts, "output/round_2/output_4_step5_qc_of_round1AND2/post_IQR2p5_UMI_counts_all_samples.csv")

# 4-3) Perform Kruskan-Wallis Test
kruskal_test <- kruskal.test(UMI_Counts ~ Sample, data = all_umi_counts)
# Print the test result
print(kruskal_test)
pairwise.wilcox.test(all_umi_counts$UMI_Counts, all_umi_counts$Sample, p.adjust.method = "bonferroni")

# Convert the matrix into a long format for ggplot
p_matrix <- pairwise.wilcox.test(all_umi_counts$UMI_Counts, all_umi_counts$Sample, p.adjust.method = "bonferroni")$p.value
p_matrix[is.na(p_matrix)] <- 1  # Replace NA values with 1 for visualization
p_long <- melt(as.matrix(p_matrix))

# Plot as heatmap
ggplot(p_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", na.value = "white") +
  labs(title = "Wilcoxon Test Pairwise P-values",
       x = "Sample", y = "Sample", fill = "p-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Benjamini-Yekutieli
pairwise.wilcox.test(all_umi_counts$UMI_Counts, all_umi_counts$Sample, p.adjust.method = "BY")

# Convert the matrix into a long format for ggplot
p_matrix <- pairwise.wilcox.test(all_umi_counts$UMI_Counts, all_umi_counts$Sample, p.adjust.method = "BY")$p.value
p_matrix[is.na(p_matrix)] <- 1  # Replace NA values with 1 for visualization
p_long <- melt(as.matrix(p_matrix))

# Plot as heatmap
ggplot(p_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", na.value = "white") +
  labs(title = "Wilcoxon Test Pairwise P-values",
       x = "Sample", y = "Sample", fill = "p-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add log-transformed UMI counts to the dataframe
all_umi_counts$log_UMI_Counts <- log1p(all_umi_counts$UMI_Counts)

# View the updated dataframe
head(all_umi_counts)


# Run pairwise Wilcoxon test with Holm correction on log-transformed data
p_matrix_log_holm <- pairwise.wilcox.test(all_umi_counts$log_UMI_Counts, all_umi_counts$Sample, p.adjust.method = "hochberg")$p.value

# Replace NA values with 1 for visualization
p_matrix_log_holm[is.na(p_matrix_log_holm)] <- 1

# View the updated p-value matrix
print(p_matrix_log_holm)

# Convert the matrix into a long format for ggplot
p_matrix <- pairwise.wilcox.test(all_umi_counts$log_UMI_Counts, all_umi_counts$Sample, p.adjust.method = "BY")$p.value
p_matrix[is.na(p_matrix)] <- 1  # Replace NA values with 1 for visualization
p_long <- melt(as.matrix(p_matrix))
print(p_matrix)
# Plot as heatmap
ggplot(p_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "white", na.value = "white") +
  labs(title = "Wilcoxon Test Pairwise P-values",
       x = "Sample", y = "Sample", fill = "p-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(all_umi_counts, aes(x = UMI_Counts, color = Sample)) +
  geom_density(alpha = 0.7, linewidth = 1) +
  theme_minimal() +
  labs(title = "Density Plot of UMI Counts Across Samples",
       x = "UMI Counts per Cell",
       y = "Density")

ggplot(all_umi_counts, aes(x = log_UMI_Counts, color = Sample)) +
  geom_density(alpha = 0.7, linewidth = 1) +
  theme_minimal() +
  labs(title = "Density Plot of UMI Counts Across Samples",
       x = "logUMI Counts per Cell",
       y = "Density")


######################
seurat_objects <- list(
  FL_22NA_72624 = FL_22NA_72624,
  FL_313_72924 = FL_313_72924,
  FL_316_72524 = FL_316_72524,
  FL_316_8724 = FL_316_8724,
  FL_323_8724 = FL_323_8724,
  FL_320_72924 = FL_320_72924
)
# Create a list to store the UMI count DataFrames
umi_count_dfs <- list()

# Loop through each Seurat object, extract UMI counts, and store in a DataFrame
for (sample_name in names(seurat_objects)) {
  umi_counts <- seurat_objects[[sample_name]]@meta.data$nCount_RNA
  
  # Create a DataFrame with cell names as row names (index) and UMI counts as a column
  umi_df <- data.frame(
    Cell = rownames(seurat_objects[[sample_name]]@meta.data),
    UMI_Counts = umi_counts
  )
  
  # Store in the list
  umi_count_dfs[[sample_name]] <- umi_df
}

# Example: View UMI count DataFrame for one sample
head(umi_count_dfs$FL_22NA_72624)

# Combine all UMI counts into a single DataFrame for Kruskal-Wallis test
all_umi_counts <- bind_rows(lapply(names(umi_count_dfs), function(sample_name) {
  umi_count_dfs[[sample_name]] %>%
    mutate(Sample = sample_name)  # Add sample name as a column
}))

# View the combined dataset
head(all_umi_counts)

# Add log-transformed UMI counts to the dataframe
all_umi_counts$log_UMI_Counts <- log1p(all_umi_counts$UMI_Counts)

# 8x10 l
ggplot(all_umi_counts, aes(x = log_UMI_Counts, color = Sample)) +
  geom_density(alpha = 0.7, linewidth = 1) +
  theme_minimal() +
  labs(title = "Density Plot of UMI Counts Across Samples",
       x = "logUMI Counts per Cell",
       y = "Density")

# In conclusion we should remove:
# FL_22NA_72624 : One of the two NFKB samples
# FL_320_72924: A SFTPB sample

#################################################################################





# 5. nCount RNA QC using: IQR x 2.5 qc

# Function to filter Seurat object based on IQR × 2.5 UMI threshold
iqr_qc_filter_umi <- function(seurat_object) {
  
  # Extract UMI counts
  umi_counts <- seurat_object$nCount_RNA
  
  # Compute IQR
  Q1 <- quantile(umi_counts, 0.25)
  Q3 <- quantile(umi_counts, 0.75)
  IQR <- Q3 - Q1
  high_umi_threshold <- Q3 + (2.5 * IQR)
  
  # Subset Seurat object to keep only valid cells
  filtered_seurat <- subset(seurat_object, subset = nCount_RNA <= high_umi_threshold)
  
  # Return the filtered Seurat object
  return(filtered_seurat)
}

# Round 1
FL_207_091423 <- iqr_qc_filter_umi(FL_207_091423)
FL_207_091423 # 7003 -> 6831

FL_31_091423 <- iqr_qc_filter_umi(FL_31_091423)
FL_31_091423 # 9089 -> 8911

FL_296_091423 <- iqr_qc_filter_umi(FL_296_091423)
FL_296_091423 # 9931 -> 9599

FL_8_091423R <- iqr_qc_filter_umi(FL_8_091423R)
FL_8_091423R # 7568 -> 7467

FL_171_091423R <- iqr_qc_filter_umi(FL_171_091423R)
FL_171_091423R # 14195 -> 13996

FL_12_091823 <- iqr_qc_filter_umi(FL_12_091823)
FL_12_091823 # 16122 -> 15899

FL_6_091823 <- iqr_qc_filter_umi(FL_6_091823)
FL_6_091823 # 9418 -> 9058

FL_625_091823 <- iqr_qc_filter_umi(FL_625_091823)
FL_625_091823 # 4012 -> 3793

FL_199_091823 <- iqr_qc_filter_umi(FL_199_091823)
FL_199_091823 # 7792 -> 7265

FL_9_091923 <- iqr_qc_filter_umi(FL_9_091923)
FL_9_091923 # 4620 -> 4468

FL_17_092023 <- iqr_qc_filter_umi(FL_17_092023)
FL_17_092023 # 7562 -> 7382

FL_45_092023 <- iqr_qc_filter_umi(FL_45_092023)
FL_45_092023 # 10408 -> 10253

FL_2_092123 <- iqr_qc_filter_umi(FL_2_092123)
FL_2_092123 # 11772 -> 11204

FL_5_092123 <- iqr_qc_filter_umi(FL_5_092123)
FL_5_092123 # 9895 -> 9538

FL_18_092123 <- iqr_qc_filter_umi(FL_18_092123)
FL_18_092123 # 16681 -> 15960

FL_296_092123 <- iqr_qc_filter_umi(FL_296_092123)
FL_296_092123 # 6081 -> 5741

FL_95_092123 <- iqr_qc_filter_umi(FL_95_092123)
FL_95_092123 # 8422 -> 7954

FL_62_092123 <- iqr_qc_filter_umi(FL_62_092123)
FL_62_092123 # 8819 -> 8250


# Round 2
FL_24_8724 <- iqr_qc_filter_umi(FL_24_8724)
FL_24_8724 # 7386 -> 7244

FL_24_72924 <- iqr_qc_filter_umi(FL_24_72924)
FL_24_72924 # 3921 -> 3773

FL_25NA_72624 <- iqr_qc_filter_umi(FL_25NA_72624)
FL_25NA_72624 # 7039 -> 6736

FL_28_72524 <- iqr_qc_filter_umi(FL_28_72524)
FL_28_72524 # 6921 -> 6830

FL_28_72624 <- iqr_qc_filter_umi(FL_28_72624)
FL_28_72624 # 10206 -> 10007

FL_33_72924 <- iqr_qc_filter_umi(FL_33_72924)
FL_33_72924 # 2031 -> 2011

FL_53_72924 <- iqr_qc_filter_umi(FL_53_72924)
FL_53_72924 # 2552 -> 2499

FL_59_72624 <- iqr_qc_filter_umi(FL_59_72624)
FL_59_72624 # 3121 -> 3009

FL_62_72524 <- iqr_qc_filter_umi(FL_62_72524)
FL_62_72524 # 8785 -> 8584

FL_62_72924 <- iqr_qc_filter_umi(FL_62_72924)
FL_62_72924 # 6763 -> 6652

FL_167_72624 <- iqr_qc_filter_umi(FL_167_72624)
FL_167_72624 # 4602 -> 4520

FL_283_72924 <- iqr_qc_filter_umi(FL_283_72924)
FL_283_72924 # 14902 -> 14177

FL_313_72924 <- iqr_qc_filter_umi(FL_313_72924)
FL_313_72924 # 6918 -> 6767

FL_314_8724 <- iqr_qc_filter_umi(FL_314_8724)
FL_314_8724 # 4316 -> 4099

FL_315_8724 <- iqr_qc_filter_umi(FL_315_8724)
FL_315_8724 # 9015 -> 8738

FL_316_72524 <- iqr_qc_filter_umi(FL_316_72524)
FL_316_72524 # 58201 -> 53555

FL_321_72624 <- iqr_qc_filter_umi(FL_321_72624)
FL_321_72624 # 3152 -> 3069

FL_323_8724 <- iqr_qc_filter_umi(FL_323_8724)
FL_323_8724 # 12821 -> 12174

FL_154_72524 <- iqr_qc_filter_umi(FL_154_72524)
FL_154_72524 # 7167 -> 6935

FL_316_8724 <- iqr_qc_filter_umi(FL_316_8724)
FL_316_8724 # 19526 -> 18757

FL_319_8724 <- iqr_qc_filter_umi(FL_319_8724)
FL_319_8724 # 9627 -> 9522





# 6. Save Data
saveRDS(FL_207_091423, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_207_091423.rds")
saveRDS(FL_31_091423, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_31_091423.rds")
saveRDS(FL_296_091423, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_296_091423.rds")
saveRDS(FL_8_091423R, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_8_091423R.rds")
saveRDS(FL_171_091423R, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_171_091423R.rds")
saveRDS(FL_12_091823, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_12_091823.rds")
saveRDS(FL_6_091823, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_6_091823.rds")
saveRDS(FL_625_091823, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_625_091823.rds")
saveRDS(FL_199_091823, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_199_091823.rds")
saveRDS(FL_9_091923, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_9_091923.rds")
saveRDS(FL_17_092023, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_17_092023.rds")
saveRDS(FL_45_092023, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_45_092023.rds")
saveRDS(FL_2_092123, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_2_092123.rds")
saveRDS(FL_5_092123, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_5_092123.rds")
saveRDS(FL_18_092123, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_18_092123.rds")
saveRDS(FL_296_092123, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_296_092123.rds")
saveRDS(FL_95_092123, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_95_092123.rds")
saveRDS(FL_62_092123, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_62_092123.rds")
saveRDS(FL_24_8724, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_24_8724.rds")
saveRDS(FL_24_72924, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_24_72924.rds")
saveRDS(FL_25NA_72624, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_25NA_72624.rds")
saveRDS(FL_28_72524, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_28_72524.rds")
saveRDS(FL_28_72624, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_28_72624.rds")
saveRDS(FL_33_72924, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_33_72924.rds")
saveRDS(FL_53_72924, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_53_72924.rds")
saveRDS(FL_59_72624, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_59_72624.rds")
saveRDS(FL_62_72524, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_62_72524.rds")
saveRDS(FL_62_72924, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_62_72924.rds")
saveRDS(FL_167_72624, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_167_72624.rds")
saveRDS(FL_283_72924, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_283_72924.rds")
saveRDS(FL_313_72924, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_313_72924.rds")
saveRDS(FL_314_8724, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_314_8724.rds")
saveRDS(FL_315_8724, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_315_8724.rds")
saveRDS(FL_316_72524, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_316_72524.rds")
saveRDS(FL_321_72624, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_321_72624.rds")
saveRDS(FL_323_8724, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_323_8724.rds")
saveRDS(FL_154_72524, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_154_72524.rds")
saveRDS(FL_316_8724, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_316_8724.rds")
saveRDS(FL_319_8724, file = "output/round_2/output_4_step5_qc_of_round1AND2/FL_319_8724.rds")





# 7. Plot
# SFTPB samples
p1 <- ggplot2::ggplot(FL_296_091423[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_296_091423") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p2 <- ggplot2::ggplot(FL_296_092123[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_296_092123") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p3 <- ggplot2::ggplot(FL_314_8724[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_314_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p4 <- ggplot2::ggplot(FL_316_72524[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_316_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p5 <- ggplot2::ggplot(FL_316_8724[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_316_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p6 <- ggplot2::ggplot(FL_319_8724[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_319_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")

# l 12.9 x 6
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3)


# Control Samples
p1 <- ggplot2::ggplot(FL_171_091423R[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_171_091423R") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p2 <- ggplot2::ggplot(FL_199_091823[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_199_091823") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p3 <- ggplot2::ggplot(FL_18_092123[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_18_092123") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p4 <- ggplot2::ggplot(FL_313_72924[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_313_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p5 <- ggplot2::ggplot(FL_315_8724[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_315_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")

# l 12.9 x 6
p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 3)


# ABCA3 samples
p1 <- ggplot2::ggplot(FL_207_091423[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_207_091423") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p2 <- ggplot2::ggplot(FL_625_091823[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_625_091823") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p3 <- ggplot2::ggplot(FL_45_092023[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_45_092023") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p4 <- ggplot2::ggplot(FL_95_092123[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_95_092123") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p5 <- ggplot2::ggplot(FL_62_092123[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_62_092123") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")

# l 12.9 x 6
p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 3)


# FARS2 and GRN samples
p1 <- ggplot2::ggplot(FL_62_72524[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_62_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p2 <- ggplot2::ggplot(FL_62_72924[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_62_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p3 <- ggplot2::ggplot(FL_9_091923[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_9_091923") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")

# l 12.9 x 3
p1 + p2 + p3 + plot_layout(ncol = 3)


# LRBA and NFKB1/NOD2 samples
p1 <- ggplot2::ggplot(FL_8_091423R[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_8_091423R") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p2 <- ggplot2::ggplot(FL_12_091823[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_12_091823") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p3 <- ggplot2::ggplot(FL_31_091423[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_31_091423") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")

# l 12.9 x 3
p1 + p2 + p3 + plot_layout(ncol = 3)


# NLRP12, SFTPC and SLC7A7 samples
p1 <- ggplot2::ggplot(FL_28_72524[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_28_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p2 <- ggplot2::ggplot(FL_28_72624[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_28_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p3 <- ggplot2::ggplot(FL_2_092123[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_2_092123") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p4 <- ggplot2::ggplot(FL_33_72924[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_33_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p5 <- ggplot2::ggplot(FL_24_72924[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_24_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p6 <- ggplot2::ggplot(FL_24_8724[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_24_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")

# l 25.3 x 3
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 6)


# SOCS1, STAT1, DICER1, FLNA samples
p1 <- ggplot2::ggplot(FL_6_091823[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_6_091823") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p2 <- ggplot2::ggplot(FL_167_72624[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_167_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p3 <- ggplot2::ggplot(FL_5_092123[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_5_092123") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p4 <- ggplot2::ggplot(FL_154_72524[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_154_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p5 <- ggplot2::ggplot(FL_321_72624[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_321_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p6 <- ggplot2::ggplot(FL_283_72924[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_283_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")

# l 25.3 x 3
p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 6)


# IKBKB, NPC2, PIK3CA, SMPD1, Unknown samples
p1 <- ggplot2::ggplot(FL_53_72924[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_53_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p2 <- ggplot2::ggplot(FL_17_092023[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_17_092023") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p3 <- ggplot2::ggplot(FL_59_72624[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_59_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p4 <- ggplot2::ggplot(FL_323_8724[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_323_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
p5 <- ggplot2::ggplot(FL_25NA_72624[[]]) +
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_25NA_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")

# l 22 x 3
p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 6)








# p 4.3x12
p1 + p2 + p3 + p4 + plot_layout(ncol = 1)

