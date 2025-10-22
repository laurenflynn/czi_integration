# Step6 SCTransform and PCA for each Seurat Object





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






# 1. Import Data
# Round 1
FL_207_091423 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_207_091423.rds")
FL_31_091423 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_31_091423.rds")
FL_296_091423 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_296_091423.rds")
FL_8_091423R <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_8_091423R.rds")
FL_171_091423R <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_171_091423R.rds")
FL_12_091823 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_12_091823.rds")
FL_6_091823 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_6_091823.rds")
FL_625_091823 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_625_091823.rds")
FL_199_091823 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_199_091823.rds")
FL_9_091923 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_9_091923.rds")
FL_17_092023 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_17_092023.rds")
FL_45_092023 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_45_092023.rds")
FL_2_092123 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_2_092123.rds")
FL_5_092123 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_5_092123.rds")
FL_18_092123 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_18_092123.rds")
FL_296_092123 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_296_092123.rds")
FL_95_092123 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_95_092123.rds")
FL_62_092123 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_62_092123.rds")
# Round 2
FL_24_8724 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_24_8724.rds")
FL_24_72924 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_24_72924.rds")
FL_25NA_72624 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_25NA_72624.rds")
FL_28_72524 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_28_72524.rds")
FL_28_72624 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_28_72624.rds")
FL_33_72924 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_33_72924.rds")
FL_53_72924 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_53_72924.rds")
FL_59_72624 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_59_72624.rds")
FL_62_72524 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_62_72524.rds")
FL_62_72924 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_62_72924.rds")
FL_167_72624 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_167_72624.rds")
FL_283_72924 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_283_72924.rds")
FL_313_72924 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_313_72924.rds")
FL_314_8724 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_314_8724.rds")
FL_315_8724 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_315_8724.rds")
FL_316_72524 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_316_72524.rds")
FL_321_72624 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_321_72624.rds")
FL_323_8724 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_323_8724.rds")
FL_154_72524 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_154_72524.rds")
FL_316_8724 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_316_8724.rds")
FL_319_8724 <- readRDS("output/round_2/output_4_step5_qc_of_round1AND2/FL_319_8724.rds")





# 2. SCT
# List of Seurat objects
seurat_objects <- list(FL_207_091423, FL_31_091423, FL_296_091423, FL_8_091423R, FL_171_091423R,
                       FL_12_091823, FL_6_091823, FL_625_091823, FL_199_091823, FL_9_091923, 
                       FL_17_092023, FL_45_092023, FL_2_092123, FL_5_092123, FL_18_092123,
                       FL_296_092123, FL_95_092123, FL_62_092123,
                       FL_24_8724, FL_24_72924, FL_25NA_72624, FL_28_72524, FL_28_72624,
                       FL_33_72924, FL_53_72924, FL_59_72624, FL_62_72524, FL_62_72924,
                       FL_167_72624, FL_283_72924, FL_313_72924, FL_314_8724, FL_315_8724,
                       FL_316_72524, FL_321_72624, FL_323_8724, FL_154_72524, FL_316_8724,
                       FL_319_8724)

# Apply SCTransform to each Seurat object
for (i in seq_along(seurat_objects)) {
  # Get Seurat object
  seurat_object <- seurat_objects[[i]]
  
  # SCTransform
  seurat_object <- SCTransform(seurat_object, vars.to.regress = "percent.mt", verbose = FALSE)
  
  # Save the transformed object back into the list
  seurat_objects[[i]] <- seurat_object
}

# adding back in
FL_207_091423 <- seurat_objects[[1]]
FL_31_091423 <- seurat_objects[[2]]
FL_296_091423 <- seurat_objects[[3]]
FL_8_091423R <- seurat_objects[[4]]
FL_171_091423R <- seurat_objects[[5]]
FL_12_091823 <- seurat_objects[[6]]
FL_6_091823 <- seurat_objects[[7]]
FL_625_091823 <- seurat_objects[[8]]
FL_199_091823 <- seurat_objects[[9]]
FL_9_091923 <- seurat_objects[[10]]
FL_17_092023 <- seurat_objects[[11]]
FL_45_092023 <- seurat_objects[[12]]
FL_2_092123 <- seurat_objects[[13]]
FL_5_092123 <- seurat_objects[[14]]
FL_18_092123 <- seurat_objects[[15]]
FL_296_092123 <- seurat_objects[[16]]
FL_95_092123 <- seurat_objects[[17]]
FL_62_092123 <- seurat_objects[[18]]
FL_24_8724 <- seurat_objects[[19]]
FL_24_72924 <- seurat_objects[[20]]
FL_25NA_72624 <- seurat_objects[[21]]
FL_28_72524 <- seurat_objects[[22]]
FL_28_72624 <- seurat_objects[[23]]
FL_33_72924 <- seurat_objects[[24]]
FL_53_72924 <- seurat_objects[[25]]
FL_59_72624 <- seurat_objects[[26]]
FL_62_72524 <- seurat_objects[[27]]
FL_62_72924 <- seurat_objects[[28]]
FL_167_72624 <- seurat_objects[[29]]
FL_283_72924 <- seurat_objects[[30]]
FL_313_72924 <- seurat_objects[[31]]
FL_314_8724 <- seurat_objects[[32]]
FL_315_8724 <- seurat_objects[[33]]
FL_316_72524 <- seurat_objects[[34]]
FL_321_72624 <- seurat_objects[[35]]
FL_323_8724 <- seurat_objects[[36]]
FL_154_72524 <- seurat_objects[[37]]
FL_316_8724 <- seurat_objects[[38]]
FL_319_8724 <- seurat_objects[[39]]





# 3. PCA
FL_207_091423 <- RunPCA(FL_207_091423, features = VariableFeatures(object = FL_207_091423))
FL_31_091423 <- RunPCA(FL_31_091423, features = VariableFeatures(object = FL_31_091423))
FL_296_091423 <- RunPCA(FL_296_091423, features = VariableFeatures(object = FL_296_091423))
FL_8_091423R <- RunPCA(FL_8_091423R, features = VariableFeatures(object = FL_8_091423R))
FL_171_091423R <- RunPCA(FL_171_091423R, features = VariableFeatures(object = FL_171_091423R))
FL_12_091823 <- RunPCA(FL_12_091823, features = VariableFeatures(object = FL_12_091823))
FL_6_091823 <- RunPCA(FL_6_091823, features = VariableFeatures(object = FL_6_091823))
FL_625_091823 <- RunPCA(FL_625_091823, features = VariableFeatures(object = FL_625_091823))
FL_199_091823 <- RunPCA(FL_199_091823, features = VariableFeatures(object = FL_199_091823))
FL_9_091923 <- RunPCA(FL_9_091923, features = VariableFeatures(object = FL_9_091923))
FL_17_092023 <- RunPCA(FL_17_092023, features = VariableFeatures(object = FL_17_092023))
FL_45_092023 <- RunPCA(FL_45_092023, features = VariableFeatures(object = FL_45_092023))
FL_2_092123 <- RunPCA(FL_2_092123, features = VariableFeatures(object = FL_2_092123))
FL_5_092123 <- RunPCA(FL_5_092123, features = VariableFeatures(object = FL_5_092123))
FL_18_092123 <- RunPCA(FL_18_092123, features = VariableFeatures(object = FL_18_092123))
FL_296_092123 <- RunPCA(FL_296_092123, features = VariableFeatures(object = FL_296_092123))
FL_95_092123 <- RunPCA(FL_95_092123, features = VariableFeatures(object = FL_95_092123))
FL_62_092123 <- RunPCA(FL_62_092123, features = VariableFeatures(object = FL_62_092123))
FL_24_8724 <- RunPCA(FL_24_8724, features = VariableFeatures(object = FL_24_8724))
FL_24_72924 <- RunPCA(FL_24_72924, features = VariableFeatures(object = FL_24_72924))
FL_25NA_72624 <- RunPCA(FL_25NA_72624, features = VariableFeatures(object = FL_25NA_72624))
FL_28_72524 <- RunPCA(FL_28_72524, features = VariableFeatures(object = FL_28_72524))
FL_28_72624 <- RunPCA(FL_28_72624, features = VariableFeatures(object = FL_28_72624))
FL_33_72924 <- RunPCA(FL_33_72924, features = VariableFeatures(object = FL_33_72924))
FL_53_72924 <- RunPCA(FL_53_72924, features = VariableFeatures(object = FL_53_72924))
FL_59_72624 <- RunPCA(FL_59_72624, features = VariableFeatures(object = FL_59_72624))
FL_62_72524 <- RunPCA(FL_62_72524, features = VariableFeatures(object = FL_62_72524))
FL_62_72924 <- RunPCA(FL_62_72924, features = VariableFeatures(object = FL_62_72924))
FL_167_72624 <- RunPCA(FL_167_72624, features = VariableFeatures(object = FL_167_72624))
FL_283_72924 <- RunPCA(FL_283_72924, features = VariableFeatures(object = FL_283_72924))
FL_313_72924 <- RunPCA(FL_313_72924, features = VariableFeatures(object = FL_313_72924))
FL_314_8724 <- RunPCA(FL_314_8724, features = VariableFeatures(object = FL_314_8724))
FL_315_8724 <- RunPCA(FL_315_8724, features = VariableFeatures(object = FL_315_8724))
FL_316_72524 <- RunPCA(FL_316_72524, features = VariableFeatures(object = FL_316_72524))
FL_321_72624 <- RunPCA(FL_321_72624, features = VariableFeatures(object = FL_321_72624))
FL_323_8724 <- RunPCA(FL_323_8724, features = VariableFeatures(object = FL_323_8724))
FL_154_72524 <- RunPCA(FL_154_72524, features = VariableFeatures(object = FL_154_72524))
FL_316_8724 <- RunPCA(FL_316_8724, features = VariableFeatures(object = FL_316_8724))
FL_319_8724 <- RunPCA(FL_319_8724, features = VariableFeatures(object = FL_319_8724))





# 4. Plot DimHeatmap 5x5 landscape
# PC1
variable_dim <- 1
DimHeatmap(FL_207_091423, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_31_091423, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_296_091423, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_8_091423R, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_171_091423R, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_12_091823, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_6_091823, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_625_091823, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_199_091823, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_9_091923, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_17_092023, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_45_092023, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_2_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_5_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_18_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_296_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_95_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_62_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_24_8724, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_24_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_25NA_72624, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_28_72524, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_28_72624, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_33_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_53_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_59_72624, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_62_72524, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_62_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_167_72624, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_283_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_313_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_314_8724, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_315_8724, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_316_72524, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_321_72624, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_323_8724, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_154_72524, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_316_8724, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_319_8724, dims = variable_dim, cells = 100, balanced = TRUE)


# PC2
variable_dim <- 2
DimHeatmap(FL_207_091423, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_31_091423, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_296_091423, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_8_091423R, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_171_091423R, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_12_091823, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_6_091823, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_625_091823, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_199_091823, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_9_091923, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_17_092023, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_45_092023, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_2_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_5_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_18_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_296_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_95_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_62_092123, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_24_8724, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_24_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_25NA_72624, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_28_72524, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_28_72624, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_33_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_53_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_59_72624, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_62_72524, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_62_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_167_72624, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_283_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_313_72924, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_314_8724, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_315_8724, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_316_72524, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_321_72624, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_323_8724, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_154_72524, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_316_8724, dims = variable_dim, cells = 100, balanced = TRUE)
DimHeatmap(FL_319_8724, dims = variable_dim, cells = 100, balanced = TRUE)





# 5. Elbow plots for Dimension Selection 4x7
# SFTPB samples
p1 <- ElbowPlot(FL_296_091423, ndims = 50)
p2 <- ElbowPlot(FL_296_092123, ndims = 50)
p3 <- ElbowPlot(FL_314_8724, ndims = 50)
p4 <- ElbowPlot(FL_316_72524, ndims = 50)
p5 <- ElbowPlot(FL_316_8724, ndims = 50)
p6 <- ElbowPlot(FL_319_8724, ndims = 50)
# Control samples
p7 <- ElbowPlot(FL_171_091423R, ndims = 50)
p8 <- ElbowPlot(FL_199_091823, ndims = 50)
p9 <- ElbowPlot(FL_18_092123, ndims = 50)

# l 9 x 17
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(ncol = 3)


# Control samples
p1 <- ElbowPlot(FL_313_72924, ndims = 50)
p2 <- ElbowPlot(FL_315_8724, ndims = 50)
# ABCA3 samples
p3 <- ElbowPlot(FL_207_091423, ndims = 50)
p4 <- ElbowPlot(FL_625_091823, ndims = 50)
p5 <- ElbowPlot(FL_45_092023, ndims = 50)
p6 <- ElbowPlot(FL_95_092123, ndims = 50)
p7 <- ElbowPlot(FL_62_092123, ndims = 50)
# FARS2
p8 <- ElbowPlot(FL_62_72524, ndims = 50)
p9 <- ElbowPlot(FL_62_72924, ndims = 50)

# l 9 x 17
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(ncol = 3)


# GRN samples
p1 <- ElbowPlot(FL_9_091923, ndims = 50)
# LRBA samples
p2 <- ElbowPlot(FL_8_091423R, ndims = 50)
p3 <- ElbowPlot(FL_12_091823, ndims = 50)
# NFKB1/NOD2 samples
p4 <- ElbowPlot(FL_31_091423, ndims = 50)
# NLRP12 samples
p5 <- ElbowPlot(FL_28_72524, ndims = 50)
p6 <- ElbowPlot(FL_28_72624, ndims = 50)
# SFTPC samples
p7 <- ElbowPlot(FL_2_092123, ndims = 50)
p8 <- ElbowPlot(FL_33_72924, ndims = 50)
# SLC7A7 samples
p9 <- ElbowPlot(FL_24_72924, ndims = 50)

# l 9 x 17
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(ncol = 3)


# SLC7A7 samples
p1 <- ElbowPlot(FL_24_8724, ndims = 50)
# SOCS1 samples
p2 <- ElbowPlot(FL_6_091823, ndims = 50)
p3 <- ElbowPlot(FL_167_72624, ndims = 50)
# STAT1 samples
p4 <- ElbowPlot(FL_5_092123, ndims = 50)
p5 <- ElbowPlot(FL_154_72524, ndims = 50)
# DICER1 samples
p6 <- ElbowPlot(FL_321_72624, ndims = 50)
# FLNA sample
p7 <- ElbowPlot(FL_283_72924, ndims = 50)
# IKBKB sample
p8 <- ElbowPlot(FL_53_72924, ndims = 50)
# NPC2 samples
p9 <- ElbowPlot(FL_17_092023, ndims = 50)

# l 9 x 17
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(ncol = 3)


# PIK3CA samples
p1 <- ElbowPlot(FL_59_72624, ndims = 50)
# SMPD1 samples
p2 <- ElbowPlot(FL_323_8724, ndims = 50)
# unknown samples
p3 <- ElbowPlot(FL_25NA_72624, ndims = 50)

# l 9x17
p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + plot_layout(ncol = 3)

# Overall up till PC15




# 6. Save Data
saveRDS(FL_207_091423, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_207_091423.rds")
saveRDS(FL_31_091423, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_31_091423.rds")
saveRDS(FL_296_091423, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_296_091423.rds")
saveRDS(FL_8_091423R, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_8_091423R.rds")
saveRDS(FL_171_091423R, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_171_091423R.rds")
saveRDS(FL_12_091823, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_12_091823.rds")
saveRDS(FL_6_091823, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_6_091823.rds")
saveRDS(FL_625_091823, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_625_091823.rds")
saveRDS(FL_199_091823, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_199_091823.rds")
saveRDS(FL_9_091923, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_9_091923.rds")
saveRDS(FL_17_092023, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_17_092023.rds")
saveRDS(FL_45_092023, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_45_092023.rds")
saveRDS(FL_2_092123, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_2_092123.rds")
saveRDS(FL_5_092123, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_5_092123.rds")
saveRDS(FL_18_092123, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_18_092123.rds")
saveRDS(FL_296_092123, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_296_092123.rds")
saveRDS(FL_95_092123, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_95_092123.rds")
saveRDS(FL_62_092123, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_62_092123.rds")
saveRDS(FL_24_8724, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_24_8724.rds")
saveRDS(FL_24_72924, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_24_72924.rds")
saveRDS(FL_25NA_72624, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_25NA_72624.rds")
saveRDS(FL_28_72524, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_28_72524.rds")
saveRDS(FL_28_72624, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_28_72624.rds")
saveRDS(FL_33_72924, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_33_72924.rds")
saveRDS(FL_53_72924, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_53_72924.rds")
saveRDS(FL_59_72624, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_59_72624.rds")
saveRDS(FL_62_72524, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_62_72524.rds")
saveRDS(FL_62_72924, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_62_72924.rds")
saveRDS(FL_167_72624, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_167_72624.rds")
saveRDS(FL_283_72924, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_283_72924.rds")
saveRDS(FL_313_72924, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_313_72924.rds")
saveRDS(FL_314_8724, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_314_8724.rds")
saveRDS(FL_315_8724, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_315_8724.rds")
saveRDS(FL_316_72524, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_316_72524.rds")
saveRDS(FL_321_72624, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_321_72624.rds")
saveRDS(FL_323_8724, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_323_8724.rds")
saveRDS(FL_154_72524, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_154_72524.rds")
saveRDS(FL_316_8724, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_316_8724.rds")
saveRDS(FL_319_8724, file = "output/round_2/output_4_step6_sctNormalization_and_PCA_for_each_object/FL_319_8724.rds")




