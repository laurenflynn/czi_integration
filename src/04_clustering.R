# Clustering
# Libraries ------
library(dplyr)
library(Seurat)
library(patchwork)
library(hdf5r)
library(ggplot2)
library(Nebulosa)
library(sctransform)
library(monocle3)
library(magrittr)
library(R.utils)
library(Matrix)
library(SeuratWrappers)
library(SingleCellExperiment)

input_file_name <-"output/03_post_integration_pca.rds"

# Read in variables from sbatch -----------------
args <- commandArgs(trailingOnly = TRUE)
if ("--dim" %in% args) {
  dimensions <- as.numeric(args[which(args == "--dim") + 1])
} else {
  dimensions <- 20
}
print(paste("Dimensions:", dimensions))
if ("--res" %in% args) {
  resolution <- as.numeric(args[which(args == "--res") + 1])
} else {
  resolution <- 0.5
}
print(paste("Resolution:", resolution))

# 1. Import Data ----
czi_combined <- readRDS(input_file_name)
czi_combined 




# 2. Processing the Data ----
## 2-1) Find Neighbors ----
czi_combined <- FindNeighbors(czi_combined, dims = dimensions)
# The FindNeighbors function constructs a K-nearest neighbor graph using the first 30 principal components.
# This KNN graph captures the similarity in gene expression profiles among cells, considering the variance and patterns represented by these 30 PCs.
# Constructing this similarity network is essential for subsequent clustering, enabling the identification of cell communities with shared expression profiles."

## 2-2) Find Clusters -----
czi_combined <- FindClusters(czi_combined, resolution = resolution)
#"The FindClusters function performs community detection on the K-nearest neighbor graph constructed from the first 30 principal components.
# By applying the Louvain algorithm (or Leiden, based on preference), it partitions the graph into distinct clusters.
# Each cluster represents a community of cells with shared expression profiles, hinting at similar biological characteristics or cell types.
# The 'resolution' parameter can be adjusted to fine-tune the granularity of the clustering process."
# Higher resolution: This will result in more clusters, potentially identifying more specific or subtle cell populations.
# This can be useful if you suspect there are many distinct cell types or states in your dataset.
# However, setting the resolution too high might over-segment the data and produce many small, potentially less meaningful clusters.
# Lower resolution: This will produce fewer clusters by merging similar cell groups. It can be useful when you want to focus on broader categories or major cell populations.
# However, setting the resolution too low might overlook some important cell types or states.





# 3. Saving Data ----
saveRDS(czi_combined, file = paste0("output/04_clustering_czi_dim_", dimensions, "_res_", resolution, ".rds"))




