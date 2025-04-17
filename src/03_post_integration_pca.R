# Script to generate pca
set.seed(123)

# Load libraries ----
library(dplyr)
#install.packages('Seurat')
library(Seurat)
library(patchwork)
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
library(data.table)
library(plyr)

input_file_name <- "output/02_integrated_biopsy_samples_metadata_updated.rds"
output_file_name <- "output/03_post_integration_pca.rds"

# Read data from integration ----
czi_combined <- readRDS(input_file_name)

## Convert active assay to 'integrated' for scaling ----
DefaultAssay(czi_combined) <- "integrated"
czi_combined
unique(czi_combined$updated_sample_id)

# Scale data ----
# During Scaling we can regress out unwanted effects shown in the code below.
# czi_combined <- ScaleData(czi_combined, verbose = TRUE, vars.to.regress = "percent.mt", "Sex", "Age")
# However during the first run we want to show the purest Data as possible, and will
# not regress out variables.
czi_combined <- ScaleData(czi_combined, verbose = TRUE)

# PCA -----
# Note we already have the VariableFeatures for the 'integrated' assay, in case we would
# want to use a another assay which does not have VariableFeatures, 
# run the 'FindVariableFeatures' function
czi_combined <- RunPCA(czi_combined, features = VariableFeatures(czi_combined))

## Elbow plots for Dimension Selection (4x7 Landscape) ----
ElbowPlot(czi_combined, ndims = 50)

## Plot DimHeatmap for 1:10 (10x15 Portrait) ----
DimHeatmap(czi_combined, dims = 1:10, cells = 1000, balanced = TRUE)
DimHeatmap(czi_combined, dims = 11:20, cells = 1000, balanced = TRUE)
DimHeatmap(czi_combined, dims = 21:30, cells = 1000, balanced = TRUE)



# Save Data ----

saveRDS(czi_combined, file = output_file_name)

