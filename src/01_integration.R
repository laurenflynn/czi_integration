# Script to integrate single cell data using Seurat


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

# Set variables ----

input_data_dir <- "data/biopsy_samples/qc_pass"
output_file_name <- "output/integrated_biopsy_samples.rds"

## Read data ------

# Preallocate a list, read in data files into this list from data dir
seurat_objects <- vector(mode = "list", length(list.files(input_data_dir)))
i <- 1
for (file in list.files(input_data_dir)) {
  print(file)
  sample <- readRDS(paste0(input_data_dir, "/", file))
  seurat_objects[[i]] <- sample
  i <- i + 1
}

# Processing Data and Integration -----

for (obj in seurat_objects) {
  print(DefaultAssay(obj))
}

features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 1750)
seurat_objects <- PrepSCTIntegration(object.list = seurat_objects, anchor.features = features)
seurat_objects <- lapply(X = seurat_objects, FUN = RunPCA, features = features)


## Find anchors -----
anchors <- FindIntegrationAnchors(object.list = seurat_objects, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:20, reduction = "rpca", 
                                  k.anchor = 5, k.filter = 500)

## Integration ------
czi_combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:20)
czi_combined
unique(czi_combined$orig.ident)
rm(anchors)
rm(seurat_objects)
gc()

# Save results ------
saveRDS(czi_combined, file = output_file_name)

print(paste("Output saved to", output_file_name))



