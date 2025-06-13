# Harmony adjustment

# Note the object that is imported is a integrated object (without 'SMPD1' sample).
# The object is post PCA & UMAP





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
library(pheatmap)
library(harmony)





# 1. Import Data
# Note the current object that is being Imported is the one that excludes the 'SMPD1' sample. 
# This is to note that the SMPD1 sample may cause batch effects and the quality of the data is concerning.
# The object is post PCA & UMAP
czi_combined <- readRDS(input_files)
czi_combined # 341046 nuclei





# 2. Run Harmony----
## 2-1) nCount_RNA----
### 2-1-1) Processing Data----
Idents(czi_combined) <- "seurat_clusters"
# Cut UMI counts into 10 bins
czi_combined$UMI_bin <- dplyr::ntile(czi_combined$nCount_RNA, 10)
# Make it a factor so Harmony knows it is a categorical
czi_combined$UMI_bin <- factor(czi_combined$UMI_bin)
# Check
table(czi_combined$UMI_bin)


### 2-1-2) RunHarmony----
czi_combined <- RunHarmony(czi_combined, 
                           group.by.vars = c("UMI_bin"))
czi_combined

# Find Neighbors
czi_combined <- FindNeighbors(czi_combined,
                              reduction = "harmony",
                              dims = 1:20)
# Find Clusters
czi_combined <- FindClusters(czi_combined,
                             resolution = 0.5)
# Run UMAP 
czi_combined <- RunUMAP(czi_combined,
                        reduction = "harmony",
                        dims = 1:20) # 'spread = 0.5' has been removed



### 2-1-3) Save Data----
saveRDS(czi_combined, file = output_file)



## 2-2) Source----
### 2-2-1) Import Data----
czi_combined <- readRDS(input_file)
czi_combined # 341046 nuclei
Idents(czi_combined) <- "seurat_clusters"
table(Idents(czi_combined))

### 2-2-2) Run Harmony----
czi_combined <- RunHarmony(czi_combined, 
                           group.by.vars = c("source"))
czi_combined

# Find Neighbors
czi_combined <- FindNeighbors(czi_combined,
                              reduction = "harmony",
                              dims = 1:20)
# Find Clusters
czi_combined <- FindClusters(czi_combined,
                             resolution = 0.5)
# Run UMAP 
czi_combined <- RunUMAP(czi_combined,
                        reduction = "harmony",
                        dims = 1:20) # 'spread = 0.5' has been removed

### 2-1-3) Save Data----
saveRDS(czi_combined, file = output_file)


