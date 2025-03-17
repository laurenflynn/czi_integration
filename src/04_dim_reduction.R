# Step 4 Non-linear Dimensional Reduction




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





# 1. Import Data ----
czi_combined <- readRDS("output/03_clustering_czi.rds")
czi_combined # 253193 nuclei





# 2. Dimensional Reduction using UMAP ----
czi_combined <- RunUMAP(czi_combined, dims = 1:20)
Idents(czi_combined) <- "seurat_clusters"
table(Idents(czi_combined))





# 3. Plot Landscape 8x7 ----
DimPlot(czi_combined, reduction = "umap", raster=FALSE, cols = "polychrome")
DimPlot(czi_combined, reduction = "umap", raster=FALSE, cols = "polychrome", label = TRUE)

##  Sample ID 10x7 ----
n_sample <- length(unique(czi_combined$orig.ident))
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
DimPlot(czi_combined, reduction = "umap", raster=FALSE, 
        cols = getPalette(n_sample), group.by = "orig.ident",
        split.by = "orig.ident")

##  Round 8x7 ----
DimPlot(czi_combined, reduction = "umap", raster=TRUE, 
        cols = c('1' = "#c9dee2", '2' = "#efcfd1"), group.by = "Round",
        split.by = "Round")

# To do: match genetics colors to han's document
##  Genetics 8x7 ----
DimPlot(czi_combined, reduction = "umap", raster=TRUE, 
        cols = c('SFTPB' = "#d77f80",
                 'Control' = "#b8d5a7",
                 'ABCA3' = "#869fbb",
                 'FARS2' = "#96bdda",
                 'GRN' = "#e5e80b",
                 'LRBA' = "#ac5379",
                 'NFKB1' = "#fff272",
                 'NLRP12' = "#edc9d5",
                 'SFTPC' = "#d5d6d1",
                 'SLC7A7' = "#72b7a2",
                 'SOCS1' = "#edc076",
                 'STAT1' = "#602d76",
                 'DICER1' = "#5e6a54",
                 'FLNA' = "#faf8cf",
                 'IKBKB' = "#b9ae95",
                 'NPC2' = "#c62f7c",
                 'PIK3CA' = "#a93337",
                 'SMPD1' = "#569ec9",
                 'Unknown' = "#201f48"), 
        group.by = "Genetics")





# 4. Saving Data ----
saveRDS(czi_combined, file = "output/04_czi_dim_reduction.rds")


