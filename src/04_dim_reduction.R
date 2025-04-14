# Step 4 Non-linear Dimensional Reduction

# Load libraries ----
library(dplyr)
# install.packages('Seurat')
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
library(RColorBrewer)


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
czi_combined <- readRDS(paste0("output/03_clustering_czi_dim_", dimensions, "_res_", resolution, ".rds"))
czi_combined 

# 2. Dimensional Reduction using UMAP ----
czi_combined <- RunUMAP(czi_combined, dims = 1:dimensions)
Idents(czi_combined) <- "seurat_clusters"
table(Idents(czi_combined))

# 3. Plot Landscape ----
DimPlot(czi_combined, reduction = "umap", raster = FALSE, cols = "polychrome")
ggsave(paste0("output/04_figures/unlabeled_umap_by_cluster_dim_", dimensions, "_res_", resolution, ".png"))
DimPlot(czi_combined, reduction = "umap", raster = FALSE, cols = "polychrome", label = TRUE)
ggsave(paste0("output/04_figures/labeled_umap_by_cluster_dim_", dimensions, "_res_", resolution, ".png"))

##  Sample ID ----
n_sample <- length(unique(czi_combined$orig.ident))
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
DimPlot(czi_combined,
  reduction = "umap", raster = FALSE,
  cols = getPalette(n_sample), group.by = "orig.ident",
  split.by = "orig.ident"
)
ggsave(paste0("output/04_figures/umap_by_id_dim_", dimensions, "_res_", resolution, ".png"))

##  Round ----
DimPlot(czi_combined,
  reduction = "umap", raster = TRUE,
  cols = c("1" = "#c9dee2", "2" = "#efcfd1"), group.by = "Round",
  split.by = "Round"
)
ggsave(paste0("output/04_figures/umap_by_round_dim_", dimensions, "_res_", resolution, ".png"))

##  Genetics ----
DimPlot(czi_combined,
  reduction = "umap", raster = TRUE,
  cols = c(
    "SFTPB" = "#d77f80",
    "Control" = "#b8d5a7",
    "ABCA3" = "#869fbb",
    "FARS2" = "#96bdda",
    "GRN" = "#e5e80b",
    "LRBA" = "#ac5379",
    "NFKB1" = "#fff272",
    "NLRP12" = "#edc9d5",
    "SFTPC" = "#d5d6d1",
    "SLC7A7" = "#72b7a2",
    "SOCS1" = "#edc076",
    "STAT1" = "#602d76",
    "DICER1" = "#5e6a54",
    "FLNA" = "#faf8cf",
    "IKBKB" = "#b9ae95",
    "NPC2" = "#c62f7c",
    "PIK3CA" = "#a93337",
    "SMPD1" = "#569ec9",
    "Unknown" = "#201f48"
  ),
  group.by = "Genetics"
)

ggsave(paste0("output/04_figures/umap_by_gene_dim_", dimensions, "_res_", resolution, ".png"))


# 4. Saving Data ----
saveRDS(czi_combined, file = paste0("output/04_czi_dim_reduction_dim_", dimensions, "_res_", resolution, ".rds"))
