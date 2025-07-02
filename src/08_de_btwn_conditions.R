# DE between conditions in cell types of interest.
# Using preprocessed file from Han sent 6/13/25 czi_combined_step10_clustering....


# Load libraries ----
library(dplyr)
# install.packages('Seurat')
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
# remotes::install_github('satijalab/seurat-wrappers')
library(Matrix)
library(SeuratWrappers)
library(SingleCellExperiment)
library(scCustomize)
library(data.table)



# 1. Import Data ----
input_file <- "output/czi_combined_step10_clustering_and_UMAP_harmony_nCountRNA.rds"
harmony_adjusted <- readRDS(input_file)
harmony_adjusted # 341046 cells
Idents(harmony_adjusted) <- "seurat_clusters"

# 2. Process Data ---
## 2-1) Subset the cell types of interest ----
# Macrophages and AT2 Cells
# Seurat clusters are not yet finalized
harmony_adjusted_at2 <- subset(harmony_adjusted, idents = c("0", "1"))
harmony_adjusted_macs <- subset(harmony_adjusted, idents = c("8", "9", "25", "27", "28"))
rm(harmony_adjusted) # Save RAM space
Idents(harmony_adjusted_at2) <- "genetics"
Idents(harmony_adjusted_macs) <- "genetics"
harmony_adjusted_at2 <- subset(harmony_adjusted_at2, idents = c("NPC2", "Control"))
harmony_adjusted_macs <- subset(harmony_adjusted_macs, idents = c("NPC2", "Control"))
# Visualize subsets
DimPlot(harmony_adjusted_at2, reduction = "umap", raster = FALSE)
DimPlot(harmony_adjusted_macs, reduction = "umap", raster = FALSE)

# 3. DE between npc2 vs non-pneumonia control for each cell type of interest ----

## 3.1) AT2 ----
de_genes_at2 <- FindMarkers(harmony_adjusted_at2,
  ident.1 = "NPC2",
  ident.2 = "Control",
  recorrect_umi = FALSE,
  only.pos = FALSE,
  min.pct = 0.01,
  logfc.threshold = 0.1,
  test.use = "wilcox"
) # Use Wilcoxon rank sum test

de_genes_macs <- FindMarkers(harmony_adjusted_macs,
  ident.1 = "NPC2",
  ident.2 = "Control",
  recorrect_umi = FALSE,
  only.pos = FALSE,
  min.pct = 0.01,
  logfc.threshold = 0.1,
  test.use = "wilcox"
) # Use Wilcoxon rank sum test


# 4 Write the results to a CSV file ----
fwrite(de_genes_at2,
  file = "output/08_de_at2_vs_ctrl.csv",
  row.names = TRUE
)

write.csv(de_genes_macs,
  file = "output/08_de_macs_vs_ctrl.csv",
  row.names = TRUE
)
