# Step2 QC for each Seurat Object round 2





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





# 3. Plot
## 3-1) Vln plots through "nFeature_RNA", "nCount_RNA", "percent.mt"
# 6x11 L
VlnPlot(FL_323_8724, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = "lightblue")



## 3-2) Vln plots of "percent.mt" and Scattered plot, pre-QC
### 3-2-1) FL_22NA_72624
# Vln
p1 <- VlnPlot(FL_22NA_72624, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_22NA_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_22NA_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-2) FL_24_8724
# Vln
p1 <- VlnPlot(FL_24_8724, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_24_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_24_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-3) FL_24_72924
# Vln
p1 <- VlnPlot(FL_24_72924, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_24_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_24_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-4) FL_25NA_72624
# Vln
p1 <- VlnPlot(FL_25NA_72624, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_25NA_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_25NA_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-5) FL_28_72524
# Vln
p1 <- VlnPlot(FL_28_72524, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_28_72524[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_28_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-6) FL_28_72624
# Vln
p1 <- VlnPlot(FL_28_72624, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_28_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_28_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-7) FL_33_72924
# Vln
p1 <- VlnPlot(FL_33_72924, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_33_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_33_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-8) FL_53_72924
# Vln
p1 <- VlnPlot(FL_53_72924, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_53_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_53_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-9) FL_59_72624
# Vln
p1 <- VlnPlot(FL_59_72624, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_59_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_59_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-10) FL_62_72524
# Vln
p1 <- VlnPlot(FL_62_72524, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_62_72524[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_62_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-11) FL_62_72924
# Vln
p1 <- VlnPlot(FL_62_72924, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_62_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_62_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-12) FL_167_72624
# Vln
p1 <- VlnPlot(FL_167_72624, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_167_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_167_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-13) FL_283_72924
# Vln
p1 <- VlnPlot(FL_283_72924, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_283_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_283_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-14) FL_313_72924
# Vln
p1 <- VlnPlot(FL_313_72924, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_313_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_313_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-15) FL_314_8724
# Vln
p1 <- VlnPlot(FL_314_8724, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_314_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_314_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-16) FL_315_8724
# Vln
p1 <- VlnPlot(FL_315_8724, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_315_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_315_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-17) FL_316_72524
# Vln
p1 <- VlnPlot(FL_316_72524, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_316_72524[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_316_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-18) FL_321_72624
# Vln
p1 <- VlnPlot(FL_321_72624, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_321_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_321_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-19) FL_323_8724
# Vln
p1 <- VlnPlot(FL_323_8724, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_323_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_323_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-20) FL_154_72524
# Vln
p1 <- VlnPlot(FL_154_72524, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_154_72524[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_154_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-21) FL_316_8724
# Vln
p1 <- VlnPlot(FL_316_8724, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_316_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_316_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-22) FL_319_8724
# Vln
p1 <- VlnPlot(FL_319_8724, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_319_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_319_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 3-2-23) FL_320_72924
# Vln
p1 <- VlnPlot(FL_320_72924, features = c("percent.mt"), y.max = 6, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_320_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_320_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)





# 4. Set 'percent.mt' threshold
## 4-1) percent.mt < 1%
x <- 1
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
# Checking
FL_22NA_72624 # 52681 nuclei
FL_24_8724 # 7386 nuclei
FL_24_72924 # 3921 nuclei
FL_25NA_72624 # 7039 nuclei
FL_28_72524 # 6921 nuclei
FL_28_72624 # 10206 nuclei
FL_33_72924 # 2031 nuclei
FL_53_72924 # 2552 nuclei
FL_59_72624 # 3121 nuclei
FL_62_72524 # 8785 nuclei
FL_62_72924 # 6763 nuclei
FL_167_72624 # 4602 nuclei
FL_283_72924 # 14902 nuclei
FL_313_72924 # 6918 nuclei
FL_314_8724 # 4316 nuclei
FL_315_8724 # 9015 nuclei
FL_316_72524 # 58201 nuclei
FL_321_72624 # 3152 nuclei
FL_323_8724 # 12821 nuclei
FL_154_72524 # 7167 nuclei
FL_316_8724 # 19526 nuclei
FL_319_8724 # 9627 nuclei
FL_320_72924 # 18996 nuclei



## 4-2) Vln plots of "percent.mt" and Scattered plot, post-QC
### 4-2-1) FL_22NA_72624
# Vln
p1 <- VlnPlot(FL_22NA_72624, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_22NA_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_22NA_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-2) FL_24_8724
# Vln
p1 <- VlnPlot(FL_24_8724, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_24_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_24_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-3) FL_24_72924
# Vln
p1 <- VlnPlot(FL_24_72924, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_24_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_24_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-4) FL_25NA_72624
# Vln
p1 <- VlnPlot(FL_25NA_72624, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_25NA_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_25NA_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-5) FL_28_72524
# Vln
p1 <- VlnPlot(FL_28_72524, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_28_72524[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_28_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-6) FL_28_72624
# Vln
p1 <- VlnPlot(FL_28_72624, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_28_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_28_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-7) FL_33_72924
# Vln
p1 <- VlnPlot(FL_33_72924, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_33_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_33_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-8) FL_53_72924
# Vln
p1 <- VlnPlot(FL_53_72924, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_53_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_53_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-9) FL_59_72624
# Vln
p1 <- VlnPlot(FL_59_72624, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_59_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_59_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-10) FL_62_72524
# Vln
p1 <- VlnPlot(FL_62_72524, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_62_72524[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_62_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-11) FL_62_72924
# Vln
p1 <- VlnPlot(FL_62_72924, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_62_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_62_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-12) FL_167_72624
# Vln
p1 <- VlnPlot(FL_167_72624, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_167_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_167_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-13) FL_283_72924
# Vln
p1 <- VlnPlot(FL_283_72924, features = c("percent.mt"),  y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_283_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_283_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-14) FL_313_72924
# Vln
p1 <- VlnPlot(FL_313_72924, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_313_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_313_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-15) FL_314_8724
# Vln
p1 <- VlnPlot(FL_314_8724, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_314_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_314_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-16) FL_315_8724
# Vln
p1 <- VlnPlot(FL_315_8724, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_315_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_315_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-17) FL_316_72524
# Vln
p1 <- VlnPlot(FL_316_72524, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_316_72524[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_316_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-18) FL_321_72624
# Vln
p1 <- VlnPlot(FL_321_72624, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_321_72624[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_321_72624") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-19) FL_323_8724
# Vln
p1 <- VlnPlot(FL_323_8724, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_323_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_323_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-20) FL_154_72524
# Vln
p1 <- VlnPlot(FL_154_72524, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_154_72524[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_154_72524") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-21) FL_316_8724
# Vln
p1 <- VlnPlot(FL_316_8724, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_316_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_316_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-22) FL_319_8724
# Vln
p1 <- VlnPlot(FL_319_8724, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_319_8724[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_319_8724") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)

### 4-2-23) FL_320_72924
# Vln
p1 <- VlnPlot(FL_320_72924, features = c("percent.mt"), y.max = 1, ncol = 1, cols = "#ff7f00")
p1 <- p1 + theme(legend.position="none")
# Scattered plot
p2 <- ggplot2::ggplot(FL_320_72924[[]]) + 
  ggplot2::geom_point(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), alpha = 0.7) +
  scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
  ggplot2::ggtitle("FL_320_72924") +
  ggplot2::labs(x = "nCount_RNA", y = "nFeature_RNA")
# 8.8 x 5
p1 + p2 + plot_layout(ncol = 2)





# 5. Count UMI Distribution with IQR cutoff

# Function to calculate IQR threshold and generate a box plot
plot_umi_boxplot <- function(seurat_object, sample_name) {
  
  # Extract UMI counts
  umi_counts <- seurat_object$nCount_RNA
  
  # Compute IQR
  Q1 <- quantile(umi_counts, 0.25)
  Q3 <- quantile(umi_counts, 0.75)
  IQR <- Q3 - Q1
  high_umi_threshold <- Q3 + (1.5 * IQR)
  
  # Count cells that will be removed
  total_cells <- length(umi_counts)
  removed_cells <- sum(umi_counts > high_umi_threshold)
  kept_cells <- total_cells - removed_cells
  
  # Create a dataframe for plotting
  df <- data.frame(Sample = sample_name, UMI_Counts = umi_counts)
  
  # Generate Box Plot
  plot <- ggplot(df, aes(x = Sample, y = UMI_Counts)) +
    geom_boxplot(fill = "skyblue", alpha = 0.6, outlier.shape = NA) + 
    geom_jitter(aes(color = UMI_Counts > high_umi_threshold), width = 0.2, size = 0.8, alpha = 0.5) +
    geom_hline(yintercept = high_umi_threshold, color = "red", linetype = "dashed", size = 1) +
    labs(title = paste("UMI Count Distribution -", sample_name),
         subtitle = paste("Total Cells:", total_cells, "| Removed:", removed_cells, "| Kept:", kept_cells),
         x = "Sample",
         y = "UMI Counts per Cell") +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), labels = c("Valid", "Outlier")) +
    theme_minimal()
  
  return(plot)
}


# Function to calculate IQR threshold and generate a violin plot
plot_umi_violin <- function(seurat_object, sample_name) {
  
  # Extract UMI counts
  umi_counts <- seurat_object$nCount_RNA
  
  # Compute IQR
  Q1 <- quantile(umi_counts, 0.25)
  Q3 <- quantile(umi_counts, 0.75)
  IQR <- Q3 - Q1
  high_umi_threshold <- Q3 + (1.5 * IQR)
  
  # Count cells that will be removed
  total_cells <- length(umi_counts)
  removed_cells <- sum(umi_counts > high_umi_threshold)
  kept_cells <- total_cells - removed_cells
  
  # Create a dataframe for plotting
  df <- data.frame(Sample = sample_name, UMI_Counts = umi_counts)
  
  # Generate Violin Plot
  plot <- ggplot(df, aes(x = Sample, y = UMI_Counts)) +
    geom_violin(fill = "skyblue", alpha = 0.6) + 
    geom_jitter(aes(color = UMI_Counts > high_umi_threshold), width = 0.2, size = 0.8, alpha = 0.5) +
    geom_hline(yintercept = high_umi_threshold, color = "red", linetype = "dashed", size = 1) +
    labs(title = paste("UMI Count Distribution -", sample_name),
         subtitle = paste("Total Cells:", total_cells, "| Removed:", removed_cells, "| Kept:", kept_cells),
         x = "Sample",
         y = "UMI Counts per Cell") +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), labels = c("Valid", "Outlier")) +
    theme_minimal()
  
  return(plot)
}

plot_umi_violin(FL_22NA_72624, "FL_22NA_72624") # 48846
plot_umi_violin(FL_24_8724, "FL_24_8724") # 7015
plot_umi_violin(FL_24_72924, "FL_24_72924") # 3616
plot_umi_violin(FL_25NA_72624, "FL_25NA_72624") # 6520
plot_umi_violin(FL_28_72524, "FL_28_72524") # 6626
plot_umi_violin(FL_28_72624, "FL_28_72624") # 9760
plot_umi_violin(FL_33_72924, "FL_33_72924") # 1956
plot_umi_violin(FL_53_72924, "FL_53_72924") # 2429
plot_umi_violin(FL_59_72624, "FL_59_72624") # 2893
plot_umi_violin(FL_62_72524, "FL_62_72524") # 8353
plot_umi_violin(FL_62_72924, "FL_62_72924") # 6469
plot_umi_violin(FL_167_72624, "FL_167_72624") # 4392
plot_umi_violin(FL_283_72924, "FL_283_72924") # 13703
plot_umi_violin(FL_313_72924, "FL_313_72924") # 6534
plot_umi_violin(FL_314_8724, "FL_314_8724") # 3955
plot_umi_violin(FL_315_8724, "FL_315_8724") # 8385
plot_umi_violin(FL_316_72524, "FL_316_72524") # 52044
plot_umi_violin(FL_321_72624, "FL_321_72624") # 2976
plot_umi_violin(FL_323_8724, "FL_323_8724") # 11683
plot_umi_violin(FL_154_72524, "FL_154_72524") # 6673
plot_umi_violin(FL_316_8724, "FL_316_8724") # 17979
plot_umi_violin(FL_319_8724, "FL_319_8724") # 9286
plot_umi_violin(FL_320_72924, "FL_320_72924") # 16548



# ** More Lenient
# Function to calculate IQR threshold and generate a violin plot
plot_umi_violin <- function(seurat_object, sample_name) {
  
  # Extract UMI counts
  umi_counts <- seurat_object$nCount_RNA
  
  # Compute IQR
  Q1 <- quantile(umi_counts, 0.25)
  Q3 <- quantile(umi_counts, 0.75)
  IQR <- Q3 - Q1
  high_umi_threshold <- Q3 + (2.5 * IQR)
  
  # Count cells that will be removed
  total_cells <- length(umi_counts)
  removed_cells <- sum(umi_counts > high_umi_threshold)
  kept_cells <- total_cells - removed_cells
  
  # Create a dataframe for plotting
  df <- data.frame(Sample = sample_name, UMI_Counts = umi_counts)
  
  # Generate Violin Plot
  plot <- ggplot(df, aes(x = Sample, y = UMI_Counts)) +
    geom_violin(fill = "#E0FFFF", alpha = 0.6) + 
    geom_jitter(aes(color = UMI_Counts > high_umi_threshold), width = 0.2, size = 0.8, alpha = 0.5) +
    geom_hline(yintercept = high_umi_threshold, color = "red", linetype = "dashed", size = 1) +
    labs(title = paste("UMI Count Distribution -", sample_name),
         subtitle = paste("Total Cells:", total_cells, "| Removed:", removed_cells, "| Kept:", kept_cells),
         x = "Sample",
         y = "UMI Counts per Cell") +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), labels = c("Valid", "Outlier")) +
    theme_minimal()
  
  return(plot)
}

plot_umi_violin(FL_22NA_72624, "FL_22NA_72624") # 50581
plot_umi_violin(FL_24_8724, "FL_24_8724") # 7244
plot_umi_violin(FL_24_72924, "FL_24_72924") # 3773
plot_umi_violin(FL_25NA_72624, "FL_25NA_72624") # 6736
plot_umi_violin(FL_28_72524, "FL_28_72524") # 6830
plot_umi_violin(FL_28_72624, "FL_28_72624") # 10007
plot_umi_violin(FL_33_72924, "FL_33_72924") # 2011
plot_umi_violin(FL_53_72924, "FL_53_72924") # 2499
plot_umi_violin(FL_59_72624, "FL_59_72624") # 3009
plot_umi_violin(FL_62_72524, "FL_62_72524") # 8584
plot_umi_violin(FL_62_72924, "FL_62_72924") # 6652
plot_umi_violin(FL_167_72624, "FL_167_72624") # 4520
plot_umi_violin(FL_283_72924, "FL_283_72924") # 14177
plot_umi_violin(FL_313_72924, "FL_313_72924") # 6767
plot_umi_violin(FL_314_8724, "FL_314_8724") # 4099
plot_umi_violin(FL_315_8724, "FL_315_8724") # 8738
plot_umi_violin(FL_316_72524, "FL_316_72524") # 53555
plot_umi_violin(FL_321_72624, "FL_321_72624") # 3069
plot_umi_violin(FL_323_8724, "FL_323_8724") # 12174
plot_umi_violin(FL_154_72524, "FL_154_72524") # 6935
plot_umi_violin(FL_316_8724, "FL_316_8724") # 18757
plot_umi_violin(FL_319_8724, "FL_319_8724") # 9522
plot_umi_violin(FL_320_72924, "FL_320_72924") # 17180


