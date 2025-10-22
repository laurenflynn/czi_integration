# Step12 Cell type annotation

# Note: We will be using Seurat v5 'Reference Mapping'.
# The reference will be the HLCA





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
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)





# 1. Set up the Reference ----
# 1-1) Import Data ----
hlca <- readRDS("output/output_4_using_SingleR_IdentificationOfCellTypes_step9/reference_data/local.rds")
hlca # 584944 cells
# Trying to find where gene symbol is stored
head(rownames(hlca), 10)
hlca[["RNA"]]
head(rownames(hlca[["RNA"]]@meta.features))
# Found that the gene symbols are stored in the following
head(hlca[["RNA"]]@meta.features$feature_name, 10)  




# 1-2) Process the Data ----
## 1-2-1) Extract feature metadata and cleaned gene symbols ----
hlca_symbols <- hlca[["RNA"]]@meta.features$feature_name
names(hlca_symbols) <- rownames(hlca)
valid_genes <- !is.na(hlca_symbols) & hlca_symbols != ""
hlca_symbols_clean <- hlca_symbols[valid_genes]
hlca_symbols_clean <- hlca_symbols_clean[!duplicated(hlca_symbols_clean)]



## 1-2-2) Subset matrices ----
counts_mat <- hlca[["RNA"]]@counts[names(hlca_symbols_clean), ]
data_mat   <- hlca[["RNA"]]@data[names(hlca_symbols_clean), ]
meta_feat  <- hlca[["RNA"]]@meta.features[names(hlca_symbols_clean), ]



## 1-2-3) Update rownames to gene symbols ----
rownames(counts_mat) <- hlca_symbols_clean
rownames(data_mat)   <- hlca_symbols_clean
rownames(meta_feat)  <- hlca_symbols_clean



## 1-2-4) Build new Seurat object from scratch ----
hlca_clean <- CreateSeuratObject(
  counts = counts_mat,
  assay = "RNA",
  meta.data = hlca@meta.data,
  project = "HLCA"
)
# Confirm
head(rownames(hlca_clean), 10)
hlca_clean # 584944 cells



## 1-2-5) Select only needed tissue befroe Dimensional reduction ----
Idents(hlca_clean) <- "tissue"
table(Idents(hlca_clean))
hlca_clean <- subset(x = hlca_clean, 
                     idents = c("lung parenchyma"))
hlca_clean # 333468
Idents(hlca_clean) <- "cell_type"
table(Idents(hlca_clean))



## 1-2-6) Down sample ----
set.seed(42)

# Target total size
target_total <- 50000

# Get current counts per cell type
celltype_counts <- table(Idents(hlca_clean))

# Calculate sampling proportions
celltype_props <- celltype_counts / sum(celltype_counts)
celltype_sample_sizes <- round(celltype_props * target_total)

# Sample from each cell type
sampled_cells <- unlist(lapply(names(celltype_sample_sizes), function(ct) {
  cells <- WhichCells(hlca_clean, idents = ct)
  n <- min(length(cells), celltype_sample_sizes[[ct]])  # Don't oversample
  sample(cells, size = n)
}))

# Subset Seurat object
hlca_clean <- subset(hlca_clean, cells = sampled_cells)
hlca_clean # 49997 cells
Idents(hlca_clean) <- "cell_type"
table(Idents(hlca_clean))



## 1-2-7) Normalization to Dimensional reduction ----
hlca_clean <- NormalizeData(hlca_clean)            # Safe even if already normalized
hlca_clean <- FindVariableFeatures(hlca_clean)
hlca_clean <- ScaleData(hlca_clean)
hlca_clean <- RunPCA(hlca_clean, npcs = 30)
hlca_clean <- RunUMAP(hlca_clean, dims = 1:30)
# Checking things worked
Idents(hlca_clean) <- "cell_type"
table(Idents(hlca_clean))
DimPlot(hlca_clean, reduction = "umap", raster = TRUE, group.by = "cell_type", cols = "polychrome", 
        label = TRUE)
# Removing the HLCA (original object to save memory)
rm(hlca)
rm(counts_mat)
rm(data_mat)
rm(meta_feat)

# When using SCT
#hlca_clean <- SCTransform(hlca_clean, assay = "RNA", new.assay.name = "SCT", verbose = TRUE)
#hlca_clean
#DefaultAssay(hlca_clean) <- "SCT"
#hlca_clean <- NormalizeData(hlca_clean)
#hlca_clean <- FindVariableFeatures(hlca_clean)
#hlca_clean <- ScaleData(hlca_clean)
#hlca_clean <- RunPCA(hlca_clean, verbose = FALSE)
#hlca_clean <- RunUMAP(hlca_clean, dims = 1:30)

# We will use log normalize
DefaultAssay(hlca_clean) <- "RNA"
hlca_clean <- NormalizeData(hlca_clean)
hlca_clean <- FindVariableFeatures(hlca_clean)
hlca_clean <- ScaleData(hlca_clean)
hlca_clean <- RunPCA(hlca_clean, npcs = 30)

# Save
saveRDS(hlca_clean, file = "data/refernce_data_4_cellAnnotation/hlca_modifed_4_celltype_annotation.rds")





# 2. Annotation for nCount-RNA Harmnony ----
# 2-1) Import Data ---- 
hlca_clean <- readRDS("data/refernce_data_4_cellAnnotation/hlca_modifed_4_celltype_annotation.rds")
hlca_clean
czi_combined <- readRDS("output/round_2/output_4_step11_harmony/harmony_nCount_RNA/czi_combined_step10_clustering_and_UMAP_harmony_nCountRNA.rds")
czi_combined # 341046 nuclei




# 2-2) Process Data ----
## 2-2-1) Set assay defaults ----
DefaultAssay(czi_combined) <- "RNA"
DefaultAssay(hlca_clean) <- "RNA"



## 2-2-2) Normalize query ----
czi_combined <- NormalizeData(czi_combined)
czi_combined <- FindVariableFeatures(czi_combined)
czi_combined <- ScaleData(czi_combined)
czi_combined <- RunPCA(czi_combined, npcs = 30)



## 2-2-3) Find the common features between the reference and query ----
shared_features <- intersect(
  VariableFeatures(hlca_clean),
  VariableFeatures(czi_combined)
)
length(shared_features) # More the better, currently n = 511



## 2-2-4) Find TransferAcnchors ----
anchors <- FindTransferAnchors(
  reference = hlca_clean,
  query = czi_combined,
  reference.assay = "RNA",
  query.assay = "RNA",
  normalization.method = "LogNormalize",
  features = shared_features,
  dims = 1:30
)



## 2-2-5) Transfer Data ----
# Re-import the query to preserve the original PCA
czi_combined <- readRDS("output/round_2/output_4_step11_harmony/harmony_nCount_RNA/czi_combined_step10_clustering_and_UMAP_harmony_nCountRNA.rds")
czi_combined # 341046 nuclei
DefaultAssay(czi_combined) <- "SCT" # This line is not needed for Label Transfer
czi_combined


## 2-2-5-1) Ann_level_3 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = hlca_clean$ann_level_3,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "hlca_ann_level_3_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "hlca_ann_level_3_labelTransfer"
table(Idents(czi_combined))
# 13 x 8 L
DimPlot(czi_combined, reduction = "umap", 
        raster = TRUE, group.by = "hlca_ann_level_3_labelTransfer", cols = "polychrome", 
        label = TRUE)


## 2-2-5-2) Ann_level_1 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = hlca_clean$ann_level_1,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "hlca_ann_level_1_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "hlca_ann_level_1_labelTransfer"
table(Idents(czi_combined))
# 8 x 6 L
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_1_labelTransfer", 
              cols = "Pastel1", 
              label = TRUE)
p2 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_1_labelTransfer", 
              cols = "Accent", 
              label = TRUE)
# 16 x 6 l
p1 + p2


## 2-2-5-3) Ann_level_2 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = hlca_clean$ann_level_2,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "hlca_ann_level_2_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "hlca_ann_level_2_labelTransfer"
table(Idents(czi_combined))
# 8 x 6 L
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_2_labelTransfer", 
              cols = "Dark2", 
              label = TRUE)
p2 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_2_labelTransfer", 
              cols = "Pastel1", 
              label = TRUE)
# 16 x 6 l
p1 + p2


## 2-2-5-4) Ann_level_4 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = hlca_clean$ann_level_4,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "hlca_ann_level_4_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "hlca_ann_level_4_labelTransfer"
table(Idents(czi_combined))
# 8 x 6 L
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_4_labelTransfer", 
              cols = "polychrome", 
              label = TRUE)
p1


## 2-2-5-5) Ann_level_5 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = hlca_clean$ann_level_5,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "hlca_ann_level_5_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "hlca_ann_level_5_labelTransfer"
table(Idents(czi_combined))
# 8 x 6 L
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_5_labelTransfer", 
              cols = "polychrome", 
              label = TRUE)
p1


## 2-2-5-6) Finest level ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = hlca_clean$ann_finest_level,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "hlca_ann_level_finest_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "hlca_ann_level_finest_labelTransfer"
table(Idents(czi_combined))
# 8 x 6 L
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_finest_labelTransfer", 
              cols = "polychrome", 
              label = TRUE)
p1


## 2-2-5-7) Cell Type (used in the HLCA atlas & paper) ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = hlca_clean$cell_type,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "hlca_ann_level_cell_type_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "hlca_ann_level_cell_type_labelTransfer"
table(Idents(czi_combined))
# 8 x 6 L
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_cell_type_labelTransfer", 
              cols = "polychrome", 
              label = TRUE)
p1




# 2-3) Save Data Note the original PCA is preserved ----
saveRDS(czi_combined, file = "output/round_2/output_4_step12_cell_type_annotation/label_transfer/reference_hlca/harmony_nCountRNA/czi_combined_step12_harmony_nCountRNA.rds")
## 2-3-1) Save Tables ----
# Step 1: Inspect available metadata columns
head(colnames(czi_combined@meta.data))
# level 1
level1_table <- table(czi_combined$genetics, czi_combined$hlca_ann_level_1_labelTransfer)
level1_table
write.csv(level1_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/reference_hlca/harmony_nCountRNA/ann_level1.csv")
# level 2
level2_table <- table(czi_combined$genetics, czi_combined$hlca_ann_level_2_labelTransfer)
level2_table
write.csv(level2_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/reference_hlca/harmony_nCountRNA/ann_level2.csv")
# level 3
level3_table <- table(czi_combined$genetics, czi_combined$hlca_ann_level_3_labelTransfer)
level3_table
write.csv(level3_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/reference_hlca/harmony_nCountRNA/ann_level3.csv")
celltype_table <- table(czi_combined$genetics, czi_combined$hlca_ann_level_cell_type_labelTransfer)
celltype_table
write.csv(celltype_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/reference_hlca/harmony_nCountRNA/ann_celltype_table.csv")




# 2-4) Plot ----
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_1_labelTransfer", 
              cols = "Accent", 
              label = TRUE)
p2 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_2_labelTransfer", 
              cols = "Pastel2", 
              label = TRUE)
p3 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_3_labelTransfer",
              cols = "polychrome",
              label = TRUE)
# 8 x 32 l
p1 + p2 + p3
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_1_labelTransfer", 
              cols = "Accent", 
              label = FALSE)
p2 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_2_labelTransfer", 
              cols = "Pastel2", 
              label = FALSE)
p3 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "hlca_ann_level_3_labelTransfer",
              cols = "polychrome",
              label = FALSE)
# 8 x 32 l
p1 + p2 + p3





# 3. Utilizing Triston's data as reference ----
# 3-1) Import Data & Set up the Reference ----
triston <- readRDS("data/dataFrom_Tristan/FinalCellsProcessedandAnnotated.RDS")
triston # 97975 cells
head(rownames(triston), 10) # Gene symbols are stored as rows.
Idents(triston) <- "annotation_lvl1"
table(Idents(triston))
Idents(triston) <- "annotation_lvl2"
table(Idents(triston))
Idents(triston) <- "annotation_lvl3"
table(Idents(triston))
p1 <- DimPlot(triston, reduction = "umap", raster = TRUE, group.by = "annotation_lvl1", cols = "polychrome", 
              label = TRUE)
p2 <- DimPlot(triston, reduction = "umap", raster = TRUE, group.by = "annotation_lvl2", cols = "polychrome", 
              label = TRUE)
p3 <- DimPlot(triston, reduction = "umap", raster = TRUE, group.by = "annotation_lvl3", cols = "polychrome", 
              label = TRUE)
# 37 x 8 l
p1+p2+p3
# Variable features and PCA is not in this dataset 
triston <- FindVariableFeatures(triston)
triston <- ScaleData(triston)
triston <- RunPCA(triston, npcs = 30)
triston
# Save
saveRDS(triston, file = "data/dataFrom_Tristan/FinalCellsProcessedandAnnotated_withVariableFeature_and_pca.rds")




# 3-2) Annotation for nCount-RNA Harmnony ----
## 3-2-1) Import Data ---- 
triston <- readRDS("data/dataFrom_Tristan/FinalCellsProcessedandAnnotated_withVariableFeature_and_pca.rds")
triston
# We are importing the step10 object (just to be catious)
czi_combined <- readRDS("output/round_2/output_4_step11_harmony/harmony_nCount_RNA/czi_combined_step10_clustering_and_UMAP_harmony_nCountRNA.rds")
czi_combined # 341046 nuclei



## 3-2-2) Set assay defaults ----
DefaultAssay(czi_combined) <- "RNA"
DefaultAssay(triston) <- "RNA"



## 3-2-3) Normalize query ----
czi_combined <- NormalizeData(czi_combined)
czi_combined <- FindVariableFeatures(czi_combined)
czi_combined <- ScaleData(czi_combined)
czi_combined <- RunPCA(czi_combined, npcs = 30)



## 3-2-4) Find the common features between the reference and query ----
shared_features <- intersect(
  VariableFeatures(triston),
  VariableFeatures(czi_combined)
)
length(shared_features) # More the better, currently n = 1368 (This is more than HLCA)



## 3-2-5) Find TransferAcnchors ----
anchors <- FindTransferAnchors(
  reference = triston,
  query = czi_combined,
  reference.assay = "RNA",
  query.assay = "RNA",
  normalization.method = "LogNormalize",
  features = shared_features,
  dims = 1:30
)



## 3-2-6) Transfer Data ----
# Re-import the query to preserve the original PCA (we are importing that has the HLCA label transfer done. This object has the original PCA)
czi_combined <- readRDS("output/round_2/output_4_step12_cell_type_annotation/label_transfer/reference_hlca/harmony_nCountRNA/czi_combined_step12_harmony_nCountRNA.rds")
czi_combined # 341046 nuclei
DefaultAssay(czi_combined) <- "SCT" # This line is not needed for Label Transfer
czi_combined


## 3-2-6-1) annotation_lvl1 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = triston$annotation_lvl1,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "t_annotation_lvl1_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "t_annotation_lvl1_labelTransfer"
table(Idents(czi_combined))


## 3-2-6-2) annotation_lvl2 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = triston$annotation_lvl2,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "t_annotation_lvl2_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "t_annotation_lvl2_labelTransfer"
table(Idents(czi_combined))


## 3-2-6-3) annotation_lvl3 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = triston$annotation_lvl3,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "t_annotation_lvl3_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "t_annotation_lvl3_labelTransfer"
table(Idents(czi_combined))




# 3-3) Plot Data ----
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "t_annotation_lvl1_labelTransfer", 
              cols = "Accent", 
              label = TRUE)
p1
p2 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "t_annotation_lvl2_labelTransfer", 
              cols = "Set3", 
              label = TRUE)
p2
p3 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "t_annotation_lvl3_labelTransfer", 
              cols = "polychrome", 
              label = TRUE)
p3
# 37 x 8 l
p1+p2+p3




# 3-4) Save Tables ----
# Step 1: Inspect available metadata columns
head(colnames(czi_combined@meta.data))
# level 1
level1_table <- table(czi_combined$genetics, czi_combined$t_annotation_lvl1_labelTransfer)
level1_table
write.csv(level1_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/triston/harmony_ncountRNA/t_ann_level1.csv")
# level 2
level2_table <- table(czi_combined$genetics, czi_combined$t_annotation_lvl2_labelTransfer)
level2_table
write.csv(level2_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/triston/harmony_ncountRNA/t_ann_level2.csv")
# level 3
level3_table <- table(czi_combined$genetics, czi_combined$t_annotation_lvl3_labelTransfer)
level3_table
write.csv(level3_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/triston/harmony_nCountRNA/t_ann_level3.csv")




# 3-5) Save Data Note the original PCA is preserved ----
czi_combined
saveRDS(czi_combined, file = "output/round_2/output_4_step12_cell_type_annotation/label_transfer/triston/harmony_nCountRNA/czi_combined_step12_harmony_nCountRNA_with_hlca_tristion_labels.rds")


czi_combined
Idents(czi_combined) <- "hlca_ann_level_1_labelTransfer"
table(Idents(czi_combined))
immune_hlca <- subset(x = czi_combined,
                      ident = c("Immune"))
Idents(czi_combined) <- "t_annotation_lvl1_labelTransfer"
table(Idents(czi_combined))
immune_tris <- subset(x = czi_combined,
                      ident = c("Immune"))

Idents(immune_hlca) <- "t_annotation_lvl1_labelTransfer"
table(Idents(immune_hlca))
Idents(immune_tris) <- "hlca_ann_level_1_labelTransfer"
table(Idents(immune_tris))





# 4. Utilizing Elizabeth's (Jim) data as reference ----
# 4-1) Import Data & Set up the Reference ----
elizabeth <- readRDS("data/elizabeth_jim_data/sketch50K_Lung_10X-R_allPhase1_June25.rds")
elizabeth # 49567 nuclei
head(rownames(elizabeth), 10) # Gene symbols are stored as rows.
Idents(elizabeth) <- "class"
table(Idents(elizabeth))
Idents(elizabeth) <- "subclass.L1"
table(Idents(elizabeth))
Idents(elizabeth) <- "subclass.L2"
table(Idents(elizabeth))
Idents(elizabeth) <- "subclass.L3"
table(Idents(elizabeth))
Idents(elizabeth) <- "subclass.L4"
table(Idents(elizabeth))
Idents(elizabeth) <- "subclass.L5"
table(Idents(elizabeth))
p1 <- DimPlot(elizabeth, reduction = "umap", raster = TRUE, group.by = "class", cols = "polychrome", 
              label = TRUE)
p2 <- DimPlot(elizabeth, reduction = "umap", raster = TRUE, group.by = "subclass.L1", cols = "polychrome", 
              label = TRUE)
p3 <- DimPlot(elizabeth, reduction = "umap", raster = TRUE, group.by = "subclass.L2", cols = "polychrome", 
              label = TRUE)
# 30 x 8 l
p1+p2+p3

# Variable features and PCA is are in this dataset, and the following is not needed. 
# elizabeth <- FindVariableFeatures(elizabeth)
# elizabeth <- ScaleData(elizabeth)
# elizabeth <- RunPCA(elizabeth, npcs = 30)




# 4-2) Annotation for nCount-RNA Harmnony ----
## 4-2-1) Import Data ---- 
elizabeth
# We are importing the step10 object (just to be cautious)
czi_combined <- readRDS("output/round_2/output_4_step11_harmony/harmony_nCount_RNA/czi_combined_step10_clustering_and_UMAP_harmony_nCountRNA.rds")
czi_combined # 341046 nuclei



## 4-2-2) Set assay defaults ----
DefaultAssay(czi_combined) <- "RNA"
DefaultAssay(elizabeth) <- "RNA"



## 4-2-3) Normalize query ----
czi_combined <- NormalizeData(czi_combined)
czi_combined <- FindVariableFeatures(czi_combined)
czi_combined <- ScaleData(czi_combined)
czi_combined <- RunPCA(czi_combined, npcs = 30)



## 4-2-4) Find the common features between the reference and query ----
shared_features <- intersect(
  VariableFeatures(elizabeth),
  VariableFeatures(czi_combined)
)
length(shared_features) # More the better, currently n = 1717 (This is more than HLCA and Triston)



## 4-2-5) Find TransferAcnchors ----
anchors <- FindTransferAnchors(
  reference = elizabeth,
  query = czi_combined,
  reference.assay = "RNA",
  query.assay = "RNA",
  normalization.method = "LogNormalize",
  features = shared_features,
  dims = 1:30
)



## 4-2-6) Transfer Data ----
# Re-import the query to preserve the original PCA (we are importing that has the HLCA label transfer done. This object has the original PCA)
czi_combined <- readRDS("output/round_2/output_4_step12_cell_type_annotation/label_transfer/triston/harmony_nCountRNA/czi_combined_step12_harmony_nCountRNA_with_hlca_tristion_labels.rds")
czi_combined # 341046 nuclei
DefaultAssay(czi_combined) <- "SCT" # This line is not needed for Label Transfer
czi_combined


## 4-2-6-1) class ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = elizabeth$class,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "e_class_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "e_class_labelTransfer"
table(Idents(czi_combined))


## 4-2-6-2) subclass.L1 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = elizabeth$subclass.L1,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "e_subclassL1_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "e_subclassL1_labelTransfer"
table(Idents(czi_combined))


## 4-2-6-3) subclass.L2 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = elizabeth$subclass.L2,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "e_subclassL2_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "e_subclassL2_labelTransfer"
table(Idents(czi_combined))


## 4-2-6-4) subclass.L3 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = elizabeth$subclass.L3,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "e_subclassL3_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "e_subclassL3_labelTransfer"
table(Idents(czi_combined))


## 4-2-6-5) subclass.L4 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = elizabeth$subclass.L4,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "e_subclassL4_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "e_subclassL4_labelTransfer"
table(Idents(czi_combined))


## 4-2-6-6) subclass.L5 ----
predictions <- TransferData(
  anchorset = anchors,
  refdata = elizabeth$subclass.L5,
  dims = 1:30,
  slot = "data"  # Use the 'data' slot, which is the log-normalized values
)
# Confirm what predictions is
head(predictions)
str(predictions)
# Just take the predicted labels
label_transfer <- predictions[, 1, drop = FALSE]
# Rename that one column
colnames(label_transfer) <- "e_subclassL5_labelTransfer"
# Make sure the rownames match the Seurat object
rownames(label_transfer) <- rownames(predictions)
# Add it to the Seurat object
czi_combined <- AddMetaData(czi_combined, metadata = label_transfer)
czi_combined
Idents(czi_combined) <- "e_subclassL5_labelTransfer"
table(Idents(czi_combined))




# 4-3) Plot Data ----
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "e_class_labelTransfer", 
              cols = "Accent", 
              label = TRUE)
p1
p2 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "e_subclassL1_labelTransfer", 
              cols = "polychrome", 
              label = TRUE)
p2
p3 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "e_subclassL2_labelTransfer", 
              cols = "polychrome", 
              label = TRUE)
p3
# 37 x 8 l
p1+p2+p3




# 4-4) Save Tables ----
# Step 1: Inspect available metadata columns
head(colnames(czi_combined@meta.data))
# e_class_labelTransfer
level1_table <- table(czi_combined$genetics, czi_combined$e_class_labelTransfer)
level1_table
write.csv(level1_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/elizabeth/e_class_labelTransfer.csv")
# level 1
level2_table <- table(czi_combined$genetics, czi_combined$e_subclassL1_labelTransfer)
level2_table
write.csv(level2_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/elizabeth/e_subclassL1_labelTransfer.csv")
# level 2
level3_table <- table(czi_combined$genetics, czi_combined$e_subclassL2_labelTransfer)
level3_table
write.csv(level3_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/elizabeth/e_subclassL2_labelTransfer.csv")
# level 3
level4_table <- table(czi_combined$genetics, czi_combined$e_subclassL3_labelTransfer)
level4_table
write.csv(level4_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/elizabeth/e_subclassL3_labelTransfer.csv")
# level 4
level5_table <- table(czi_combined$genetics, czi_combined$e_subclassL4_labelTransfer)
level5_table
write.csv(level5_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/elizabeth/e_subclassL4_labelTransfer.csv")
# level 5
level6_table <- table(czi_combined$genetics, czi_combined$e_subclassL5_labelTransfer)
level6_table
write.csv(level5_table, "output/round_2/output_4_step12_cell_type_annotation/label_transfer/elizabeth/e_subclassL5_labelTransfer.csv")


# 3-5) Save Data Note the original PCA is preserved ----
czi_combined
saveRDS(czi_combined, file = "output/round_2/output_4_step12_cell_type_annotation/label_transfer/triston/harmony_nCountRNA/czi_combined_step12_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels.rds")
























