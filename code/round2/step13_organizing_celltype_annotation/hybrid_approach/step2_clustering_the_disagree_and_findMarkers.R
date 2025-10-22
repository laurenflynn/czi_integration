# Step2 Clustering the disagree and finding markers





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
library(tidyr)





# 1. Import Data ----
obj_dis <- readRDS("output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/obj_disagree_step13_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")
obj_dis # 17,152 nuclei




# 2. Process Data ----
DefaultAssay(obj_dis) <- "RNA" # just to be sure




# 2-1) Identify variable features ----
obj_dis <- FindVariableFeatures(obj_dis, selection.method = "vst", nfeatures = 3000)




# 2-2) Scale & PCA ----
obj_dis <- ScaleData(obj_dis, vars.to.regress = c("nCount_RNA"), verbose = FALSE) # To be consistent with the global object
obj_dis <- RunPCA(obj_dis, features = VariableFeatures(obj_dis))
ElbowPlot(obj_dis, ndims = 50)




# 2-3) Clustering ----
obj_dis <- FindNeighbors(obj_dis, dims = 1:30)   # adjust based on elbow
obj_dis <- FindClusters(obj_dis, resolution = 0.2)  # try 0.2–0.8
obj_dis <- RunUMAP(obj_dis, dims = 1:30)
# 4x6 l
DimPlot(obj_dis, group.by = "seurat_clusters", raster = TRUE,
        label = TRUE, repel = FALSE, cols = "polychrome") +
  ggtitle("Disagree subset: UMAP by new clusters")
DimPlot(obj_dis, group.by = "seurat_clusters", raster = TRUE,
        label = TRUE, repel = TRUE, cols = "polychrome", split.by = "dip_or_lympho") +
  ggtitle("Disagree subset: UMAP by new clusters")





# 3. To help manual annotation ----
# 3-1) Make module and get scores ----
epi    <- c("EPCAM","KRT18","MUC1","CLDN18","SFTPC")
endo   <- c("PECAM1","VWF","CDH5","CLDN5","KDR")
imm    <- c("PTPRC","LYZ","TYROBP","HLA-DRA","CD3D")
stroma <- c("COL1A1","COL1A2","DCN","PDGFRA","PDGFRB")
DefaultAssay(obj_dis) <- "RNA"
obj_dis <- AddModuleScore(obj_dis, list(epi),    name = "score_epi")
obj_dis <- AddModuleScore(obj_dis, list(endo),   name = "score_endo")
obj_dis <- AddModuleScore(obj_dis, list(imm),    name = "score_imm")
obj_dis <- AddModuleScore(obj_dis, list(stroma), name = "score_stroma")




# 3-2) Plot ----
## 3-2-1) Vln plot
# 6x8 l
VlnPlot(
  obj_dis,
  features = c("score_epi1","score_endo1","score_imm1","score_stroma1"),
  group.by = "seurat_clusters", pt.size = 0, ncol = 2
)



## 3-2-2) Heatmap of mean expression ----
scores <- FetchData(
  obj_dis,
  vars = c("score_epi1","score_endo1","score_imm1","score_stroma1","seurat_clusters")
)

# Disagree: cluster scoring → labels → heatmap
library(dplyr)
library(tibble)
library(pheatmap)
library(viridis)

# --- params ---
margin_thresh <- 0.08   # tune 0.05–0.10

# --- 1) Fetch per-cell scores once ---
scores <- FetchData(
  obj_dis,
  vars = c("score_epi1","score_endo1","score_imm1","score_stroma1","seurat_clusters")
)

# --- 2) Compute cluster means (one row per cluster) ---
cl_means <- scores %>%
  group_by(seurat_clusters) %>%
  summarize(
    epi    = mean(score_epi1,    na.rm = TRUE),
    endo   = mean(score_endo1,   na.rm = TRUE),
    imm    = mean(score_imm1,    na.rm = TRUE),
    stroma = mean(score_stroma1, na.rm = TRUE),
    .groups = "drop"
  )

# --- 3) Derive top class and margin (no list-cols, no repeats) ---
mat <- as.matrix(dplyr::select(cl_means, epi, endo, imm, stroma))
top_idx <- max.col(mat, ties.method = "first")
margin  <- apply(mat, 1, function(x){ s <- sort(x, decreasing = TRUE); s[1] - s[2] })
label_map <- c("Epithelial","Endothelial","Immune","Stroma")

cl_means <- cl_means %>%
  mutate(
    top_idx     = top_idx,
    top_label   = label_map[top_idx],
    margin      = margin,
    final_label = ifelse(margin >= margin_thresh, top_label, "not_sure")
  ) %>%
  arrange(seurat_clusters)

# quick peek
print(cl_means, n = Inf)

# --- 4) Heatmap (row-scaled) of cluster means ---
mat_means <- cl_means %>%
  dplyr::select(seurat_clusters,
                Epithelial = epi,
                Endothelial = endo,
                Immune = imm,
                Stroma = stroma) %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  column_to_rownames("seurat_clusters") %>%
  as.matrix()

stopifnot(is.matrix(mat_means), nrow(mat_means) > 0)

# 3x8 p
pheatmap(
  mat_means,
  color = viridis(256),
  cluster_rows = TRUE, cluster_cols = FALSE,
  scale = "row",
  main = sprintf("Module score means per cluster (row-scaled), margin ≥ %.2f", margin_thresh)
)





# 4. Save ----
saveRDS(obj_dis, file = "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step2_clustering_the_disagree_and_findMarkers/obj_disagree_step13_2_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")
obj_dis <- readRDS("output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step2_clustering_the_disagree_and_findMarkers/obj_disagree_step13_2_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")
obj_dis # 17,152 nuclei




# 5. Find Markers
markers_dis <- FindAllMarkers(obj_dis, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_dis, 
          file = "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step2_clustering_the_disagree_and_findMarkers/de_markers_4_clusters.csv", 
          row.names = TRUE)



