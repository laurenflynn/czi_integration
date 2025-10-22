# Step1 Find the commonl in level1

# Note: We will be using Seurat v5 'Reference Mapping'.
# The reference will be Tristan & Elizabeth
# We will be identifying the population that agrees using the two datasets.





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
czi_combined <- readRDS("output/round_2/output_4_step12_cell_type_annotation/label_transfer/triston/harmony_nCountRNA/czi_combined_step12_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")
czi_combined # 341046 nuclei
table(Idents(czi_combined))
Idents(czi_combined) <- "updated_sample_id"
table(Idents(czi_combined))





# 2. Process Data ----
# 2-1) Organize reference data (Elizabeth) ----
# Elizabeth data in level one the mesenchyme cell type is divided into two cell types
# They need to be combined
Idents(czi_combined) <- "e_class_labelTransfer"
table(Idents(czi_combined))

# Create a new column in the metadata to store the cell type annotations
czi_combined$e_annotation_lvl1 <- Idents(czi_combined)
czi_combined$e_annotation_lvl1 <- plyr::mapvalues(czi_combined$e_annotation_lvl1, 
                                        from = c("Epithelial", "Endothelial", "Connective", "Immune", "PNS"),
                                        to = c("Epithelial", "Endothelial", "Mesenchyme", "Immune", "Mesenchyme"))
Idents(czi_combined) <- "e_annotation_lvl1"
table(Idents(czi_combined))




# 2-2) Save ----
saveRDS(czi_combined, file = "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/czi_combined_step13_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")





# 3. Identify agreement & disagreement ----
# 3-1) Setup harmonized labels (making the names the same) ---- 
library(purrr)
library(stringr)

# our object with both label-transfer results already stored:
# - t_annotation_lvl1_labelTransfer  (Tristan)
# - e_annotation_lvl1                (Elizabeth)

lvl1_map <- c(
  "Epithelium"   = "Epithelial",
  "Epithelial"   = "Epithelial",
  "Endothelium"  = "Endothelial",
  "Endothelial"  = "Endothelial",
  "Mesenchyme"   = "Stroma",
  "Stroma"       = "Stroma",
  "Immune"       = "Immune"
)

czi_combined$L1_tristan   <- unname(lvl1_map[ as.character(czi_combined$t_annotation_lvl1_labelTransfer) ])
czi_combined$L1_elizabeth <- unname(lvl1_map[ as.character(czi_combined$e_annotation_lvl1) ])

# sanity
table(Tristan   = czi_combined$L1_tristan,   useNA = "ifany")
table(Elizabeth = czi_combined$L1_elizabeth, useNA = "ifany")
sum(is.na(czi_combined$L1_elizabeth))



## 3-1-1) Elizabeth as baseline ----
# Cross-tab (rows = Elizabeth, cols = Tristan)
tab_L1 <- table(Elizabeth = czi_combined$L1_elizabeth,
                Tristan   = czi_combined$L1_tristan)
# Row-normalize (%)
prop_eliz <- round(prop.table(tab_L1, 1) * 100, 1)  # row-normalized (% by Elizabeth call)
prop_eliz
# Convert to data frame for ggplot
df_heat <- as.data.frame(as.table(round(prop_eliz, 1)))
head(df_heat)



## 3-1-2) Tristan as baseline ----
tab_L1 <- table(Tristan   = czi_combined$L1_tristan,
                Elizabeth = czi_combined$L1_elizabeth)
# Normalize by columns = Tristan
prop_tris <-round(prop.table(tab_L1, 1) * 100, 1)
df_heat2 <- as.data.frame(as.table(round(prop_tris, 1)))



## 3-1-3) Plot ---
p1 <- ggplot(df_heat, aes(x = Tristan, y = Elizabeth, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Freq),
                color = Freq > 70), size = 4) +
  scale_color_manual(values = c("TRUE"="black","FALSE"="white"), guide = "none") +
  scale_fill_viridis_c(option = "inferno", limits = c(0,100)) +
  labs(title = "Elizabeth → Tristan",
       x = "Tristan annotation", y = "Elizabeth annotation", fill = "%") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(df_heat2, aes(x = Elizabeth, y = Tristan, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f", Freq),
                color = Freq > 70), size = 4) +
  scale_color_manual(values = c("TRUE"="black","FALSE"="white"), guide = "none") +
  scale_fill_viridis_c(option = "inferno", limits = c(0,100)) +
  labs(title = "Tristan → Elizabeth",
       x = "Elizabeth annotation", y = "Tristan annotation", fill = "%") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 4 x 9.2 l
p1 + p2  # side by side

# Overall exact-match rate
agree_mask <- czi_combined$L1_elizabeth == czi_combined$L1_tristan
mean(agree_mask)  # e.g., 0.82 = 82% agreement

# Per-class agreement when Elizabeth = "Epithelial"
epi_mask <- czi_combined$L1_elizabeth == "Epithelial"
mean(czi_combined$L1_tristan[epi_mask] == "Epithelial")
# result shows 0.9497077




# 3-2) Lock in Level-1 by agreement ----
czi_combined$L1_status <- ifelse(agree_mask, "Agree", "Disagree")
czi_combined$L1_final  <- ifelse(agree_mask, czi_combined$L1_elizabeth, NA)
table(czi_combined$L1_final, useNA = "ifany")




# 3-3) Plot the UMAP of Agree & Disagree ----
cellColors = c("Epithelial" = "#c0b1ce",  "Endothelial" = "#88c97f", 
               "Immune" = "#eabe9e", "Stroma" = "#f4ed9e")
p1 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "L1_tristan",
              split.by = "L1_status",
              cols = cellColors, 
              label = TRUE)
# 5x10l
p1
p2 <- DimPlot(czi_combined, 
              reduction = "umap", 
              raster = TRUE, 
              group.by = "L1_elizabeth",
              split.by = "L1_status",
              cols = cellColors, 
              label = TRUE)
# 5x10l
p2
# 5 x 20 l
p1 | p2





# 4. Subset out the agree & plot markers ----
czi_combined # 341046 nuclei
table(Idents(czi_combined))
Idents(czi_combined) <- "L1_status"
table(Idents(czi_combined))
czi_combined




# 4-1) Subset the agreed ----
obj_agree <- subset(x = czi_combined, idents = c("Agree"))
obj_agree
Idents(obj_agree) <- "L1_elizabeth"
table(Idents(obj_agree))
obj_agree
Idents(obj_agree) <- "L1_tristan"
table(Idents(obj_agree))
obj_agree
Idents(obj_agree) <- "L1_final"
table(Idents(obj_agree))
# Set up the order
obj_agree$L1_final <- factor(obj_agree$L1_final, 
                             levels = c("Epithelial","Stroma","Endothelial", "Immune"))
obj_agree_table <- table(obj_agree$L1_final, obj_agree$genetics)
obj_agree_table
write.csv(obj_agree_table, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/l1_agree_vs_genetics.csv")
total_table <- table(Idents(obj_agree))
total_table
write.csv(obj_agree_table, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/l1_agree_vs_genetics_total.csv")




# 4-2) Plot marker expression ----
library(dplyr)
library(tidyr)
Idents(obj_agree) <- "L1_final"
table(Idents(obj_agree))
epi    <- c("EPCAM","KRT8","KRT18","KRT19","MUC1")
endo   <- c("PECAM1","CLDN5","KDR","VWF","CDH5")
immune <- c("PTPRC","LYZ","TYROBP","HLA-DRA","CD3D")
stroma <- c("COL1A1","COL1A2","DCN","PDGFRA","PDGFRB")
genes  <- c(epi, endo, immune, stroma)

dp <- DotPlot(obj_agree, features = genes, group.by = "L1_final")$data
# dp columns include: features.plot (gene), id (group), pct.exp, avg.exp.scaled, avg.exp
# Wide table: rows = disease, cols = genes, values = % expressing
pct_table <- dp %>%
  dplyr::select(disease = id, gene = features.plot, pct.exp) %>%
  tidyr::pivot_wider(names_from = gene, values_from = pct.exp)
pct_table
write.csv(pct_table, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/obj_agree_4_dotplot_percentage.csv")

# Re-sort
epi    <- c("MUC1","KRT19","EPCAM","KRT18","KRT8")
endo   <- c("PECAM1","VWF","CLDN5","CDH5","KDR")
immune <- c("PTPRC","HLA-DRA","TYROBP","LYZ","CD3D")
stroma <- c("COL1A2","COL1A1","PDGFRA","DCN","PDGFRB")
genes  <- c(epi, stroma, endo, immune)



## 4-2-1) Agree ----
p_agree <- DotPlot(obj_agree, features = genes, cols = c("grey90","firebrick"), dot.scale = 4) +
  scale_colour_gradient2(low = "#000004", mid = "#bc3754", high = "#fde725", midpoint = 0) +
  RotatedAxis() + ggtitle("Agreement L1")
p_agree
stackVln <- VlnPlot(obj_agree,
                    features = genes,
                    cols = cellColors,
                    log = FALSE,
                    assay = "RNA",
                    stack = TRUE, 
                    fill.by = "ident") + NoLegend()
# l 4x13
p_agree + stackVln 



## 4-2-2) Tristan ----
czi_combined$L1_tristan <- factor(czi_combined$L1_tristan, 
                             levels = c("Epithelial","Stroma","Endothelial", "Immune"))
Idents(czi_combined) <- "L1_tristan"
table(Idents(czi_combined))
p_tris <- DotPlot(czi_combined, features = genes, cols = c("grey90","firebrick"), dot.scale = 4) +
  scale_colour_gradient2(low = "#000004", mid = "#bc3754", high = "#fde725", midpoint = 0) +
  RotatedAxis() + ggtitle("L1_tristan")
p_tris
stackVln <- VlnPlot(czi_combined,
                    features = genes,
                    cols = cellColors,
                    log = FALSE,
                    assay = "RNA",
                    stack = TRUE, 
                    fill.by = "ident") + NoLegend()
# l 4x13
p_tris + stackVln 



## 4-2-3) Elizabeth ----
czi_combined$L1_elizabeth <- factor(czi_combined$L1_elizabeth, 
                                  levels = c("Epithelial","Stroma","Endothelial", "Immune"))
Idents(czi_combined) <- "L1_elizabeth"
table(Idents(czi_combined))
p_eliz <- DotPlot(czi_combined, features = genes, cols = c("grey90","firebrick"), dot.scale = 4) +
  scale_colour_gradient2(low = "#000004", mid = "#bc3754", high = "#fde725", midpoint = 0) +
  RotatedAxis() + ggtitle("L1_elizabeth")
p_eliz
stackVln <- VlnPlot(czi_combined,
                    features = genes,
                    cols = cellColors,
                    log = FALSE,
                    assay = "RNA",
                    stack = TRUE, 
                    fill.by = "ident") + NoLegend()
# l 4x13
p_eliz + stackVln 




# 4-3) Save agreed object ----
saveRDS(obj_agree, file = "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/l1_temp_agreed_step13_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")
rm(obj_agree) # To save memory





# 5. Subset the disagreed ----
czi_combined # 341046 nuclei
table(Idents(czi_combined))
Idents(czi_combined) <- "L1_status"
table(Idents(czi_combined))
czi_combined




# 5-1) Subset the disagreed ----
obj_dis <- subset(x = czi_combined, idents = c("Disagree"))
obj_dis # 17152 cells



## 5-1-1) Tristan ----
Idents(obj_dis) <- "L1_tristan"
table(Idents(obj_dis))
p_tris <- DotPlot(obj_dis, features = genes, cols = c("grey90","firebrick"), dot.scale = 4) +
  scale_colour_gradient2(low = "#000004", mid = "#bc3754", high = "#fde725", midpoint = 0) +
  RotatedAxis() + ggtitle("L1_tristan_disagree")
p_tris
stackVln <- VlnPlot(obj_dis,
                    features = genes,
                    cols = cellColors,
                    log = FALSE,
                    assay = "RNA",
                    stack = TRUE, 
                    fill.by = "ident") + NoLegend()
# l 4x13
p_tris + stackVln 



## 5-1-2) Elizabeth ----
Idents(obj_dis) <- "L1_elizabeth"
table(Idents(obj_dis))
p_eliz <- DotPlot(obj_dis, features = genes, cols = c("grey90","firebrick"), dot.scale = 4) +
  scale_colour_gradient2(low = "#000004", mid = "#bc3754", high = "#fde725", midpoint = 0) +
  RotatedAxis() + ggtitle("L1_elizabeth_disagree")
p_eliz
stackVln <- VlnPlot(obj_dis,
                    features = genes,
                    cols = cellColors,
                    log = FALSE,
                    assay = "RNA",
                    stack = TRUE, 
                    fill.by = "ident") + NoLegend()
# l 4x13
p_eliz + stackVln 





# 6. Investigate if the the disagree cells are enriched in 'genetic' conditions ----
Idents(obj_dis) <- "genetics"
table(Idents(obj_dis))
Idents(obj_dis) <- "updated_sample_id"
table(Idents(obj_dis))
obj_disagree_table <- table(obj_dis$genetics, obj_dis$updated_sample_id)
obj_disagree_table
write.csv(obj_disagree_table, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/l1_disagree_sample_vs_genetics.csv")
Idents(obj_dis) <- "genetics"
table(Idents(obj_dis))
total_table <- table(Idents(obj_dis))
total_table
write.csv(total_table, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/l1_disagree_sample_vs_genetics_total.csv")
# Save Data
saveRDS(obj_dis, file = "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/obj_disagree_step13_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")




# 6-1) Are some genetics enriched for disagreement? ----
library(dplyr)
# Rows: Status = Agree or Disagree
# Columns: Genetics = one coulumn per genotype
tab <- table(Status = czi_combined$L1_status, Genetics = czi_combined$genetics)
fisher_df <- lapply(colnames(tab), function(g){
  # 2×2 matrix: rows = Disagree/Agree, cols = g vs not g
  # Build a 2x2 for "genetics = g" vs "all other genetics"
  m <- matrix(c(tab["Disagree", g], tab["Agree", g],
                sum(tab["Disagree", ]) - tab["Disagree", g],
                sum(tab["Agree", ]) - tab["Agree", g]),
              nrow=2, byrow=TRUE)
  ft <- fisher.test(m)
  data.frame(
    genetics = g,
    OR  = as.numeric(ft$estimate),
    lo  = ft$conf.int[1],
    hi  = ft$conf.int[2],
    p   = ft$p.value,
    disagree_in_genetics = tab["Disagree", g],
    total_in_genetics    = sum(tab[, g])
  )
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p, method="fdr"),
         rate_disagree = disagree_in_genetics / total_in_genetics) %>%
  arrange(desc(OR)) %>%
  tibble::as_tibble()   # <— force tibble class
# Now dplyr::select() will work fine
fisher_df %>%
  dplyr::select(genetics, OR, lo, hi, p, FDR, rate_disagree) %>%
  print(n = Inf)




# 6-2) Plot ----
# Forest plot
library(ggplot2)
ggplot(fisher_df, aes(x = reorder(genetics, OR), y = OR)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_pointrange(aes(ymin = lo, ymax = hi)) +
  coord_flip() +
  scale_y_log10() +
  labs(y = "Odds ratio (log scale)", x = "Genetics",
       title = "Enrichment of disagreement by genetics (Fisher exact)") +
  theme_minimal(base_size = 12)




# 6-3) Investigate if the disagree compared to agree have different quality ----
meta <- czi_combined@meta.data
p1 <- ggplot(meta, aes(L1_status, nCount_RNA)) + geom_violin() + geom_boxplot(width=0.1)
p2 <- ggplot(meta, aes(L1_status, nFeature_RNA)) + geom_violin() + geom_boxplot(width=0.1)
p3 <- ggplot(meta, aes(L1_status, percent.mt)) + geom_violin() + geom_boxplot(width=0.1)
# 3x 15 l
p1 | p2 | p3
pvals <- c(
  wilcox.test(nCount_RNA ~ L1_status, data = meta)$p.value,
  wilcox.test(nFeature_RNA ~ L1_status, data = meta)$p.value,
  wilcox.test(percent.mt ~ L1_status, data = meta)$p.value
)
p.adjust(pvals, method = "BH")

# Cliff’s Delta (non-parametric effect size)
library(effsize)
cliff.delta(nCount_RNA ~ L1_status, data = meta)
cliff.delta(nFeature_RNA ~ L1_status, data = meta)
cliff.delta(percent.mt ~ L1_status, data = meta)
# Ranges: -1 to 1.
# Rules of thumb: |δ| < 0.147 → negligible |δ| < 0.33 → small |δ| < 0.474 → medium otherwise → large

# log normalized based
meta$log_nCount_RNA   <- log1p(meta$nCount_RNA)
meta$log_nFeature_RNA <- log1p(meta$nFeature_RNA)

qc_vars <- c("nCount_RNA","nFeature_RNA","log_nCount_RNA","log_nFeature_RNA","percent.mt")

qc_results <- lapply(qc_vars, function(v) {
  form <- as.formula(paste(v, "~ L1_status"))
  wt <- wilcox.test(form, data = meta)
  cd <- effsize::cliff.delta(form, data = meta)
  
  data.frame(
    metric = v,
    p_value = wt$p.value,
    cliff_delta = cd$estimate
  )
}) %>% bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

qc_results



# To find where the disagreement is the most
dsub <- czi_combined@meta.data |> 
  filter(L1_status=="Disagree") |>
  transmute(Elizabeth=L1_elizabeth, Tristan=L1_tristan)

swap_tab <- table(Elizabeth=dsub$Elizabeth, Tristan=dsub$Tristan)
swap_tab
# Heatmap:
library(pheatmap)
pheatmap(swap_tab, cluster_rows=FALSE, cluster_cols=FALSE,
         color=colorRampPalette(c("white","firebrick3"))(101),
         main="Disagreement pairs (counts)")

library(viridis)
# 4x4.3 l
pheatmap(
  swap_tab,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = plasma(101),  # instead of white → firebrick3
  main = "Disagreement pairs (counts)"
)

df <- as.data.frame(as.table(swap_tab))
colnames(df) <- c("Elizabeth","Tristan","value")

# Find the max (the yellow square in plasma)
max_val <- max(df$value, na.rm = TRUE)
# 4x4.3 l
ggplot(df, aes(Tristan, Elizabeth, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value,
                color = (value == max_val))) +  # TRUE only for the yellow square
  scale_fill_viridis(option = "plasma", direction = 1) +  # yellow = high
  scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
  coord_fixed() +
  labs(title = "Disagreement pairs (counts)", x = "Tristan", y = "Elizabeth") +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank())

# Aditional check
DimPlot(obj_dis, reduction = "umap", group.by = "genetics", raster = "TRUE", cols = "polychrome")
DimPlot(obj_dis, reduction = "umap", group.by = "genetics", split.by = "genetics", raster = "TRUE", cols = "polychrome", ncol = 4) +
  ggtitle("Disagreement UMAP split by genetics")

meta <- as.data.frame(meta)
tab1 <- meta |>
  dplyr::count(L1_status) |>
  dplyr::mutate(pct = 100 * n / sum(n))
tab1

library(dplyr)
library(tidyr)
# Start from metadata
meta <- as.data.frame(czi_combined@meta.data)
# 1) Build a consensus Level-1 label (no bias)
meta$consensus_lvl1 <- ifelse(meta$L1_status == "Agree" & !is.na(meta$L1_final),
                              as.character(meta$L1_final),
                              "Uncertain")
# (Optional) order the levels for nicer display
lvl_order <- c("Epithelial","Stroma","Endothelial","Immune","Uncertain")
meta$consensus_lvl1 <- factor(meta$consensus_lvl1, levels = lvl_order)
# 2) Your tab2: counts and row-wise % by consensus level
meta <- as.data.frame(meta)
tab2_eliz <- as.data.frame(czi_combined@meta.data) %>%
  mutate(L1_elizabeth = as.character(L1_elizabeth),
         L1_status    = as.character(L1_status)) %>%
  dplyr::count(L1_elizabeth, L1_status, .drop = FALSE) %>%
  dplyr::group_by(L1_elizabeth) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  tidyr::pivot_wider(names_from = L1_status, values_from = c(n, pct), values_fill = 0) %>%
  dplyr::arrange(dplyr::desc(pct_Disagree))
tab2_eliz
tab2_tris <- as.data.frame(czi_combined@meta.data) %>%
  mutate(L1_tristan = as.character(L1_tristan),
         L1_status  = as.character(L1_status)) %>%
  dplyr::count(L1_tristan, L1_status, .drop = FALSE) %>%
  dplyr::group_by(L1_tristan) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  tidyr::pivot_wider(names_from = L1_status, values_from = c(n, pct), values_fill = 0) %>%
  dplyr::arrange(dplyr::desc(pct_Disagree))
tab2_tris
tab_agree <- as.data.frame(czi_combined@meta.data) %>%
  dplyr::filter(L1_status == "Agree") %>%
  dplyr::count(L1_final) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::arrange(desc(n))
tab_agree
tab_genetics <- as.data.frame(czi_combined@meta.data) %>%
  dplyr::count(genetics, L1_status) %>%
  dplyr::group_by(genetics) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  tidyr::pivot_wider(names_from = L1_status, values_from = c(n, pct), values_fill = 0) %>%
  dplyr::arrange(desc(pct_Disagree))
tab_genetics
write.csv(tab_genetics, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/tab_genetics.csv")
tab_sample <- as.data.frame(czi_combined@meta.data) %>%
  dplyr::count(updated_sample_id, L1_status) %>%
  dplyr::group_by(updated_sample_id) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  tidyr::pivot_wider(names_from = L1_status, values_from = c(n, pct), values_fill = 0)
tab_sample
write.csv(tab_sample, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/tab_samples.csv")




# 6-4) Investigate if there is a difference in age ----
library(dplyr)
library(ggplot2)

obj_dis <- readRDS("output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/obj_disagree_step13_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")
obj_dis # 17152 nuclei
# Disagreement-only metadata (one row per sample to avoid overcounting cells)
meta_dis <- obj_dis@meta.data %>%
  dplyr::select(updated_sample_id, genetics, age) %>%
  distinct()

# Define enriched vs stable groups from your forest plot
enriched <- c("ABCA3","SFTPC","FLNA","STAT1","Control_p")
stable   <- c("SFTPB","GRN","NLRP12","FARS2","NPC2", "SLC7A7", "Unknown", "PIK3CA", 
              "IKBKB", "DICER", "LRBA", "Control", "SOCS1")

meta_dis <- meta_dis %>%
  mutate(disagree_group = case_when(
    genetics %in% enriched ~ "Enriched",
    genetics %in% stable   ~ "Stable",
    TRUE                   ~ "Other"
  ))

# ---- Test only Enriched vs Stable ----
df_test <- meta_dis %>% filter(disagree_group %in% c("Enriched","Stable"))

# Wilcoxon test (non-parametric, safer with small n)
wilcox.test(age ~ disagree_group, data = df_test)

# Or a t-test if age looks roughly normal
t.test(age ~ disagree_group, data = df_test)

# ---- Plot ----
# 3x7 l
ggplot(df_test, aes(x = disagree_group, y = age, fill = disagree_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  labs(x = "", y = "Age (days)",
       title = "Age of samples contributing to disagreement cells",
       subtitle = "Comparing enriched vs stable genetics") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")



## 6-4-1) This is on global. Disagree vs age ----
meta_sample <- as.data.frame(czi_combined@meta.data) %>%
  group_by(updated_sample_id, genetics, age) %>%
  summarise(
    total = n(),
    disagree = sum(L1_status == "Disagree"),
    rate_disagree = disagree / total,
    .groups = "drop"
  )
# Correlation
cor.test(meta_sample$rate_disagree, meta_sample$age, method = "spearman")
genetics_palette <- c(
  "SFTPB"     = "#f38187",
  "Control"   = "#b6dbab",
  "Control_p" = "#009444",
  "ABCA3"     = "#80a6cc",
  "FARS2"     = "#8cc6ec",
  "GRN"       = "#e7eb94",
  "LRBA"      = "#c75186",
  "NFKB_NOD2" = "#fff450",
  "NLRP12"    = "#f9cbdf",
  "SFTPC"     = "#d9dad9",
  "SLC7A7"    = "#59c1ad",
  "SOCS1"     = "#fec369",
  "STAT1"     = "#782b90",
  "DICER"     = "#64735a",
  "FLNA"      = "#faf8cd",
  "IKBKB"     = "#c2b59b",
  "NPC2"      = "#ec008c",
  "PIK3CA"    = "#cb1d2e",
  "Unknown"   = "#262262"
)
# 5x6.3 l
ggplot(meta_sample, aes(x = age, y = rate_disagree, color = genetics)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = genetics_palette, drop = FALSE) +
  labs(
    x = "Age",
    y = "Disagreement rate (%)",
    title = "Per-sample disagreement vs donor age",
    subtitle = "Each point = one sample; dashed line = linear fit",
    color = "Genetics"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
# For QC metric
meta_sample_qc <- as.data.frame(czi_combined@meta.data) %>%
  group_by(updated_sample_id, genetics, age) %>%
  summarise(
    total       = n(),
    disagree    = sum(L1_status == "Disagree"),
    rate_disagree = disagree / total,
    median_nCount   = median(nCount_RNA),
    median_nFeature = median(nFeature_RNA),
    median_mt       = median(percent.mt),
    .groups = "drop"
  )
cor.test(meta_sample_qc$rate_disagree, meta_sample_qc$median_nCount,   method = "spearman")
cor.test(meta_sample_qc$rate_disagree, meta_sample_qc$median_nFeature, method = "spearman")
cor.test(meta_sample_qc$rate_disagree, meta_sample_qc$median_mt,       method = "spearman")
# Example for nCount
p1 <- ggplot(meta_sample_qc, aes(x = median_nCount, y = rate_disagree, color = genetics)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(values = genetics_palette, drop = FALSE) +
  labs(
    x = "Median nCount_RNA per sample",
    y = "Disagreement rate (%)",
    title = "Per-sample disagreement vs sequencing depth",
    subtitle = "Each point = one sample; dashed line = linear fit",
    color = "Genetics"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")





# 7. Are some samples enrichced for disagreement? ---- 
# 7-1) Re-do 'Are some genetics enriched for disagreement?' ----
library(dplyr)
fmtp  <- function(x) ifelse(x < 1e-16, "<1e-16", formatC(x, format="e", digits=2))
stars <- function(x) case_when(x < 0.001 ~ "***", x < 0.01 ~ "**", x < 0.05 ~ "*", TRUE ~ "")

# our existing fisher_df code (unchanged)
tab <- table(Status = czi_combined$L1_status, Genetics = czi_combined$genetics)
fisher_df <- lapply(colnames(tab), function(g){
  m <- matrix(c(tab["Disagree", g], tab["Agree", g],
                sum(tab["Disagree", ]) - tab["Disagree", g],
                sum(tab["Agree", ])    - tab["Agree", g]),
              nrow=2, byrow=TRUE)
  ft <- fisher.test(m)
  data.frame(
    genetics = g,
    OR  = as.numeric(ft$estimate),
    lo  = ft$conf.int[1],
    hi  = ft$conf.int[2],
    p   = ft$p.value,
    disagree_in_genetics = tab["Disagree", g],
    total_in_genetics    = sum(tab[, g])
  )
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p, "fdr"),
         rate_disagree = disagree_in_genetics / total_in_genetics) %>%
  tibble::as_tibble()

# Add log2-scale columns + labels
fisher_df_plot <- fisher_df %>%
  mutate(
    log2OR = log2(OR), log2lo = log2(lo), log2hi = log2(hi),
    p_lab = fmtp(p), FDR_lab = fmtp(FDR), star = stars(FDR)
  )
library(ggplot2)
ggplot(fisher_df_plot, aes(x = reorder(genetics, log2OR), y = log2OR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = log2lo, ymax = log2hi), width = 0.2) +
  geom_text(aes(label = star), nudge_y = 0.10, size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(y = "log2(OR)  (0 = no enrichment)", x = "Genetics",
       title = "Genetics-level disagreement enrichment",
       subtitle = "95% CI on log2 scale; stars = FDR") +
  theme_minimal(base_size = 13)




# 7-2) Per sample
# Build per-sample table (your code)
tab_sample <- table(Status = czi_combined$L1_status,
                    Sample = czi_combined$updated_sample_id)

# (Optional) continuity-corrected Fisher (guards zeros)
safe_fisher <- function(a,b,c,d){
  if(min(a,b,c,d) == 0){
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  }
  fisher.test(matrix(c(a,b,c,d), nrow=2, byrow=TRUE))
}

fisher_sample <- lapply(colnames(tab_sample), function(s){
  a <- tab_sample["Disagree", s]
  b <- tab_sample["Agree",    s]
  c <- sum(tab_sample["Disagree", ]) - a
  d <- sum(tab_sample["Agree",    ]) - b
  ft <- safe_fisher(a,b,c,d)
  data.frame(
    sample   = s,
    OR       = as.numeric(ft$estimate),
    lo       = ft$conf.int[1],
    hi       = ft$conf.int[2],
    p        = ft$p.value,
    disagree_in_sample = a,
    total_in_sample    = a + b
  )
}) %>% dplyr::bind_rows() %>%
  dplyr::mutate(FDR = p.adjust(p, "fdr"))

# Join genetics per sample
sample_meta <- czi_combined@meta.data %>%
  dplyr::select(updated_sample_id, genetics) %>% distinct()

fisher_sample <- fisher_sample %>%
  dplyr::left_join(sample_meta, by = c("sample" = "updated_sample_id"))

# Colors
library(ggsci)
genetics_levels <- unique(fisher_sample$genetics)
genetics_colors <- setNames(ggsci::pal_d3("category20")(length(genetics_levels)),
                            genetics_levels)

# Add log2 columns + labels
fisher_sample_plot <- fisher_sample %>%
  mutate(
    log2OR = log2(OR), log2lo = log2(lo), log2hi = log2(hi),
    p_lab = fmtp(p), FDR_lab = fmtp(FDR), star = stars(FDR)
  )

ggplot(fisher_sample_plot, aes(x = reorder(sample, log2OR), y = log2OR, color = genetics)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = log2lo, ymax = log2hi), width = 0.2) +
  geom_text(aes(label = star), nudge_y = 0.10, size = 4, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_color_manual(values = genetics_colors) +
  labs(y = "log2(OR)  (0 = no enrichment)", x = "Sample",
       title = "Sample-level disagreement enrichment",
       subtitle = "95% CI on log2 scale; stars = FDR",
       color = "Genetics") +
  theme_minimal(base_size = 13)




# 7-3) Save ----
saveRDS(czi_combined, file = "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/czi_combined_step13_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")





# 8. Plot to match colors ----
# 8-1) Import Data ----
czi_combined <- readRDS("output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/czi_combined_step13_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")
czi_combined # 341046 nuclei
Idents(czi_combined) <- "genetics"
table(Idents(czi_combined))
# Color map:
# SFTPB: #f38187, Control: #b6dbab, Control_p: #009444, ABCA3: #80a6cc, FARS2: #8cc6ec, GRN: #e7eb94,
# LRBA: #c75186, NFKB_NOD2: #fff450, NLRP12: #f9cbdf, SFTPC: #d9dad9, SLC7A7: #59c1ad, SOCS1: #fec369
# STAT1: #782b90, DICER: #64735a, FLNA: #faf8cd, IKBKB: #c2b59b,NPC2: #ec008c, PIK3CA: #cb1d2e, Unknown: #262262 




# 8-2) Set up color ----
genetics_palette <- c(
  "SFTPB"     = "#f38187",
  "Control"   = "#b6dbab",
  "Control_p" = "#009444",
  "ABCA3"     = "#80a6cc",
  "FARS2"     = "#8cc6ec",
  "GRN"       = "#e7eb94",
  "LRBA"      = "#c75186",
  "NFKB_NOD2" = "#fff450",
  "NLRP12"    = "#f9cbdf",
  "SFTPC"     = "#d9dad9",
  "SLC7A7"    = "#59c1ad",
  "SOCS1"     = "#fec369",
  "STAT1"     = "#782b90",
  "DICER"     = "#64735a",
  "FLNA"      = "#faf8cd",
  "IKBKB"     = "#c2b59b",
  "NPC2"      = "#ec008c",
  "PIK3CA"    = "#cb1d2e",
  "Unknown"   = "#262262"
)

# sanity check (optional): warn if any levels lack a color
check_levels <- function(x, palette) {
  missing <- setdiff(unique(x), names(palette))
  if (length(missing)) warning("Missing colors for: ", paste(missing, collapse=", "))
}




# 8-3) Re-plot ----
# add log2 + labels (as you had)
## 8-3-1) Re-do 'Are some genetics enriched for disagreement?' ----
library(dplyr)
fmtp  <- function(x) ifelse(x < 1e-16, "<1e-16", formatC(x, format="e", digits=2))
stars <- function(x) case_when(x < 0.001 ~ "***", x < 0.01 ~ "**", x < 0.05 ~ "*", TRUE ~ "")

# our existing fisher_df code (unchanged)
tab <- table(Status = czi_combined$L1_status, Genetics = czi_combined$genetics)
fisher_df <- lapply(colnames(tab), function(g){
  m <- matrix(c(tab["Disagree", g], tab["Agree", g],
                sum(tab["Disagree", ]) - tab["Disagree", g],
                sum(tab["Agree", ])    - tab["Agree", g]),
              nrow=2, byrow=TRUE)
  ft <- fisher.test(m)
  data.frame(
    genetics = g,
    OR  = as.numeric(ft$estimate),
    lo  = ft$conf.int[1],
    hi  = ft$conf.int[2],
    p   = ft$p.value,
    disagree_in_genetics = tab["Disagree", g],
    total_in_genetics    = sum(tab[, g])
  )
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p, "fdr"),
         rate_disagree = disagree_in_genetics / total_in_genetics) %>%
  tibble::as_tibble()

fisher_df_plot <- fisher_df %>%
  mutate(
    log2OR = log2(OR), log2lo = log2(lo), log2hi = log2(hi),
    p_lab = fmtp(p), FDR_lab = fmtp(FDR), star = stars(FDR),
    genetics = factor(genetics, levels = names(genetics_palette))
  )

check_levels(fisher_df_plot$genetics, genetics_palette)

library(ggplot2)
p_genetics <- ggplot(fisher_df_plot,
                     aes(x = reorder(genetics, log2OR), y = log2OR, color = genetics)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = log2lo, ymax = log2hi), width = 0.2) +
  geom_text(aes(label = star), nudge_y = 0.10, size = 4, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_color_manual(values = genetics_palette, breaks = names(genetics_palette)) +
  labs(y = "log2(OR)  (0 = no enrichment)", x = "Genetics",
       title = "Genetics-level disagreement enrichment",
       subtitle = "95% CI on log2 scale; stars = FDR",
       color = "Genetics") +
  theme_minimal(base_size = 13)
p_genetics



## 8-3-2) Per sample ----
# Build per-sample table (your code)
tab_sample <- table(Status = czi_combined$L1_status,
                    Sample = czi_combined$updated_sample_id)

# (Optional) continuity-corrected Fisher (guards zeros)
safe_fisher <- function(a,b,c,d){
  if(min(a,b,c,d) == 0){
    a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5
  }
  fisher.test(matrix(c(a,b,c,d), nrow=2, byrow=TRUE))
}

fisher_sample <- lapply(colnames(tab_sample), function(s){
  a <- tab_sample["Disagree", s]
  b <- tab_sample["Agree",    s]
  c <- sum(tab_sample["Disagree", ]) - a
  d <- sum(tab_sample["Agree",    ]) - b
  ft <- safe_fisher(a,b,c,d)
  data.frame(
    sample   = s,
    OR       = as.numeric(ft$estimate),
    lo       = ft$conf.int[1],
    hi       = ft$conf.int[2],
    p        = ft$p.value,
    disagree_in_sample = a,
    total_in_sample    = a + b
  )
}) %>% dplyr::bind_rows() %>%
  dplyr::mutate(FDR = p.adjust(p, "fdr"))

# Join genetics per sample
sample_meta <- czi_combined@meta.data %>%
  dplyr::select(updated_sample_id, genetics) %>% distinct()

fisher_sample <- fisher_sample %>%
  dplyr::left_join(sample_meta, by = c("sample" = "updated_sample_id"))

fisher_sample_plot <- fisher_sample %>%
  mutate(
    log2OR = log2(OR), log2lo = log2(lo), log2hi = log2(hi),
    p_lab = fmtp(p), FDR_lab = fmtp(FDR), star = stars(FDR),
    genetics = factor(genetics, levels = names(genetics_palette))
  )

check_levels(fisher_sample_plot$genetics, genetics_palette)

p_samples <- ggplot(fisher_sample_plot,
                    aes(x = reorder(sample, log2OR), y = log2OR, color = genetics)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = log2lo, ymax = log2hi), width = 0.2) +
  geom_text(aes(label = star), nudge_y = 0.10, size = 4, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  scale_color_manual(values = genetics_palette, breaks = names(genetics_palette)) +
  labs(y = "log2(OR)  (0 = no enrichment)", x = "Sample",
       title = "Sample-level disagreement enrichment",
       subtitle = "95% CI on log2 scale; stars = FDR",
       color = "Genetics") +
  theme_minimal(base_size = 13)
p_samples

# 13x8 l
p_genetics | p_samples






# additional: To get the number of cells of lvl 1 in agree and disagree ----
# Agree
obj_agree <- readRDS("output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/l1_temp_agreed_step13_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")
obj_agree # 
Idents(obj_agree) <- "L1_final"
agree_table <- table(Idents(obj_agree))
agree_table
write.csv(agree_table, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/l1_agree_cell_number.csv")
rm(obj_agree)
obj_agree_table <- table(obj_agree$updated_sample_id, obj_agree$L1_final)
obj_agree_table
write.csv(obj_agree_table, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/l1_agree_vs_samples_cell_number.csv")


# Disagree
obj_disagree <- readRDS("output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/obj_disagree_step13_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")
obj_disagree # 17152 nuclei

# Tristan
Idents(obj_disagree) <- "L1_tristan"
disagree_table <- table(Idents(obj_disagree))
disagree_table
write.csv(disagree_table, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/l1_disagree_tristan_cell_number.csv")

# Elizabeth
Idents(obj_disagree) <- "L1_elizabeth"
disagree_table <- table(Idents(obj_disagree))
disagree_table
write.csv(disagree_table, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/l1_disagree_elizabeth_cell_number.csv")







