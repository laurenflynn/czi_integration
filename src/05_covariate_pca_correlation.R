# Load libraries ----
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


#########
test_pc_metadata_association <- function(seurat_obj, pc_dim = 20, meta_vars = NULL) {
  # Extract metadata and PCA embeddings
  meta <- seurat_obj@meta.data
  pcs <- Embeddings(seurat_obj, "pca")[, 1:pc_dim]
  colnames(pcs) <- paste0("PC", 1:pc_dim)
  
  # Use all metadata columns if none specified
  if (is.null(meta_vars)) meta_vars <- colnames(meta)
  
  # Filter valid metadata variables
  meta_vars <- intersect(meta_vars, colnames(meta))
  if (length(meta_vars) == 0) stop("No valid metadata variables found.")
  
  # Initialize result matrices
  pval_matrix <- matrix(NA, nrow = length(meta_vars), ncol = pc_dim)
  effect_matrix <- matrix(NA, nrow = length(meta_vars), ncol = pc_dim)
  rownames(pval_matrix) <- rownames(effect_matrix) <- meta_vars
  colnames(pval_matrix) <- colnames(effect_matrix) <- colnames(pcs)
  
  for (meta_var in meta_vars) {
    var_data <- meta[[meta_var]]
    
    for (i in 1:pc_dim) {
      pc_vals <- pcs[, i]
      
      if (all(is.na(var_data))) next
      
      if (is.numeric(var_data)) {
        # Pearson correlation
        test_result <- cor.test(var_data, pc_vals)
        r <- test_result$estimate
        r2 <- r^2  # effect size = R-squared
        pval_matrix[meta_var, i] <- test_result$p.value
        effect_matrix[meta_var, i] <- r2
      } else {
        # ANOVA
        df <- data.frame(pc = pc_vals, group = as.factor(var_data))
        aov_model <- aov(pc ~ group, data = df)
        aov_summary <- summary(aov_model)
        pval_matrix[meta_var, i] <- aov_summary[[1]]$`Pr(>F)`[1]
        
        # Effect size = η² = SSB / SST
        ssb <- aov_summary[[1]]$`Sum Sq`[1]  # Between-group SS
        sst <- sum(aov_model$residuals^2) + ssb  # Total SS
        eta2 <- ssb / sst
        effect_matrix[meta_var, i] <- eta2
      }
    }
  }
  
  return(list(pval = pval_matrix, effect_size = effect_matrix))
}

# Run correlation func ----
results <- test_pc_metadata_association(czi_combined, pc_dim = 30,
                                        meta_vars = c("percent.mt", "Round", "Library_conc_nM", "Source", 
                                                      "DIP_or_Lympho", "nCount_RNA", "Ancestry", "Sample_type",
                                                      "Lobe", "Sex", "PCR_cycle", "Age", 
                                                      "B_Location_L_or_R", "Genetics", "updated_sampleID", "Sample_type"))

pvals <- results$pval
effects <- results$effect_size

# Visualize -log10(p) -----
# Visualize -log10(p-values)
logp <- -log10(pvals)
logp[is.infinite(logp)] <- max(logp[is.finite(logp)], na.rm = TRUE)

# 4.5x10 l
# Heatmap with Ward's hierarchical clustering
library(pheatmap)
pheatmap(
  logp,
  main = "-log10(p-values): Metadata-PC Associations",
  clustering_method = "ward.D2",        # This sets Ward clustering
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(100)
)


# Visualize effect sizes (R² or η²) ----
pheatmap(effects, main = "Effect Sizes (R² or η²)",
         clustering_method = "ward.D2",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "blue"))(100)
)

# Heatmap with Ward's hierarchical clustering ----
library(pheatmap)
pheatmap(
  logp,
  main = "-log10(p-values): Metadata-PC Associations",
  clustering_method = "ward.D2",        # This sets Ward clustering
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(100)
)


## Visualize effect sizes (R² or η²) ----
pheatmap(effects, main = "Effect Sizes (R² or η²)",
         clustering_method = "ward.D2",
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("white", "blue"))(100)
)



## Use Euclidean distance and Ward's method ----
dist_matrix <- dist(logp_matrix, method = "euclidean")
hc <- hclust(dist_matrix, method = "ward.D2")  # "ward.D2" is recommended over "ward.D"

# Plot dendrogram
plot(hc, main = "Hierarchical Clustering of Metadata Variables (Ward.D2)", cex = 0.8)


