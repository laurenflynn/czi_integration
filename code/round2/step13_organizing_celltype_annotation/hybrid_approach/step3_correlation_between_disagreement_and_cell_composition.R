# Step 3 correlation betwee disaagreement and cell type composition ----





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
czi_combined <- readRDS("output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step1_find_common_in_lvl1/czi_combined_step13_harmony_nCountRNA_with_hlca_tristion_elizabeth_labels_lognormalized.rds")
czi_combined # 341046 nuclei
Idents(czi_combined) <- "genetics"
table(Idents(czi_combined))





# 2. Process Data ----
# Set up meta data
meta <- as.data.frame(czi_combined@meta.data)
# Set up color
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




# 2-1) per-sample disagreement and Epithelial proportion (agree only) ----
meta_sample_test <- meta %>%
  group_by(updated_sample_id, genetics) %>%
  summarise(
    total_cells     = n(),
    disagree_cells  = sum(L1_status == "Disagree"),
    rate_disagree   = disagree_cells / total_cells,
    agree_total     = sum(L1_status == "Agree"),
    agree_epi_cells = sum(L1_final == "Epithelial", na.rm = TRUE),
    pct_epi_agree   = ifelse(agree_total > 0, agree_epi_cells / agree_total, NA_real_),
    .groups = "drop"
  )
# Spearman correlation
ct <- cor.test(meta_sample_test$pct_epi_agree, meta_sample_test$rate_disagree,
               method = "spearman", use = "pairwise.complete.obs")
ct # -> use ct$estimate for rho and ct$p.value in our caption
# Plot  5x6.3 l
p1 <- ggplot(meta_sample_test,
       aes(x = pct_epi_agree * 100, y = rate_disagree * 100, color = genetics)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = genetics_palette, drop = FALSE) +
  labs(
    x = "Epithelial (%) within agreed cells",
    y = "Disagreement rate (%)",
    title = "Epithelial composition vs disagreement (per sample)",
    subtitle = paste0("Spearman’s rho = ",
                      sprintf("%.2f", unname(ct$estimate)),
                      ", p = ", formatC(ct$p.value, format = "g", digits = 3)),
    color = "Genetics"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
p1




# 2-2) per-sample disagreement and Endothelial proportion (agree only) ----
meta_sample_test <- meta %>%
  group_by(updated_sample_id, genetics) %>%
  summarise(
    total_cells     = n(),
    disagree_cells  = sum(L1_status == "Disagree"),
    rate_disagree   = disagree_cells / total_cells,
    agree_total     = sum(L1_status == "Agree"),
    agree_epi_cells = sum(L1_final == "Endothelial", na.rm = TRUE),
    pct_epi_agree   = ifelse(agree_total > 0, agree_epi_cells / agree_total, NA_real_),
    .groups = "drop"
  )
# Spearman correlation
ct <- cor.test(meta_sample_test$pct_epi_agree, meta_sample_test$rate_disagree,
               method = "spearman", use = "pairwise.complete.obs")
ct # -> use ct$estimate for rho and ct$p.value in our caption
# Plot  5x6.3 l
p2 <- ggplot(meta_sample_test,
             aes(x = pct_epi_agree * 100, y = rate_disagree * 100, color = genetics)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = genetics_palette, drop = FALSE) +
  labs(
    x = "Endothelial (%) within agreed cells",
    y = "Disagreement rate (%)",
    title = "Endothelial composition vs disagreement (per sample)",
    subtitle = paste0("Spearman’s rho = ",
                      sprintf("%.2f", unname(ct$estimate)),
                      ", p = ", formatC(ct$p.value, format = "g", digits = 3)),
    color = "Genetics"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
p2




# 2-3) per-sample disagreement and Immune proportion (agree only) ----
meta_sample_test <- meta %>%
  group_by(updated_sample_id, genetics) %>%
  summarise(
    total_cells     = n(),
    disagree_cells  = sum(L1_status == "Disagree"),
    rate_disagree   = disagree_cells / total_cells,
    agree_total     = sum(L1_status == "Agree"),
    agree_epi_cells = sum(L1_final == "Immune", na.rm = TRUE),
    pct_epi_agree   = ifelse(agree_total > 0, agree_epi_cells / agree_total, NA_real_),
    .groups = "drop"
  )
# Spearman correlation
ct <- cor.test(meta_sample_test$pct_epi_agree, meta_sample_test$rate_disagree,
               method = "spearman", use = "pairwise.complete.obs")
ct # -> use ct$estimate for rho and ct$p.value in our caption
# Plot  5x6.3 l
p3 <- ggplot(meta_sample_test,
             aes(x = pct_epi_agree * 100, y = rate_disagree * 100, color = genetics)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = genetics_palette, drop = FALSE) +
  labs(
    x = "Immune (%) within agreed cells",
    y = "Disagreement rate (%)",
    title = "Immune composition vs disagreement (per sample)",
    subtitle = paste0("Spearman’s rho = ",
                      sprintf("%.2f", unname(ct$estimate)),
                      ", p = ", formatC(ct$p.value, format = "g", digits = 3)),
    color = "Genetics"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
p3




# 2-4) per-sample disagreement and Stroma proportion (agree only) ----
meta_sample_test <- meta %>%
  group_by(updated_sample_id, genetics) %>%
  summarise(
    total_cells     = n(),
    disagree_cells  = sum(L1_status == "Disagree"),
    rate_disagree   = disagree_cells / total_cells,
    agree_total     = sum(L1_status == "Agree"),
    agree_epi_cells = sum(L1_final == "Stroma", na.rm = TRUE),
    pct_epi_agree   = ifelse(agree_total > 0, agree_epi_cells / agree_total, NA_real_),
    .groups = "drop"
  )
# Spearman correlation
ct <- cor.test(meta_sample_test$pct_epi_agree, meta_sample_test$rate_disagree,
               method = "spearman", use = "pairwise.complete.obs")
ct # -> use ct$estimate for rho and ct$p.value in our caption
# Plot  5x6.3 l
p4 <- ggplot(meta_sample_test,
             aes(x = pct_epi_agree * 100, y = rate_disagree * 100, color = genetics)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = genetics_palette, drop = FALSE) +
  labs(
    x = "Stroma (%) within agreed cells",
    y = "Disagreement rate (%)",
    title = "Stroma composition vs disagreement (per sample)",
    subtitle = paste0("Spearman’s rho = ",
                      sprintf("%.2f", unname(ct$estimate)),
                      ", p = ", formatC(ct$p.value, format = "g", digits = 3)),
    color = "Genetics"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")
p4
# 5x25.2 L
p1 | p2 | p3 | p4





# 3. Dominance-Disagreement Analysis ----
meta <- as.data.frame(czi_combined@meta.data)
meta_sample_agree <- meta %>%
  group_by(updated_sample_id, genetics) %>%
  summarise(
    total_cells     = n(),
    disagree_cells  = sum(L1_status == "Disagree"),
    agree_cells     = sum(L1_status == "Agree"),
    # agreed-only proportions (guard against divide-by-zero)
    prop_Epithelial_agree  = ifelse(agree_cells > 0,
                                    sum(L1_status == "Agree" & L1_final == "Epithelial") / agree_cells, NA_real_),
    prop_Endothelial_agree = ifelse(agree_cells > 0,
                                    sum(L1_status == "Agree" & L1_final == "Endothelial") / agree_cells, NA_real_),
    prop_Stroma_agree      = ifelse(agree_cells > 0,
                                    sum(L1_status == "Agree" & L1_final == "Stroma") / agree_cells, NA_real_),
    prop_Immune_agree      = ifelse(agree_cells > 0,
                                    sum(L1_status == "Agree" & L1_final == "Immune") / agree_cells, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    rate_disagree = disagree_cells / total_cells,
    dominance     = pmax(prop_Epithelial_agree, prop_Endothelial_agree,
                         prop_Stroma_agree,     prop_Immune_agree, na.rm = TRUE)
  ) %>%
  # keep one row per sample and drop samples with no agreed cells
  distinct(updated_sample_id, .keep_all = TRUE) %>%
  filter(is.finite(dominance), is.finite(rate_disagree))

# Spearman correlation
ct <- cor.test(meta_sample_agree$dominance, meta_sample_agree$rate_disagree,
               method = "spearman", use = "pairwise.complete.obs")
ct
dplyr::count(meta_sample_agree, updated_sample_id) %>% filter(n > 1)  # should be empty
# 5x6.2 l
ggplot(meta_sample_agree, aes(x = dominance, y = rate_disagree, color = genetics)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = genetics_palette, drop = FALSE) +
  labs(
    x = "Dominance (max Level-1 fraction within agreed cells)",
    y = "Disagreement rate (%)",
    title = "Do samples dominated by one cell class disagree less?",
    subtitle = paste0("Spearman’s rho = ",
                      sprintf("%.2f", unname(ct$estimate)),
                      ", p = ", formatC(ct$p.value, format = "g", digits = 3)),
    color = "Genetics"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")





# 4. To get the cell composition table per sample ---- 
# 0) Force a plain data.frame ----
meta_df <- as.data.frame(czi_combined@meta.data)

# 1) Per-sample totals (context) ----
totals <- meta_df %>%
  dplyr::group_by(updated_sample_id, genetics) %>%
  dplyr::summarise(
    total_all     = dplyr::n(),
    total_agree   = sum(L1_status == "Agree"),
    rate_disagree = 100 * (1 - total_agree / total_all),
    .groups = "drop"
  )

# 2) Agreed-only composition; avoid select() by subsetting columns with base R ----
agree_long <- meta_df %>%
  dplyr::filter(L1_status == "Agree", !is.na(L1_final)) %>%
  dplyr::count(updated_sample_id, genetics, L1_final, name = "n_agree") %>%
  dplyr::group_by(updated_sample_id, genetics) %>%
  dplyr::mutate(pct = 100 * n_agree / sum(n_agree)) %>%
  dplyr::ungroup()

agree_comp <- tidyr::pivot_wider(
  agree_long[, c("updated_sample_id","genetics","L1_final","pct")],  # <-- no dplyr::select
  names_from  = L1_final,
  values_from = pct,
  values_fill = 0
)

# 3) Nice column names + merge with totals ----
agree_comp <- dplyr::rename_with(agree_comp, ~ paste0("pct_", .x),
                                 -c(updated_sample_id, genetics))

agree_comp_table <- totals %>%
  dplyr::left_join(agree_comp, by = c("updated_sample_id","genetics")) %>%
  dplyr::select(updated_sample_id, genetics, total_all, total_agree, rate_disagree,
                dplyr::starts_with("pct_")) %>%
  dplyr::arrange(genetics, updated_sample_id)

# sanity check: one row per sample
stopifnot(!any(duplicated(agree_comp_table$updated_sample_id)))
agree_comp_table
write.csv(agree_comp_table, "output/round_2/step13_organizing_celltype_annotation/hybrid_approach/output_4_step3_correlation_between_disagreement_and_cell_composition/agree_celltype_percentages_per_sample.csv", row.names = FALSE)






