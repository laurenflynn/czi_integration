# Step1 Create Seurat Object for each sample of round 2





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
dataPath <- "/Users/seunghanbaek/projects/snRNAseq_CZI_pipeline/data/round2_seq_data_10xformat/"
FL_22NA_72624 <- Read10X(data.dir = paste(dataPath, "L-22NA-72624/", sep = ""))
FL_24_8724 <- Read10X(data.dir = paste(dataPath, "L-24-8724/", sep = ""))
FL_24_72924 <- Read10X(data.dir = paste(dataPath, "L-24-72924/", sep = ""))
FL_25NA_72624 <- Read10X(data.dir = paste(dataPath, "L-25NA-72624/", sep = ""))
FL_28_72524 <- Read10X(data.dir = paste(dataPath, "L-28-72524/", sep = ""))
FL_28_72624 <- Read10X(data.dir = paste(dataPath, "L-28-72624/", sep = ""))
FL_33_72924 <- Read10X(data.dir = paste(dataPath, "L-33-72924/", sep = ""))
FL_53_72924 <- Read10X(data.dir = paste(dataPath, "L-53-72924/", sep = ""))
FL_59_72624 <- Read10X(data.dir = paste(dataPath, "L-59-72624/", sep = ""))
FL_62_72524 <- Read10X(data.dir = paste(dataPath, "L-62-72524/", sep = ""))
FL_62_72924 <- Read10X(data.dir = paste(dataPath, "L-62-72924/", sep = ""))
FL_167_72624 <- Read10X(data.dir = paste(dataPath, "L-167-72624/", sep = ""))
FL_283_72924 <- Read10X(data.dir = paste(dataPath, "L-283-72924/", sep = ""))
FL_313_72924 <- Read10X(data.dir = paste(dataPath, "L-313-72924/", sep = ""))
FL_314_8724 <- Read10X(data.dir = paste(dataPath, "L-314-8724/", sep = ""))
FL_315_8724 <- Read10X(data.dir = paste(dataPath, "L-315-8724/", sep = ""))
FL_316_8724 <- Read10X(data.dir = paste(dataPath, "L-316-8724/", sep = ""))
FL_316_72524 <- Read10X(data.dir = paste(dataPath, "L-316-72524/", sep = ""))
FL_319_8724 <- Read10X(data.dir = paste(dataPath, "L-319-8724/", sep = ""))
FL_320_72924 <- Read10X(data.dir = paste(dataPath, "L-320-72924", sep = ""))
FL_321_72624 <- Read10X(data.dir = paste(dataPath, "L-321-72624/", sep = ""))
FL_323_8724 <- Read10X(data.dir = paste(dataPath, "L-323-8724/", sep = ""))
FL_154_72524 <- Read10X(data.dir = paste(dataPath, "L-154-72524/", sep = ""))





# 2. Create Seurat Object
# min.cells = 3 : include genes detected in at least 3 nuclei (cells).
# min.features = 100 : include nucleus (cells) with at lease 100 genes detected.
FL_22NA_72624 <- CreateSeuratObject(counts = FL_22NA_72624, project = "FL_22NA_72624", min.cells = 3, min.features = 100)
FL_22NA_72624 # 52688 nuclei

FL_24_8724 <- CreateSeuratObject(counts = FL_24_8724, project = "FL_24_8724", min.cells = 3, min.features = 100)
FL_24_8724 # 7389 nuclei

FL_24_72924 <- CreateSeuratObject(counts = FL_24_72924, project = "FL_24_72924", min.cells = 3, min.features = 100)
FL_24_72924 # 3924 nuclei

FL_25NA_72624 <- CreateSeuratObject(counts = FL_25NA_72624, project = "FL_25NA_72624", min.cells = 3, min.features = 100)
FL_25NA_72624 # 7054 nuclei

FL_28_72524 <- CreateSeuratObject(counts = FL_28_72524, project = "FL_28_72524", min.cells = 3, min.features = 100)
FL_28_72524 # 6938 nuclei

FL_28_72624 <- CreateSeuratObject(counts = FL_28_72624, project = "FL_28_72624", min.cells = 3, min.features = 100)
FL_28_72624 # 10226 nuclei

FL_33_72924 <- CreateSeuratObject(counts = FL_33_72924, project = "FL_33_72924", min.cells = 3, min.features = 100)
FL_33_72924 # 2041 nuclei

FL_53_72924 <- CreateSeuratObject(counts = FL_53_72924, project = "FL_53_72924", min.cells = 3, min.features = 100)
FL_53_72924 # 2587 nuclei

FL_59_72624 <- CreateSeuratObject(counts = FL_59_72624, project = "FL_59_72624", min.cells = 3, min.features = 100)
FL_59_72624 # 3123 nuclei

FL_62_72524 <- CreateSeuratObject(counts = FL_62_72524, project = "FL_62_72524", min.cells = 3, min.features = 100)
FL_62_72524 # 8790 nuclei

FL_62_72924 <- CreateSeuratObject(counts = FL_62_72924, project = "FL_62_72924", min.cells = 3, min.features = 100)
FL_62_72924 # 6776 nuclei

FL_167_72624 <- CreateSeuratObject(counts = FL_167_72624, project = "FL_167_72624", min.cells = 3, min.features = 100)
FL_167_72624 # 4614 nuclei

FL_283_72924 <- CreateSeuratObject(counts = FL_283_72924, project = "FL_283_72924", min.cells = 3, min.features = 100)
FL_283_72924 # 15046 nuclei

FL_313_72924 <- CreateSeuratObject(counts = FL_313_72924, project = "FL_313_72924", min.cells = 3, min.features = 100)
FL_313_72924 # 6922 nuclei

FL_314_8724 <- CreateSeuratObject(counts = FL_314_8724, project = "FL_314_8724", min.cells = 3, min.features = 100)
FL_314_8724 # 4367 nuclei

FL_315_8724 <- CreateSeuratObject(counts = FL_315_8724, project = "FL_315_8724", min.cells = 3, min.features = 100)
FL_315_8724 # 9037 nuclei

FL_316_72524 <- CreateSeuratObject(counts = FL_316_72524, project = "FL_316_72524", min.cells = 3, min.features = 100)
FL_316_72524 # 59587 nuclei

FL_321_72624 <- CreateSeuratObject(counts = FL_321_72624, project = "FL_321_72624", min.cells = 3, min.features = 100)
FL_321_72624 # 3154 nuclei

FL_323_8724 <- CreateSeuratObject(counts = FL_323_8724, project = "FL_323_8724", min.cells = 3, min.features = 100)
FL_323_8724 # 12973 nuclei

FL_154_72524 <- CreateSeuratObject(counts = FL_154_72524, project = "FL_154_72524", min.cells = 3, min.features = 100)
FL_154_72524 # 7174 nuclei

FL_316_8724 <- CreateSeuratObject(counts = FL_316_8724, project = "FL_316_8724", min.cells = 3, min.features = 100)
FL_316_8724 # 20401 nuclei

FL_319_8724 <- CreateSeuratObject(counts = FL_319_8724, project = "FL_319_8724", min.cells = 3, min.features = 100)
FL_319_8724 # 9737 nuclei

FL_320_72924 <- CreateSeuratObject(counts = FL_320_72924, project = "FL_320_72924", min.cells = 3, min.features = 100)
FL_320_72924 # 20348 nuclei





# 3. Save Data
saveRDS(FL_22NA_72624, file = "output/round_2/output_4_step1_createSeuratObjects/FL_22NA_72624.rds")
saveRDS(FL_24_8724, file = "output/round_2/output_4_step1_createSeuratObjects/FL_24_8724.rds")
saveRDS(FL_24_72924, file = "output/round_2/output_4_step1_createSeuratObjects/FL_24_72924.rds")
saveRDS(FL_25NA_72624, file = "output/round_2/output_4_step1_createSeuratObjects/FL_25NA_72624.rds")
saveRDS(FL_28_72524, file = "output/round_2/output_4_step1_createSeuratObjects/FL_28_72524.rds")
saveRDS(FL_28_72624, file = "output/round_2/output_4_step1_createSeuratObjects/FL_28_72624.rds")
saveRDS(FL_33_72924, file = "output/round_2/output_4_step1_createSeuratObjects/FL_33_72924.rds")
saveRDS(FL_53_72924, file = "output/round_2/output_4_step1_createSeuratObjects/FL_53_72924.rds")
saveRDS(FL_59_72624, file = "output/round_2/output_4_step1_createSeuratObjects/FL_59_72624.rds")
saveRDS(FL_62_72524, file = "output/round_2/output_4_step1_createSeuratObjects/FL_62_72524.rds")
saveRDS(FL_62_72924, file = "output/round_2/output_4_step1_createSeuratObjects/FL_62_72924.rds")
saveRDS(FL_167_72624, file = "output/round_2/output_4_step1_createSeuratObjects/FL_167_72624.rds")
saveRDS(FL_283_72924, file = "output/round_2/output_4_step1_createSeuratObjects/FL_283_72924.rds")
saveRDS(FL_313_72924, file = "output/round_2/output_4_step1_createSeuratObjects/FL_313_72924.rds")
saveRDS(FL_314_8724, file = "output/round_2/output_4_step1_createSeuratObjects/FL_314_8724.rds")
saveRDS(FL_315_8724, file = "output/round_2/output_4_step1_createSeuratObjects/FL_315_8724.rds")
saveRDS(FL_316_72524, file = "output/round_2/output_4_step1_createSeuratObjects/FL_316_72524.rds")
saveRDS(FL_321_72624, file = "output/round_2/output_4_step1_createSeuratObjects/FL_321_72624.rds")
saveRDS(FL_323_8724, file = "output/round_2/output_4_step1_createSeuratObjects/FL_323_8724.rds")
saveRDS(FL_154_72524, file = "output/round_2/output_4_step1_createSeuratObjects/FL_154_72524.rds")
saveRDS(FL_316_8724, file = "output/round_2/output_4_step1_createSeuratObjects/FL_316_8724.rds")
saveRDS(FL_319_8724, file = "output/round_2/output_4_step1_createSeuratObjects/FL_319_8724.rds")
saveRDS(FL_320_72924, file = "output/round_2/output_4_step1_createSeuratObjects/FL_320_72924.rds")

# End of Step 1





