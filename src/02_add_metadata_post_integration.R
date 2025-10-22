# Step 2 Add metadata to integrated R object

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
library(data.table)


input_file_name <- "output/01_integrated_biopsy_samples.rds"
metadata_file_name <- "output/00_metadata.csv"
output_file_name <- "output/02_integrated_biopsy_samples_metadata_updated.rds"

# 1. Import Data ----

czi_combined <- readRDS(input_file_name) # Not with any metadata
czi_combined # 341046 nuclei

metadata <- fread(metadata_file_name)



# 2. Update Meta Data ----
# Just in case checking
table(Idents(czi_combined))
Idents(czi_combined) <- "orig.ident"
table(Idents(czi_combined))

## Add updated sample ID ----
czi_combined$updated_sample_id <- Idents(czi_combined)
czi_combined$updated_sample_id <- plyr::mapvalues(
  czi_combined$updated_sample_id,
  from = metadata$old_sample_id,
  to = metadata$updated_sample_id
)
Idents(czi_combined) <- "updated_sample_id"
table(Idents(czi_combined))

# Map the names of the current ident to the new column of idents we are providing
# Loop through metadata
for (col in colnames(metadata)) {
  print(col)
  if (col != "updated_sample_id") {
    czi_combined[[col]] <- Idents(czi_combined)
    czi_combined[[col]] <- plyr:::mapvalues(
      czi_combined$updated_sample_id,
      from = metadata$updated_sample_id,
      to = metadata[[col]],
      warn_missing = FALSE
    )
  }
  if(is.numeric(metadata[[col]])){
    czi_combined@meta.data[[col]] <- as.numeric(as.character(czi_combined@meta.data[[col]]))
  }
  Idents(czi_combined) <- col
  print(table(Idents(czi_combined)))
}


# Save
saveRDS(czi_combined, file = output_file_name)

