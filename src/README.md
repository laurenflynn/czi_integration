## Files

- 00_generate_metadata.R Generates output/00_metadata.csv; all metadata has been verified (excluded from Git)

- 01_integration.R Integrates biopsy samples that have passed QC, generates output/01_integrated_biopsy_samples.rds

- 01a_run_integration.sbatch Slurm script to run 01_integration.R

- 02_add_metadata_post_integration.R Adds updated metadata to integrated R object, generates 02_integrated_biopsy_samples_metadata_updated.rds

- 03_post_integration_pca.R Runs principal component analysis on integrated data, generates 03_post_integration_pca.rds

- 04_clustering.R Uses KNN to find clusters, output in output/04_clustering

- 05_dim_reduction.R Generates UMAPs, output in output/05_dim_reduction/ and output/05_figures/

- 05a_run_varying_dim_res.sbatch Runs 04_clustering.R and 05_dim_reduction.R

