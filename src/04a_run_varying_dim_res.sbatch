#!/bin/bash

#SBATCH --partition=bch-compute          
#SBATCH --time=1:00:00         
#SBATCH --job-name=umap_vary_dim_res        
#SBATCH --output=output/logs/03_04_output_%J.txt      
#SBATCH --error=output/logs/03_04_error_%J.err # Name of the error file
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks=1              # Number of cpus
#SBATCH --mem=24G


source /programs/biogrids.shrc
export R_X=4.4.2

for dim in 15 20 25 30
do
  echo "Number of dimensions: $dim"
  for res in 0.5 0.25
  do
    echo "Resolution: $res"
    Rscript src/03_clustering.R  --dim $dim --res $res
    Rscript src/04_dim_reduction.R  --dim $dim --res $res
  done
done
