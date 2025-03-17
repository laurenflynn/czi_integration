#! /bin/bash

#SBATCH --partition=bch-compute
#SBATCH --time=04:00:00             # hours-minutes-seconds
#SBATCH --mem=300G
#SBATCH --job-name=integration
#SBATCH --output=output/logs/integration_%J.out       
#SBATCH --error=output/logs/integration_%J.err    
#SBATCH --nodes=1
#SBATCH --ntasks=1

source /programs/biogrids.shrc
export R_X=4.4.2

Rscript src/01_integration.R