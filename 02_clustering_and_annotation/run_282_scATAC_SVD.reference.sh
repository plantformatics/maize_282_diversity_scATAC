#!/bin/bash

## submission properties

#SBATCH --partition=highmem_p
#SBATCH --job-name=REFERENCE_scATAC_cluster
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-06:00:00
#SBATCH --mem=450gb
#SBATCH --output=DO_282_maize_REFERENCE_cluster.%j.log
#SBATCH --error=DO_282_maize_REFERENCE_cluster.%j.err

# set env
source ~/.zshrc
cd $SLURM_SUBMIT_DIR

# modules
module load R/4.0.0-foss-2019b

# run
Rscript SVD_analysis_scATAC.reference.R
