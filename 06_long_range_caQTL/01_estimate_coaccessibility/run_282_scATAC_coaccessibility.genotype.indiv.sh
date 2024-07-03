#!/bin/bash

## submission properties

#SBATCH --partition=highmem_p
#SBATCH --job-name=COACCESS_scATAC_cluster
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=6-23:00:00
#SBATCH --mem=450gb
#SBATCH --output=DO_282_maize_COACCESS_cluster.%j.log
#SBATCH --error=DO_282_maize_COACCESS_cluster.%j.err

# set env
source ~/.zshrc
cd $SLURM_SUBMIT_DIR

# modules

# run
Rscript co_accessibility_genotypes.indiv.R
Rscript eFDR_coACRs.R maize_282.co_accessibility.per_nucleus_genotypes.indiv.rds maize_282.co_accessibility.permuted.per_nucleus_genotypes.indiv.rds maize_282.per_nucleus_genotypes.indiv
