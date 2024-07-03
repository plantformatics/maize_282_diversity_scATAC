#!/bin/bash

## submission properties

#SBATCH --partition=hugemem_p
#SBATCH --job-name=UMAP_markers
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=7-00:00:00
#SBATCH --mem=950g
#SBATCH --output=ASSESS_marker_accessibility.%j.log
#SBATCH --error=ASSESS_marker_accessibility.%j.err

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc

# vars
sparseGB=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/sparse/genes/all_pools.raw.genes.sparse.rds
all_rd=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step0_QC_data/maize_282.v8.4.NOHARMONY_CELLs.reduced_dimensions.txt
all_meta=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step0_QC_data/maize_282.v8.4.NOHARMONY_CELLs.metadata.txt

# run for each
Rscript plot_marker_accessibility.R $all_meta $sparseGB $all_rd markers.filtered.v5.bed 30 maize_282.v8.4.NOHARMONY_CELLs
