#!/bin/bash

## submission properties

#SBATCH --partition=highmem_p
#SBATCH --job-name=pseudobulk_DAR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=00-10:00:00
#SBATCH --mem=300g
#SBATCH --output=DO_DAR_detection.%j.log
#SBATCH --error=DO_DAR_detection.%j.err

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc

# vars
meta=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step0_QC_data/maize_282.v8.3.ALL_CELLs.metadata.leiden.pruned.subclustered.txt
sparse=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/sparse/genes/all_pools.raw.genes.sparse.rds
markers=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step3_differential_accessibility/genes/markers.filtered.v5.bed

# run
Rscript pseudobulk_DAR_analysis.v2.R $meta $sparse $markers
Rscript PlotClusterZscore.R $meta $sparse $markers
