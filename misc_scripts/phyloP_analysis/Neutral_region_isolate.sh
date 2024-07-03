#!/bin/bash
#SBATCH --job-name=Neutral_regions       # Job name
#SBATCH --partition=batch               # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=48                  # Number of CPU cores per task
#SBATCH --mem=100gb                          # Job memory request
#SBATCH --time=96:00:00                     # Time limit hrs:min:sec
#SBATCH --output=nt_region_ID.%j.out         # Standard output log
#SBATCH --error=nt_region_ID.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=john.mendieta@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

cd $SLURM_SUBMIT_DIR

source /apps/lmod/lmod/init/zsh

ml Python/3.8.2-GCCcore-8.3.
ml STAR
ml BEDTools
ml SAMtools 
ml bioawk 
ml Trimmomatic
ml snakemake
ml BLAST+
ml ucsc

snakemake -rs ID_neutral_regions.snake --configfile ID_neutral_regions_config.yaml --cores 48 #-R take_bed_depth_hist



