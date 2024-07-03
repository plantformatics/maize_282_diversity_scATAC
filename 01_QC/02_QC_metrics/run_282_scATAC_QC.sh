#!/bin/bash

## submission properties

#SBATCH --partition=hugemem_p
#SBATCH --job-name=ANALYZE_scATAC_QC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=1-00:00:00
#SBATCH --mem=1750gb
#SBATCH --output=DO_282_maize_QC.%j.log
#SBATCH --error=DO_282_maize_QC.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marand@uga.edu

# set env
cd $SLURM_SUBMIT_DIR

# modules
module load parallel/20200422-GCCcore-8.3.0
module load R/4.0.0-foss-2019b
module load MACS2/2.2.7.1-foss-2019b-Python-3.7.4

# global vars
threads=30

# input variables
bedfiles=$(ls /scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/bamfiles/raw_bam/pool*.bed.gz )

# functions
runQC(){

	# input
	bed=$1
	id=$( basename $bed )
	name=$( echo $id | cut -d'.' -f1 )
	ref=/scratch/apm25309/reference_genomes/Zmays/v5/Zm-B73-REFERENCE-NAM-5.0.fa.fai.sorted
	ann=/scratch/apm25309/reference_genomes/Zmays/v5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.sorted.gff3
	
	# verbose
	echo " - running Socrates scATAC-seq QC analysis for $name ..."

	# run 
	Rscript QC_scATAC_data.v2.R $bed $name $ann $ref
	Rscript filter_low_qual_cells.v2.R $name.raw.soc.rds $name

}
export -f runQC

parallel -j $threads runQC {} ::: $bedfiles
