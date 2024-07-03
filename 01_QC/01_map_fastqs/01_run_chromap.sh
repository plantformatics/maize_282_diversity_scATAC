#!/bin/bash

## submission properties

#SBATCH --partition=batch
#SBATCH --job-name=chromap_pool1-40
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=7-00:00:00
#SBATCH --mem=115g
#SBATCH --output=DO_chromap_pool1-2.%j.log
#SBATCH --error=DO_chromap_pool1-2.%j.err
#SBATCH --array=1-40

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc

# load modules

# variables
fqdir=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/raw_data/fastq

# functions
mapReads(){

	# vars
	id=$1
	fqdir=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/raw_data/fastq
	index=/scratch/apm25309/reference_genomes/Zmays/v5/Zm-B73-REFERENCE-NAM-5.0.index
	ref=/scratch/apm25309/reference_genomes/Zmays/v5/Zm-B73-REFERENCE-NAM-5.0.fa
	whitelist=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/bamfiles/raw_bam/737K-cratac-v1.txt
	read1=$fqdir/pool${id}_R1_001.fastq.gz
	read2=$fqdir/pool${id}_R3_001.fastq.gz
	bc=$fqdir/pool${id}_R2_001.fastq.gz
	out=pool${id}

	# run
	chromap -l 2000 \
	-q 1 \
	--remove-pcr-duplicates-at-cell-level \
	--low-mem \
	--trim-adapters \
	-x $index \
	-r $ref \
	-1 $read1 \
	-2 $read2 \
	-o $out.raw.bam \
	-b $bc \
	--SAM \
	--barcode-whitelist $whitelist \
	--bc-error-threshold 2 \
	-t 48 

}
export -f mapReads

# run
mapReads $SLURM_ARRAY_TASK_ID

