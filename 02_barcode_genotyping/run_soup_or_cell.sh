#!/bin/bash

## submission properties
#SBATCH --partition=highmem_p
#SBATCH --job-name=INFORMED_souporcell_pool
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=00-16:00:00
#SBATCH --mem=490g
#SBATCH --output=DO_INFORMED_souporcell.%j.log
#SBATCH --error=DO_INFORMED_souporcell.%j.err

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc

# load modules
conda activate souporcell
ml parallel/20200422-GCCcore-8.3.0

# inputs
threads=3
bamdir=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/bamfiles/raw_bam
vcfdir=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/vcf_282
vcffile=$vcfdir/282_maize.biallelic.maf05.172.GP.dp_filt.PHASED.unzipped.vcf

# functions
runDemux(){

	# inputs
	i=$1
	bamdir=$2
	vcf=$3

	# new inputs
	bam=$bamdir/pool$i.updatedBC.mq10.rmdup.bam
	id=pool$i
	out=$4/$id

	# sample info
	samples=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/vcf_282/genotypes_173exp.B73_AGPv5.info.10_2_21.txt
	meta=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/bamfiles/raw_bam/All_pools.metadata.updated_BCs.500.txt
	fasta=/scratch/apm25309/reference_genomes/Zmays/v5/Zm-B73-REFERENCE-NAM-5.0.fa

	# make output dir
	if [ ! -d $out ]; then
		mkdir $out
	fi

	# subset genotypes
	if [ ! -f $out/$id.INPUT.genotypes.txt ]; then

		if [ $id == "pool21" ] || [ $id == "pool22" ]; then
			awk -F'\t' -v poolID='pool1' '$2==poolID' $samples | cut -f1 - | uniq - > $out/$id.INPUT.genotypes.txt
		elif [ $id == "pool40" ]; then
			awk -F'\t' -v poolID='pool39_mix' '$2==poolID' $samples | cut -f1 - | uniq - > $out/$id.INPUT.genotypes.txt
		else
			awk -F'\t' -v poolID=$id '$2==poolID' $samples | cut -f1 - | uniq - > $out/$id.INPUT.genotypes.txt
		fi
	fi

	# subset barcodes
	if [ ! -f $out/$id.BARCODES.txt ]; then
		awk -F'\t' -v poolID=$id '$7==poolID' $meta | cut -f2 - | uniq - > $out/$id.BARCODES.txt
		sed -i 's/CB:Z://' $out/$id.BARCODES.txt
	fi

	# filter VCF
	if [ ! -f $out/$id.filtered.vcf ]; then
		perl filterVCF.v3.pl $vcf $out/$id.INPUT.genotypes.txt > $out/$id.filtered.vcf
	fi

	# all genotypes
	cut -f1 $samples | sort -k1,1 - | uniq - | grep -v "Genotype" - > $out/all_genotypes.txt
	
	# save genotypes to var
	genos=$(while read i; do
		echo -n "$i "
		done < $out/$id.INPUT.genotypes.txt)
	num_k=$(wc -l < $out/$id.INPUT.genotypes.txt)
	echo " - k = $num_k | $genos"

	# run demuxlet -- change call rate to 0.9
	python ~/Software/souporcell/souporcell_pipeline.py -i $bam \
		--fasta $fasta \
		--threads 30 \
		--out_dir $out \
		-k $num_k \
		--min_alt 5 \
		--min_ref 5 \
		--skip_remap True \
		--known_genotypes $out/$id.filtered.vcf \
		--known_genotypes_sample_names $genos \
		-b $out/$id.BARCODES.txt \
		--no_umi True 

	# allow all genotypes
#	python ~/Software/souporcell/souporcell_pipeline.py -i $bam \
#                --fasta $fasta \
#                --threads 30 \
#                --out_dir ${out}_FDR \
#                -k $all_k \
#		--min_alt 5 \
#		--min_ref 5 \
#                --skip_remap True \
#                --known_genotypes $vcf \
#		--known_genotypes_sample_names $allgenos \
#                -b $out/$id.BARCODES.txt \
#		--no_umi True

}
export -f runDemux

# run demuxlet
pools=$( seq 1 40 )
#parallel -j $threads runDemux {} $bamdir $vcffile $PWD ::: $pools
runDemux 40 $bamdir $vcffile $PWD
