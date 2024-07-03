#!/bin/bash

#SBATCH --partition=highmem_p
#SBATCH --job-name=process_sc_bams1
#SBATCH --ntasks=40
#SBATCH --time=7-00:00:00
#SBATCH --mem=490g
#SBATCH --output=LOGS_BATCH_1_PROCESS_BAMs.%j.log
#SBATCH --error=LOGS_BATCH_1_PROCESS_BAMs.%j.err

# set env
cd $SLURM_SUBMIT_DIR

# threads
threads=5
qual=10

# load modules
ml picard/2.16.0-Java-1.8.0_144
ml parallel/20200422-GCCcore-8.3.0

# set dir
source ~/.zshrc

# functions
doCall(){

	# input
	base=pool$1
	qual=$2
	threads=$3

	# fixing mate pair information
	echo "fixing read pairs ..."
	if [ ! -f $base.raw.csort.bam ]; then
		samtools sort -n $base.raw.bam > $base.raw.csort.bam
	fi
	samtools fixmate $base.raw.csort.bam $base.raw.csort.fm.bam
	samtools sort $base.raw.csort.fm.bam > $base.raw.sort.bam

	# filter for mapped reads only
	echo "retaining only mapped reads ..."
	samtools view -@ $threads -bhq $qual -f 3 $base.raw.sort.bam > $base.mq$qual.bam

	# fix header
	echo "fixing SAM header ..."
	samtools view -H $base.mq$qual.bam | perl editHeader.pl - > $base.newheader
	samtools reheader $base.newheader $base.mq$qual.bam > $base.mq$qual.hdr.bam
	mv $base.mq$qual.hdr.bam $base.mq$qual.bam

	# run picard
	echo "removing dups - $base ..."
	java -Xmx100g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		REMOVE_DUPLICATES=true \
		METRICS_FILE=$base.metrics \
		I=$base.mq$qual.bam \
		O=$base.mq$qual.rmdup.bam \
		BARCODE_TAG=CB \
		ASSUME_SORT_ORDER=coordinate \
		USE_JDK_DEFLATER=true \
		USE_JDK_INFLATER=true

	# fix barcodes
	echo "fixing barcodes and removing multi-mapped reads ..."
	perl fixBC.pl $base.mq$qual.rmdup.bam $base | samtools view -bhS - > $base.BC.mq$qual.rmdup.bam

	# make Tn5 bed files
	echo "making Tn5 bed files ..."
	perl makeTn5bed.pl $base.BC.mq$qual.rmdup.bam | sort -k1,1 -k2,2n - | uniq - > $base.tn5.mq$qual.bed
	gzip $base.tn5.mq$qual.bed

}
export -f doCall

# set-up run
pools=$( seq 1 40 )

# run
parallel -j 2 doCall {} $qual $threads ::: $pools
