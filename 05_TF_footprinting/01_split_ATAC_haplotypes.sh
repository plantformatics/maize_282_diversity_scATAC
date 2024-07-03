#!/bin/bash

## submission properties

#SBATCH --partition=highmem_p
#SBATCH --job-name=partition_BAM_4_footprinting
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=500gb
#SBATCH --output=LOG_tobias_preprocessing.%j.log
#SBATCH --error=LOG_tobias_preprocessing.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=marand@uga.edu
#SBATCH --array=1-20

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc

# files
config=celltypes.txt

# functions
processACR(){

	# variables
	celltype=$1
	outfile=$2/$celltype/$celltype
	acrs=footprinting_ACRs.bed
	bam=all_pools.updatedBC.mq10.rmdup.caQTL-ACRs.bam
	meta=maize_282.v8.3.ALL_CELLs.metadata.leiden.pruned.subclustered.txt
	groups=ACR_haplotype_groups.txt

	# create output dir
	if [ ! -d $celltype ]; then
                mkdir $celltype
        else
                rm -rf $celltype
                mkdir $celltype
        fi

	# split metadata
	awk -v CTid=$celltype '$26==CTid {print;}' $meta > $meta.$celltype

	# stream
	samtools view -L $acrs $bam | 
		perl splitBAMs.pl $meta.$celltype $groups $outfile - 

	# clean-up
	rm $meta.$celltype

}

# iterate over each cell types
sID=$(awk -v ArrayID=$SLURM_ARRAY_TASK_ID '$2==ArrayID {print $1}' $config)
processACR $sID $PWD
