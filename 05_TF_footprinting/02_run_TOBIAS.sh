#!/bin/bash

## submission properties

#SBATCH --partition=iob_p
#SBATCH --job-name=TOBIAS_footprinting
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=7-00:00:00
#SBATCH --mem=500gb
#SBATCH --output=LOG_tobias_preprocessing.%j.log
#SBATCH --error=LOG_tobias_preprocessing.%j.err
#SBATCH --array=1-20

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc

# modules
ml GLib/2.75.0-GCCcore-12.2.0
source .venv/bin/activate

# files
config=celltypes.txt

# functions
TOBIAScor(){

	# verbose
        echo "- running TOBIAS ATACorrect for $1..."

	# variables
	celltype=$1
	dir=$2/$celltype
	acrs=footprinting_ACRs.bed
	ref1=/scratch/apm25309/reference_genomes/Zmays/v5/Zm-B73-REFERENCE-NAM-5.0.fa.hap1
	ref2=/scratch/apm25309/reference_genomes/Zmays/v5/Zm-B73-REFERENCE-NAM-5.0.fa.hap2

	# run tobias
	TOBIAS ATACorrect --bam $dir/$celltype.haplotype1.bam \
		--genome $ref1 \
		--peaks $acrs \
		--cores 30 \
		--outdir $dir

	TOBIAS ATACorrect --bam $dir/$celltype.haplotype2.bam \
                --genome $ref2 \
                --peaks $acrs \
                --cores 30 \
                --outdir $dir
}
export -f TOBIAScor

TOBIASscore(){

	# verbose
	echo "- running TOBIAS footprint scores for $1..."

	# variables
	celltype=$1
        dir=$2/$celltype
        acrs=footprinting_ACRs.bed

	# run footprint scores
	TOBIAS ScoreBigwig --signal $dir/$celltype.haplotype1_corrected.bw \
		--regions $acrs\
		--output $dir/${celltype}.haplotype1_footprints.bw \
		--cores 30 

        TOBIAS ScoreBigwig --signal $dir/$celltype.haplotype2_corrected.bw \
                --regions $acrs\
                --output $dir/${celltype}.haplotype2_footprints.bw \
                --cores 30


}
export -f TOBIASscore

TOBIASbin(){

	        # verbose
        echo "- running TOBIAS BINDetect for $1..."

        # variables
        celltype=$1
        dir=$2/$celltype
        motifs=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step11_allelic_ATAC_footprinting/JASPAR2022_CORE_plants_redundant_pfms_jaspar.jaspar
	ref=/scratch/apm25309/reference_genomes/Zmays/v5/Zm-B73-REFERENCE-NAM-5.0.fa
	acrs=footprinting_ACRs.bed

	# run TFBS prediction
	TOBIAS BINDetect --motifs $motifs \
		--signals $dir/${celltype}.haplotype1_footprints.bw $dir/${celltype}.haplotype2_footprints.bw \
		--genome $ref \
		--peaks $acrs \
		--outdir $dir/BINDetect_output \
		--cond_names $celltype.hap1 $celltype.hap2 \
		--cores 30
		
}
export -f TOBIASbin

# iterate over each cell types
sID=$(awk -v ArrayID=$SLURM_ARRAY_TASK_ID '$2==ArrayID {print $1}' $config)
TOBIAScor $sID $PWD
TOBIASscore $sID $PWD
TOBIASbin $sID $PWD
