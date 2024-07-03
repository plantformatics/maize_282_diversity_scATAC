#!/bin/bash

## submission properties

#SBATCH --partition=iob_p
#SBATCH --job-name=TOBIAS_footprinting
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=50gb
#SBATCH --output=LOG_tobias_postprocessing.%j.log
#SBATCH --error=LOG_tobias_postprocessing.%j.err
#SBATCH --array=1-20

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc

# files
config=celltypes.txt

# functions
process(){

	# verbose
	echo " - processing footprints for $1 ..."

	# variables
	celltype=$1
	motifs=motifs.txt
	wd=$2	

	# iterate over motifs
	while read i; do
		echo " - process footprints for motif $i (cell type = $celltype)..."
		hap1=$wd/$celltype/BINDetect_output/${i}_${i}/beds/${i}_${i}_${celltype}.hap1_bound.bed
		hap2=$wd/$celltype/BINDetect_output/${i}_${i}/beds/${i}_${i}_${celltype}.hap2_bound.bed
		bedtools intersect -a $hap1 -b $hap2 -v > $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap1_specific.bed
		bedtools intersect -a $hap2 -b $hap1 -v > $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap2_specific.bed
		bedtools intersect -a $hap1 -b $hap2 -u > $wd/$celltype/BINDetect_output/$celltype.$i.haplotype_shared.bed
		perl -ne 'chomp;my $newl = "$_\t" . "hap1\n"; print "$newl";' $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap1_specific.bed \
			> $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap1_specific.bed1
		perl -ne 'chomp;my $newl = "$_\t" . "hap2\n"; print "$newl";' $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap2_specific.bed \
			> $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap2_specific.bed1
		mv $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap1_specific.bed1 $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap1_specific.bed
		mv $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap2_specific.bed1 $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap2_specific.bed
		cat $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap1_specific.bed $wd/$celltype/BINDetect_output/${i}_${i}/beds/$celltype.$i.hap2_specific.bed > \
			$wd/$celltype/BINDetect_output/$celltype.$i.haplotype_specific.bed
	done < $motifs

	# merge motifs
	cat $wd/$celltype/BINDetect_output/$celltype.*.haplotype_specific.bed | sort -k1,1 -k2,2n - > $wd/$celltype/BINDetect_output/$celltype.All_motifs.haplotype_specific.bed
	cat $wd/$celltype/BINDetect_output/$celltype.*.haplotype_shared.bed | sort -k1,1 -k2,2n - > $wd/$celltype/BINDetect_output/$celltype.All_motifs.haplotype_shared.bed
}
export -f process


# iterate over each cell types
sID=$(awk -v ArrayID=$SLURM_ARRAY_TASK_ID '$2==ArrayID {print $1}' $config)
process $sID $PWD

# finished
