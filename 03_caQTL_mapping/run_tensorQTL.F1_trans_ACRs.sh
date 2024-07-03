#!/bin/bash

## submission properties

#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:A100:1
#SBATCH --job-name=TENSOR_caQTL_mapping
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=500gb
#SBATCH --output=LOG_282_maize_tensorQTL.%j.log
#SBATCH --error=LOG_282_maize_tensorQTL.%j.err

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc
source ~/Software/tensorqtl/venv/bin/activate

# modules
ml GLib/2.72.1-GCCcore-11.3.0
ml tabix/0.2.6-GCCcore-11.3.0

# input files
path=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step10_caQTL_celltypes/fastQTL_testing
config=$path/celltype_BED_files.txt

# functions
DOcaQTL(){

	# vars
	acr=$1
	output1=$( basename $acr )
	output=$( echo $output1 | cut -d'.' -f1 )
	dir=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step10_caQTL_celltypes/fastQTL_testing/celltype_results_09_2023
	cov=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step10_caQTL_celltypes/matrices/$output.covariates.txt.PCs_only
	snp=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/vcf_282/raw_VCF/hmp321_agpv5.v1.imputed_ds.sorted.MAF_hets
	transACRs=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step22_plot_hybrid_data/B73_trans.282_ACRs.bed

	# filter ACRs by cis-QTL
	zcat $acr | grep '#' - > $acr.header
	zcat $acr | bedtools intersect -a - -b $transACRs -u | cat $acr.header - > $output.new.bed
	bgzip $output.new.bed
	tabix -p bed $output.new.bed.gz

	# trans
        python3 -m tensorqtl \
                $snp $output.new.bed.gz $output.tensorQTL_trans_p_0.0001 \
                -o $dir/$output.TRANS_F1 \
                --covariates $cov \
		--dosages \
		--pval_threshold 0.0001 \
                --mode trans

	# clean-up
	rm $acr.header
	rm $output.new.bed.gz
	rm $output.new.bed.gz.tbi

	
}
export -f DOcaQTL

# run
DOcaQTL /scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step10_caQTL_celltypes/matrices/Merged_Genotype.INT_ACRs.08.08.2023.fastQTL.txt.bed.gz

# finished
echo " "
echo "Finished running FAST_caQTL mapping ..."
echo " "
