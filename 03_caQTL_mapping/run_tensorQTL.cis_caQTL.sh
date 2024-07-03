#!/bin/bash

## submission properties

#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:A100:1
#SBATCH --job-name=TENSOR_caQTL_mapping
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=500gb
#SBATCH --output=LOG_282_maize_tensorQTL.%j.log
#SBATCH --error=LOG_282_maize_tensorQTL.%j.err
#SBATCH --array=1-34

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc
source ~/Software/tensorqtl/venv/bin/activate

# modules
ml GLib/2.72.1-GCCcore-11.3.0

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

	# make dir
	if [ ! -d $dir/$output ]; then
		mkdir $dir/$output
	else
		rm -rf $dir/$output
		mkdir $dir/$output
	fi

	# run with permutations
	python3 -m tensorqtl \
		$snp $acr $output.tensorQTL_perms \
		-o $dir/$output \
		--covariates $cov \
		--mode cis \
		--qvalue_lambda 0.85 \
		--dosages

	# get nominal p-values for all ACR-snp pairs
	python3 -m tensorqtl \
                $snp $acr $output.tensorQTL_nominal \
		-o $dir/$output \
                --covariates $cov \
                --mode cis_nominal \
                --dosages

	# conditional cis-QTL
	python3 -m tensorqtl \
                $snp $acr $output.tensorQTL_conditional \
		-o $dir/$output \
                --covariates $cov \
		--cis_output $dir/$output.tensorQTL_perms.cis_qtl.txt.gz \
                --mode cis_independent \
                --dosages

	# trans
        python3 -m tensorqtl \
                $snp $acr $output.tensorQTL_trans \
                -o $dir/$output \
                --covariates $cov \
		--dosages \
                --mode trans

	
}
export -f DOcaQTL

# run
sID=$(awk -v ArrayID=$SLURM_ARRAY_TASK_ID '$2==ArrayID {print $1}' $config)
DOcaQTL $sID

# finished
echo " "
echo "Finished running tensorQTL caQTL mapping ..."
echo " "
