#!/bin/bash

## submission properties

#SBATCH --partition=gpu_p
#SBATCH --gres=gpu:A100:1
#SBATCH --job-name=TENSOR_caQTL_mapping
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:00
#SBATCH --mem=950gb
#SBATCH --output=LOG_282_maize_tensorQTL.%j.log
#SBATCH --error=LOG_282_maize_tensorQTL.%j.err

# set env
cd $SLURM_SUBMIT_DIR
source ~/.zshrc
source ~/Software/tensorqtl/venv/bin/activate

# modules
#ml GLib/2.72.1-GCCcore-11.3.0
ml GLib/2.75.0-GCCcore-12.2.0

# vars
acr=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step6_co_accessibility/Merged_Genotype.coACR.sorted.hic_hichip_support.05.14.2024.bed
output=Merged_Genotype_TRANS
dir=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step6_co_accessibility/map_coACRs
cov=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/single_cell_analysis/step10_caQTL_celltypes/matrices/Merged_Genotype.covariates.txt.PCs_only
snp=/scratch/apm25309/single_cell/ATACseq/multiplexed_genotypes/vcf_282/raw_VCF/hmp321_agpv5.v1.imputed_ds.sorted.MAF_hets

# make dir
if [ ! -d $dir/$output ]; then
	mkdir $dir/$output
else
	rm -rf $dir/$output
	mkdir $dir/$output
fi

# get trans qtl
python3 -m tensorqtl \
                $snp $acr $output.05.14.2024.tensorQTL_trans_p_0.0001 \
                -o $dir/$output \
                --covariates $cov \
                --dosages \
                --pval_threshold 0.0001 \
		--batch_size 1000 \
                --mode trans

# finished
echo " "
echo "Finished running coACR mapping ..."
echo " "
