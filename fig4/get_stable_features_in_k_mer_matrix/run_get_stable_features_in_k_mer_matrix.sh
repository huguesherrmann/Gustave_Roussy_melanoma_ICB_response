#!/bin/bash

#SBATCH --job-name=get_stable_feature_in_k_mer_matrix
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --cpus-per-task=6
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/get_stable_feature_in_k_mer_matrix.out

module load gcc
module load r/4.1.1
module load python/3.9.1

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/config/gr1234_w_purity.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 6
