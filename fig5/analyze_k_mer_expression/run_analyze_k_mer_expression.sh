#!/bin/bash

#SBATCH --job-name=analyze_k_mer_expression
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/analyze_k_mer_expression.out

module load gcc
module load r/4.1.1
module load samtools/1.9

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_k_mer_expression/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_k_mer_expression/config/gr1234_w_purity.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1
