#!/bin/bash

#SBATCH --job-name=compare_gene_detection
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=8
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/compare_gene_detection_markovitz.out

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/h_herrmann/.environment_conda/kallisto_tximport
module load gcc
module load r/4.1.1

SMK="/mnt/beegfs/userdata/h_herrmann/mela_ici_response/compare_gene_detection/Snakefile"
CONFIG="/mnt/beegfs/userdata/h_herrmann/mela_ici_response/compare_gene_detection/config/markovitz.json"

/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake -s $SMK --configfile $CONFIG --cores 1 --unlock
/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake -s $SMK --configfile $CONFIG --cores 8
