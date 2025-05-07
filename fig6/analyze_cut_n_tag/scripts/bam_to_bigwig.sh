#!/bin/bash

#SBATCH --job-name=bam_to_bigwig
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --cpus-per-task=6
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/bam_to_bigwig.out

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/h_herrmann/.environment_conda/deeptools


# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/star/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/star/config/gr1234.json
# ----- Parameters ----- #


