#!/bin/bash

#SBATCH --job-name=compare_alignment_metrics
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --cpus-per-task=6
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/compare_alignment_metrics_markovitz.out

module load gcc
module load r/4.1.1
source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/h_herrmann/.environment_conda/meta_bam
module load star/2.7.1a


# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/compare_alignment_metrics/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/compare_alignment_metrics/config/markovitz.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 6
