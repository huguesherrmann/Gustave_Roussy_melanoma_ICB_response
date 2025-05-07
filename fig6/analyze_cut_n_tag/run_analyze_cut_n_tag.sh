#!/bin/bash

#SBATCH --job-name=analyze_cut_n_tag
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/analyze_cut_n_tag.out

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/h_herrmann/.environment_conda/analyze_cut_n_tag
module load gcc
module load r/4.1.1

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_cut_n_tag/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_cut_n_tag/config/gr1234_h3k4me1.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 4
