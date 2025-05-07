#!/bin/bash

#SBATCH --job-name=analyze_wes
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/analyze_wes_gr1234.out

module load gcc
module load r/4.1.1

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_wes/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_wes/config/gr1234.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1
