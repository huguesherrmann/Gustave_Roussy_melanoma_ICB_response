#!/bin/bash

#SBATCH --job-name=compare_intron_distrib
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --cpus-per-task=6
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/compare_intron_distrib.out

module load python/3.9.1

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/compare_read_distribution/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/compare_read_distribution/config/gr1234.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 6
