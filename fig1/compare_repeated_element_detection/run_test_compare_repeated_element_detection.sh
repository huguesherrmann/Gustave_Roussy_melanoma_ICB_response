#!/bin/bash

#SBATCH --job-name=compare_repeated_element_detection
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8gb
#SBATCH --cpus-per-task=12
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/compare_repeated_element_detection_test.out

module load gcc
module load r/4.1.1

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/compare_repeated_element_detection/test_snk
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/compare_repeated_element_detection/config/test.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 12
