#!/bin/bash

#SBATCH --job-name=quantify_intron_expression
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/quantify_intron_expression_test.out


# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/quantify_intron_expression/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/quantify_intron_expression/config/test.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1
