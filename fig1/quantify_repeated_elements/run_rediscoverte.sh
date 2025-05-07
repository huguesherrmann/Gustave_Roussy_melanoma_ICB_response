#!/bin/bash

#SBATCH --job-name=quantify_repeated_elements
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6gb
#SBATCH --cpus-per-task=18
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/quantify_repeated_elements_gr1234.out

module load gcc
module load r/4.1.1

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SNAKEFILE=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/quantify_repeated_elements/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/quantify_repeated_elements/config/gr1234.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SNAKEFILE --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SNAKEFILE --configfile $CONFIG --cores 18
