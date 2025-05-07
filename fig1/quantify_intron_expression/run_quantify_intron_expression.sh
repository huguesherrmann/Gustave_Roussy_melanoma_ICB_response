#!/bin/bash

#SBATCH --job-name=quantify_intron_expression
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem=50gb
#SBATCH --cpus-per-task=8
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/quantify_intron_expression_gr1234.out

module load singularity/3.6.3
module load python/3.9.1

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/quantify_intron_expression/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/quantify_intron_expression/config/gr1234.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 8
