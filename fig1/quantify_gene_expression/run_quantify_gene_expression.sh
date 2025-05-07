#!/bin/bash

#SBATCH --job-name=quantify_gene_expression
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH --cpus-per-task=9
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/quantify_gene_expression_ines.out

module load gcc
module load r/4.1.1

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/quantify_gene_expression/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/quantify_gene_expression/config/ines.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 9
