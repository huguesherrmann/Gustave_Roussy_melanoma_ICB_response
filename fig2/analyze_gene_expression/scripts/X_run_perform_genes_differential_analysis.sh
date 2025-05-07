#!/bin/bash

#SBATCH --job-name=perform_genes_differential_analysis
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=4
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/perform_genes_differential_analysis.out

module load gcc
module load r/4.1.1

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SCRIPT=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/scripts/perform_genes_differential_analysis.R
# ----- Parameters ----- #

Rscript $SCRIPT
