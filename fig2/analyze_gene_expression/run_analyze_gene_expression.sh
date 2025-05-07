#!/bin/bash

#SBATCH --job-name=analyze_gene_expression
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=2
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/analyze_gene_expression.out

module load gcc
module load r/4.1.1

# ----- Parameters ----- #
SNAKEMAKE=/mnt/beegfs/software/snakemake/5.23.0/bin/snakemake
SMK=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/Snakefile
CONFIG=/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_gene_expression/config/gr1234.json
# ----- Parameters ----- #

$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 1 --unlock
$SNAKEMAKE -s $SMK --configfile $CONFIG --cores 2
