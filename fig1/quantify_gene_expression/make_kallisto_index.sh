#!/bin/bash

#SBATCH --job-name=kallisto_index
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --mem=110gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/kallisto_index.out

KALLISTO="/mnt/beegfs/userdata/h_herrmann/tools_software/kallisto/build/src/kallisto"
#$KALLISTO index -i /mnt/beegfs/userdata/h_herrmann/preprocess_RNAseq/kallisto_index/transcriptome_Hg38v42.idx /mnt/beegfs/userdata/h_herrmann/preprocess_RNAseq/transcriptome_reference/transcriptome_Hg38v42.fa.gz

$KALLISTO index -i /mnt/beegfs/userdata/h_herrmann/preprocess_RNAseq/kallisto_index/intron_transcriptome_Hg38v24.idx /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_gene_expression/gr1234/introns/cDNA_introns.fa
