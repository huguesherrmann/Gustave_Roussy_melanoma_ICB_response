#!/bin/bash

#SBATCH --job-name=query
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/query_markovits.out

# ----- Parameters ----- #
WORK_DIR="/mnt/beegfs/scratch/h_herrmann/kmtricks/mgh/"
IDX_DIR=${WORK_DIR}index_kmers_1_2_cpb_mgh
FASTA=${WORK_DIR}all_ALR_and_intergenics_sequences.fa
OUT_KMERS=${WORK_DIR}query_all_ALR_and_intergenics_sequences.fa.txt
SINGULARITY="/mnt/beegfs/software/singularity/3.6.3/bin/singularity"
IMG="/mnt/beegfs/userdata/h_herrmann/tools_software/kamrat/KaMRaT_26_02_2024.sif"
# ----- Parameters ----- #

$SINGULARITY exec --bind $WORK_DIR $IMG kamrat query -idxdir $IDX_DIR -fasta $FASTA -toquery mean -outpath $OUT_KMERS
