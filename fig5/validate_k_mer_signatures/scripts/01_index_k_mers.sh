#!/bin/bash

#SBATCH --job-name=idx
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=200gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/index.out

# ----- Parameters ----- #
WORK_DIR="/mnt/beegfs/scratch/h_herrmann/kmtricks/markovits/"
IDX_DIR=${WORK_DIR}index_kmers_1_3_markovitz
KMERS=${WORK_DIR}kmers_1_3_markovitz.txt
SINGULARITY="/mnt/beegfs/software/singularity/3.6.3/bin/singularity"
IMG="/mnt/beegfs/userdata/h_herrmann/tools_software/kamrat/KaMRaT_26_02_2024.sif"
# ----- Parameters ----- #

mkdir $IDX_DIR
$SINGULARITY exec --bind $WORK_DIR $IMG kamrat index -intab $KMERS -outdir $IDX_DIR -klen 31 -unstrand
