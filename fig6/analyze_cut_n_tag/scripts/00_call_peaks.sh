#!/bin/bash

#SBATCH --job-name=call_peaks
#SBATCH --ntasks=1
#SBATCH --partition=shortq
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/call_peaks.out

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/h_herrmann/.environment_conda/macs2

BAM_DIR="/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/bam/H3K27me3/"
OUT_DIR="/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/macs2/"

macs2 callpeak --nomodel --extsize 200 --keep-dup all -f BAM -p 0.001 -t ${BAM_DIR}D1332T43_H3K27me3.sorted.bam -n ${OUT_DIR}D1332T43_H3K27me3
macs2 callpeak --nomodel --extsize 200 --keep-dup all -f BAM -p 0.001 -t ${BAM_DIR}D1332T46_H3K27me3.sorted.bam -n ${OUT_DIR}D1332T46_H3K27me3
macs2 callpeak --nomodel --extsize 200 --keep-dup all -f BAM -p 0.001 -t ${BAM_DIR}D437T40_H3K27me3.sorted.bam -n ${OUT_DIR}D437T40_H3K27me3
macs2 callpeak --nomodel --extsize 200 --keep-dup all -f BAM -p 0.001 -t ${BAM_DIR}D437T81_H3K27me3.sorted.bam -n ${OUT_DIR}D437T81_H3K27me3
