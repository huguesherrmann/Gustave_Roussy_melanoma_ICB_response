#!/bin/bash

#SBATCH --job-name=quantify_peaks
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/quantify_peaks.out

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/h_herrmann/.environment_conda/analyze_cut_n_tag
module load gcc
module load r/4.1.1


featureCounts -T 4 \
    -a /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/macs2/common_peaks.bed.saf \
    -F SAF \
    -
    -o /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/counts/common_peak_feature_counts.txt \
    /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/bam/H3K4me1/D1332T35_H3K4me1.sorted.bam \
    /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/bam/H3K4me1/D1332T43_H3K4me1.sorted.bam \
    /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/bam/H3K4me1/D1332T46_H3K4me1.sorted.bam \
    /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/bam/H3K4me1/D437T40_H3K4me1.sorted.bam \
    /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/bam/H3K4me1/D437T73_H3K4me1.sorted.bam \
    /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/bam/H3K4me1/D437T81_H3K4me1.sorted.bam
