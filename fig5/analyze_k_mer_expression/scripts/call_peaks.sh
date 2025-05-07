#!/bin/bash

#SBATCH --job-name=macs
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/userdata/h_herrmann/macs3/macs.out

# ----- Parameters ----- #
WORK_DIR="/mnt/beegfs/userdata/h_herrmann/"
BAM=${WORK_DIR}mela_ici_response_results/compare_sequencing_depth/gr1234/metabam_non_responders/metabam_non_responders_200m.bam
C_BAM=${WORK_DIR}mela_ici_response_results/compare_sequencing_depth/gr1234/metabam_responders/metabam_responders_200m.bam
OUT_DIR=${WORK_DIR}macs3/
SINGULARITY="/mnt/beegfs/software/singularity/3.6.3/bin/singularity"
IMG="/mnt/beegfs/userdata/h_herrmann/tools_software/MACS/macs3_3_0_0a6.sif"
# ----- Parameters ----- #

$SINGULARITY exec --bind $WORK_DIR $IMG macs3 callpeak -t $BAM -c $C_BAM -f BAMPE -g hs --outdir ${WORK_DIR}macs3/try2 -q 0.001

#chr start end length absolute_peak_position pileup_height_at_peak_posisition -log10(pvalue) fold_enrichment -log10(qvalue)
