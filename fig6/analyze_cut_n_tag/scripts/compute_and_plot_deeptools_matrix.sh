#!/bin/bash

#SBATCH --job-name=compute_and_plot_deeptools_matrix
#SBATCH --ntasks=1
#SBATCH --partition=shortq
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/compute_and_plot_deeptools_matrix.out

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/h_herrmann/.environment_conda/deeptools

# ----- Parameters ----- #
BW_DIR=/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/bw/
BED_DIR=/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/bed/
MATRICES_DIR=/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/deeptools_matrices/
PLOT_DIR=/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/profile_plots/
MARK="H3K4me1"
# ----- Parameters ----- #

SAMPLES="D437T73 D1332T43 D1332T35 D1332T46 D437T81"

for sample in $SAMPLES;
do
    echo "Computing and plotting deeptools matrix sample: "$sample
    computeMatrix scale-regions -S ${BW_DIR}${sample}_${MARK}.bw -R ${BED_DIR}${sample}_expressed_and_TSS_regions.bed -m 2000 -o ${MATRICES_DIR}matrix_${MARK}_${sample}_expressed_and_TSS_regions.gz
    computeMatrix scale-regions -S ${BW_DIR}${sample}_${MARK}.bw -R ${BED_DIR}${sample}_not_expressed_and_TSS_regions.bed -m 2000 -o ${MATRICES_DIR}matrix_${MARK}_${sample}_not_expressed_and_TSS_regions.gz

    plotProfile -m ${MATRICES_DIR}/matrix_${MARK}_${sample}_expressed_and_TSS_regions.gz -out ${PLOT_DIR}profile_${MARK}_${sample}_expressed_and_TSS_regions.pdf --startLabel "TSS -5K" --endLabel "TSS +5K" -T "H3K4me1 - expressed regions" --yMin 0 --yMax 200
    plotProfile -m ${MATRICES_DIR}/matrix_${MARK}_${sample}_not_expressed_and_TSS_regions.gz -out ${PLOT_DIR}profile_${MARK}_${sample}_not_expressed_and_TSS_regions.pdf --startLabel "TSS -5K" --endLabel "TSS +5K" -T "H3K4me1 - not expressed regions" --yMin 0 --yMax 200
    plotProfile -m ${MATRICES_DIR}/matrix_${MARK}_${sample}_expressed_and_TSS_regions.gz -out ${PLOT_DIR}profile_${MARK}_${sample}_expressed_and_TSS_regions.png --startLabel "TSS -5K" --endLabel "TSS +5K" -T "H3K4me1 - expressed regions" --yMin 0 --yMax 250
    plotProfile -m ${MATRICES_DIR}/matrix_${MARK}_${sample}_not_expressed_and_TSS_regions.gz -out ${PLOT_DIR}profile_${MARK}_${sample}_not_expressed_and_TSS_regions.png --startLabel "TSS -5K" --endLabel "TSS +5K" -T "H3K4me1 - not expressed regions" --yMin 0 --yMax 250
    echo -e "Done\n"
done
