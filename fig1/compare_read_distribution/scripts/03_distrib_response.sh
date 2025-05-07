#!/bin/bash

#SBATCH --job-name=022_get_read_distribution
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/022_get_read_distribution.out

module load python/3.9.1


read_distribution.py -i /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_metabam/gr1234/metabam_responders/metabam_responders_200m.final.out \
    -r /mnt/beegfs/userdata/h_herrmann/mela_ici_response/misc/scripts/022_get_read_distribution/gencode_hg38v24.bed12 \
    > /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/misc/022_get_read_distribution/gr1234_resp.txt

read_distribution.py -i /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_metabam/gide/metabam_non_responders/metabam_non_responders_200m.final.out \
    -r /mnt/beegfs/userdata/h_herrmann/mela_ici_response/misc/scripts/022_get_read_distribution/gencode_hg38v24.bed12 \
    > /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/misc/022_get_read_distribution/gr1234_non_resp.txt
