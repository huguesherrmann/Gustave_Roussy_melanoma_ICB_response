#!/bin/bash

#SBATCH --job-name=022_get_read_distribution
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --error=/mnt/beegfs/scratch/h_herrmann/log/022_get_read_distribution.out

module load python/3.9.1


read_distribution.py -i /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_metabam/gr1234/metabam_for_metagene/metabam_114m.bam \
    -r /mnt/beegfs/userdata/h_herrmann/mela_ici_response/misc/scripts/022_get_read_distribution/gencode_hg38v24.bed12 \
    > /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/misc/022_get_read_distribution/gr1234_distrib.txt

read_distribution.py -i /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_metabam/gide/metabam_for_metagene/metabam_69m.bam \
    -r /mnt/beegfs/userdata/h_herrmann/mela_ici_response/misc/scripts/022_get_read_distribution/gencode_hg38v24.bed12 \
    > /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/misc/022_get_read_distribution/gide_distrib.txt

read_distribution.py -i /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_metabam/liu/metabam_for_metagene/metabam_49m.bam \
    -r /mnt/beegfs/userdata/h_herrmann/mela_ici_response/misc/scripts/022_get_read_distribution/gencode_hg38v24.bed12 \
    > /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/misc/022_get_read_distribution/liu_distrib.txt

read_distribution.py -i /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_metabam/hugo/metabam_for_metagene/metabam_71m.bam \
    -r /mnt/beegfs/userdata/h_herrmann/mela_ici_response/misc/scripts/022_get_read_distribution/gencode_hg38v24.bed12 \
    > /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/misc/022_get_read_distribution/hugo_distrib.txt

read_distribution.py -i /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_metabam/riaz/metabam_for_metagene/metabam_55m.bam \
    -r /mnt/beegfs/userdata/h_herrmann/mela_ici_response/misc/scripts/022_get_read_distribution/gencode_hg38v24.bed12 \
    > /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/misc/022_get_read_distribution/riaz_distrib.txt

read_distribution.py -i /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_metabam/mgh/metabam_for_metagene/metabam_72m.bam \
    -r /mnt/beegfs/userdata/h_herrmann/mela_ici_response/misc/scripts/022_get_read_distribution/gencode_hg38v24.bed12 \
    > /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/misc/022_get_read_distribution/mgh_distrib.txt

read_distribution.py -i /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/make_metabam/markovits/metabam_for_metagene/metabam_29m.bam \
    -r /mnt/beegfs/userdata/h_herrmann/mela_ici_response/misc/scripts/022_get_read_distribution/gencode_hg38v24.bed12 \
    > /mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/misc/022_get_read_distribution/markovits_distrib.txt
