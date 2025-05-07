# ......................................................
# Get introns
# 02/08/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(BUSpaRse))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
gtf <- "/mnt/beegfs/userdata/h_herrmann/preprocess_RNAseq/annotation_genome_ref/gencode_hg38v24_annotation.gtf.gz"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/quantify_intron_expression/introns/"
read_length <- 100


# ......................................................
#
#   GET INTRONS ----
#
# ......................................................
get_velocity_files(X = gtf, 
                   L = read_length, 
                   Genome = BSgenome.Hsapiens.UCSC.hg38,
                   out_path = out_dir, 
                   isoform_action = "separate")
