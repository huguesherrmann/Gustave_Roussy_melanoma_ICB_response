# ......................................................
# Make BED file of expressed and manually-selected regions 
# 27/09/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_cut_n_tag/scripts/functions_for_analyzing_cut_n_tag.R")
set.seed(2024)


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
# parser <- ArgumentParser(description = "Make BED file of expressed and manually-selected regions.")
# parser$add_argument("--metadata", type = "character", help = "Metadata table. include chromosome, start, end, strand and sample expression.")
# parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
# args <- parser$parse_args()
# 
# metadata <- args$metadata
# out_dir <- args$out_dir
gtf <- "/mnt/beegfs/userdata/h_herrmann/preprocess_RNAseq/annotation_genome_ref/gencode_hg38v24_primary_assembly_annotation.gtf"
out_dir <-"/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/"
n_random <- 1000


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
gtf <- fread("preprocess_RNAseq/annotation_genome_ref/gencode_hg38v24_primary_assembly_annotation.gtf", skip = 5, header = FALSE)


# ......................................................
#
#   SAMPLE GENES ----
#
# ......................................................
gtf <- gtf %>% filter(V3 == "gene") %>%
   mutate(Name = row_number()) %>%
   mutate(TSS = 1) %>%
   select(V1, V4, V5, Name, TSS, V7) %>%
   mutate(V1 = str_remove("chr", V1))
colnames(gtf) <- c("Chromosome", "Start", "End", "Name", "TSS", "Strand")

sample_gtf <- gtf %>% sample_n(n_random)


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(sample_gtf, paste0(out_dir, "metadata/random_genes.bed"), sep = "\t", quote = FALSE, row.names = FALSE)
