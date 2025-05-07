# ......................................................
# Get number of coding and non-coding genes detected
# 09/11/22
# Hugues HERRMANN
# ......................................................
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Get number of coding and non-coding genes detected.")
parser$add_argument("--tximport_dir", type = "character", help = "Path to directory containing all tximport directories.")
parser$add_argument("--annotation", type = "character", help = "Path to gene annotation file.")
parser$add_argument("--cohort", type = "character", help = "Name of the cohort.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

tximport_dir <- args$tximport_dir
annotation <- args$annotation
cohort <- args$cohort
out_dir <- args$out_dir
# tximport_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_gene_detection/test/tximport/"
# annotation <- "/mnt/beegfs/userdata/h_herrmann/results_and_data/gene_level_analysis/annot_lookup.tsv"
# cohort <- "test"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_gene_detection/test/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
annotation <- annotation <- read_tsv(annotation, show_col_types = FALSE)
pcg <- annotation %>% filter(gene_biotype == "protein_coding") %>% pull(ensembl_gene_id)
non_coding <- annotation %>% filter(gene_biotype != "protein_coding") %>% pull(ensembl_gene_id)


# ......................................................
#
#   LOAD TXIMPORT FILES ----
#
# ......................................................
sample_dirs <- list.dirs(tximport_dir, recursive = FALSE)

full_df <- data.frame(Sample_ID = character(),
                      Coding_expressed = integer(),
                      Non_coding_expressed = integer(),
                      N_reads = character(),
                      Cohort = character())
for (dir in sample_dirs) {
   n_reads <- str_remove(basename(dir), "sample_")
   
   tmp <- read.table(paste0(dir, "/gene_counts.tsv"), row.names = 1)
   pcg_tmp <- tmp[pcg, ]
   non_coding_tmp <- tmp[non_coding, ]
   
   df_tmp <- data.frame(Sample_ID = colnames(pcg_tmp),
                        Coding_expressed = colSums(pcg_tmp > 0),
                        Non_coding_expressed = colSums(non_coding_tmp > 0)) %>%
      mutate(N_reads = n_reads) %>%
      mutate(Cohort = cohort)
   
   full_df <- full_df %>% add_row(df_tmp)
}


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(full_df, paste0(out_dir, "tximport/expressed_genes_per_sampling.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)