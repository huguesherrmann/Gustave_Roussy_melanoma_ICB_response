# ......................................................
# Make BED file of expressed and manually-selected regions 
# 27/09/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/analyze_cut_n_tag/scripts/functions_for_analyzing_cut_n_tag.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Make BED file of expressed and manually-selected regions.")
parser$add_argument("--metadata", type = "character", help = "Metadata table. include chromosome, start, end, strand and sample expression.")
parser$add_argument("--quantif", type = "character", help = "TSS quantification table.")
parser$add_argument("--n_reads", type = "character", help = "number of lines in fastq table.")
parser$add_argument("--offset_plot", type = "integer", default = 5000, help = "Number of nucleotides around the region of interest for the plot.")
parser$add_argument("--offset_signal", type = "integer", default = 2000, help = "Number of nucleotides around TSS to estimate the signal (for computing TSS enrichment score).")
parser$add_argument("--offset_noise", type = "integer", default = 100, help = "Number of nucleotides at the edges of the TSS signal region to estimate the noise (for computing TSS enrichment score).")
parser$add_argument("--threshold", type = "double", default = 0, help = "Minimal expression to consider the region as expressed.")
parser$add_argument("--mark", type = "character", help = "Mark studied. useful for the name of the output.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

metadata <- args$metadata
quantif <- args$quantif
n_reads <- args$n_reads
offset_plot <- args$offset_plot
offset_signal <- args$offset_signal
offset_noise <- args$offset_noise
threshold <- args$threshold
mark <- args$mark
out_dir <- args$out_dir
# metadata <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/metadata/metadata_regions.tsv"
# quantif <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/featureCounts/H3K4me1/quantif_tss_regions.txt"
# n_reads <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/compare_alignment_metrics/gr1234/n_reads/n_lines_fastq.txt"
# out_dir <-"/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_cut_n_tag/gr1234/"
# offset_plot <- 5000
# offset_signal <- 2000
# offset_noise <- 100
# threshold <- 0.15
# mark <- "H3K4me1"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
metadata <- read_tsv(metadata, show_col_types = FALSE)
quantif <- read_tsv(quantif, skip = 1, show_col_types = FALSE)
n_reads <- read_table(n_reads, show_col_types = FALSE)


# ......................................................
#
#   FORMAT DATA AND GET CPM ----
#
# ......................................................
# Format column names
quantif <- quantif %>% column_to_rownames("Geneid") %>%
   select(all_of(starts_with("/")))
colnames(quantif) <- str_remove(basename(colnames(quantif)), ".bam")

samples <- colnames(quantif)

cpm_quantif <- normalize_cpm(quantif, n_reads, base = 1e6) %>%
   rownames_to_column("Name")
#table(as.matrix(cpm_quantif %>% select(-Name)) > 0.2)

metadata <- inner_join(metadata %>% select(Chromosome, Start, End, Name, TSS, Strand), cpm_quantif, by = "Name")


# ......................................................
#
#   MAKE BED ----
#
# ......................................................
# Make 6 bed files for each sample, the first 2 bed files are used for plotting TSS enrichment:
# 1- one with the coordinates of TSS (+/- @offset) when the region is expressed (>= expression @threshold)
# 2- one with the coordinates of TSS (+/- @offset) when the region is not expressed (< expression @threshold)
# the last 4 bed files are used for computing TSS enrichment score: 3 and 4 are for estimating the signal and 5 and 6 are used for estimating the noise
# TSS enrichment score is computed according to https://www.encodeproject.org/data-standards/terms/#enrichment
# 3- one with the coordinates of TSS (+/- 2000bp) when the region is expressed (>= expression @threshold)
# 4- one with the coordinates of TSS (+/- 2000bp) when the region is not expressed (< expression @threshold)
# 5- one with the coordinates of the edges of the +/- 2000bp window (+/- 100bp) when the region is expressed (>= expression @threshold)
# 6- one with the coordinates of the edges of the +/- 2000bp window (+/- 100bp) when the region is not expressed (< expression @threshold)
#
#                 TSS
#                  |
#  ------------------------------------------
#        |==================|
#      -2000              +2000     => signal
#        |=|              |=|
#      -100               +100      => noise
#
for (sample in samples) {
   # Bed for plotting TSS enrichment
   bed_expression_plot <- make_bed(metadata, 
                                   sample, 
                                   offset_plot = offset_plot, 
                                   expression_threshold = threshold,
                                   type = "expression",
                                   objective = "plot")
   bed_expression_plot <- plyranges::as_granges(bed_expression_plot)
   plyranges::write_bed(bed_expression_plot, paste0(out_dir, "bed/", mark, "/", sample, "_expressed_and_TSS_regions.bed"))
   
   bed_no_expression_plot <- make_bed(metadata,
                                      sample,
                                      offset_plot = offset_plot, 
                                      expression_threshold = threshold,
                                      type = "no_expression",
                                      objective = "plot")
   bed_no_expression_plot <- plyranges::as_granges(bed_no_expression_plot)
   plyranges::write_bed(bed_no_expression_plot, paste0(out_dir, "bed/", mark, "/", sample, "_not_expressed_and_TSS_regions.bed"))
   
   # Bed for computing TSS enrichment score
   bed_expression_signal_score <- make_bed(metadata, 
                                           sample, 
                                           offset_signal = offset_signal, 
                                           expression_threshold = threshold,
                                           type = "expression",
                                           objective = "score_signal")
   write.table(bed_expression_signal_score, 
               paste0(out_dir, "bed/", mark, "/", sample, "_expressed_and_TSS_regions_for_scoring_signal.bed.saf"), 
               sep = "\t", quote = FALSE, row.names = FALSE)
   
   bed_expression_noise_score <- make_bed(metadata, 
                                           sample, 
                                           offset_noise = offset_noise, 
                                           expression_threshold = threshold,
                                           type = "expression",
                                           objective = "score_noise")
   write.table(bed_expression_noise_score, 
               paste0(out_dir, "bed/", mark, "/", sample, "_expressed_and_TSS_regions_for_scoring_noise.bed.saf"), 
               sep = "\t", quote = FALSE, row.names = FALSE)

   bed_no_expression_signal_score <- make_bed(metadata, 
                                              sample, 
                                              offset_signal = offset_signal, 
                                              expression_threshold = threshold,
                                              type = "no_expression",
                                              objective = "score_signal")
   write.table(bed_no_expression_signal_score, 
               paste0(out_dir, "bed/", mark, "/", sample, "_not_expressed_and_TSS_regions_for_scoring_signal.bed.saf"), 
               sep = "\t", quote = FALSE, row.names = FALSE)
   
   bed_no_expression_noise_score <- make_bed(metadata, 
                                             sample, 
                                             offset_noise = offset_noise, 
                                             expression_threshold = threshold,
                                             type = "no_expression",
                                             objective = "score_noise")
   write.table(bed_no_expression_noise_score, 
               paste0(out_dir, "bed/", mark, "/", sample, "_not_expressed_and_TSS_regions_for_scoring_noise.bed.saf"), 
               sep = "\t", quote = FALSE, row.names = FALSE)
}
