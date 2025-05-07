# ......................................................
# Get intergenic and ALR annotations and sequences
# 06/06/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

get_annotation_and_sequences_per_class <- function(differential, annotation, correspondence, class, condition, label_fasta = NULL) {
   if (is.null(label_fasta)) label_fasta <- class
   
   filtered_differential <- differential %>% filter(Class == class & Condition == condition)
   filtered_annotation <- merge(filtered_differential, correspondence %>% select(tag, Regulon), by = "Regulon") %>%
      merge(., annotation, by = "tag")
   
   fasta <- character(nrow(filtered_annotation) * 2)
   fasta[c(TRUE, FALSE)] <- paste0(">", label_fasta, "_", 1:nrow(filtered_annotation))
   fasta[c(FALSE, TRUE)] <- filtered_annotation$contig
   
   return(list("annotation" = filtered_annotation, "fasta" = fasta))
}


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Get intergenic and ALR annotations and sequences.")
parser$add_argument("--annotation", type = "character", help = "Path to regulon annotation table.")
parser$add_argument("--correspondence", type = "character", help = "Path to tag - regulon correspondence table.")
parser$add_argument("--differential", type = "character", help = "Path to differential regulons.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

annotation <- args$annotation
correspondence <- args$correspondence
differential <- args$differential
out_dir <- args$out_dir
# correspondence <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234_w_purity/regulons/all_correspondence_contigs_regulons.tsv"
# annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/merge_and_annotate/gr1234/annotation_merged_contigs_gide_masked_7_5_gr1234.tsv.gz"
# differential <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/definitive_differential/filtered_differential_regulons.tsv"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/analyze_k_mer_expression/gr1234_w_purity/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
differential <- read_tsv(differential, show_col_types = FALSE)

correspondence <- fread(correspondence)

annotation <- fread(annotation, select = c("contig", "tag", "chromosome", "start", "end"))


# --------------------------------------------
#
#   GET ANNOTATION AND SEQUENCES
#
# --------------------------------------------
alr <-  get_annotation_and_sequences_per_class(differential %>% filter(grepl("ALR", Regulon)), 
                                               annotation, 
                                               correspondence, 
                                               class = "repeat", 
                                               condition = "NR")
nr_intergenics <-  get_annotation_and_sequences_per_class(differential,
                                                          annotation, 
                                                          correspondence, 
                                                          class = "intergenic", 
                                                          label_fasta = "nr_intergenic",
                                                          condition = "NR")
r_intergenics <-  get_annotation_and_sequences_per_class(differential,
                                                         annotation, 
                                                         correspondence, 
                                                         class = "intergenic", 
                                                         label_fasta = "r_intergenic",
                                                         condition = "R")

id_to_sequence <- data.frame(ID = c(paste0(">repeat_", 1:nrow(alr$annotation)),
                                    paste0(">nr_intergenic_", 1:nrow(nr_intergenics$annotation)),
                                    paste0(">r_intergenic_", 1:nrow(r_intergenics$annotation))),
                             Contig = c(alr$annotation$contig, nr_intergenics$annotation$contig, 
                                        r_intergenics$annotation$contig),
                             Tag = c(alr$annotation$tag, nr_intergenics$annotation$tag, 
                                     r_intergenics$annotation$tag))

all_fasta <- c(alr$fasta, nr_intergenics$fasta, r_intergenics$fasta)
all_annotations <- alr$annotation %>% add_row(nr_intergenics$annotation) %>%
   add_row(r_intergenics$annotation)
 

# --------------------------------------------
#
#   EXPORT
#
# --------------------------------------------
writeLines(alr$fasta, paste0(out_dir, "sequences/NR_ALR_sequences.fa"))
writeLines(nr_intergenics$fasta, paste0(out_dir, "sequences/NR_intergenic_sequences.fa"))
writeLines(r_intergenics$fasta, paste0(out_dir, "sequences/R_intergenic_sequences.fa"))
writeLines(all_fasta, paste0(out_dir, "sequences/all_ALR_and_intergenics_sequences.fa"))

write.table(alr$annotation, paste0(out_dir, "annotations/NR_ALR_annotation.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
write.table(nr_intergenics$annotation, paste0(out_dir, "annotations/NR_intergenics_annotation.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
write.table(r_intergenics$annotation, paste0(out_dir, "annotations/R_intergenics_annotation.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
write.table(all_annotations, paste0(out_dir, "annotations/all_ALR_and_intergenics_annotation.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

write.table(id_to_sequence, paste0(out_dir, "sequences/id_to_sequences.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
