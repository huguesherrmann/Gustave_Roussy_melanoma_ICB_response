# ......................................................
# Write files with contigs to query
# 05/08/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
stable_regulons <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/stability/all_stable_regulons.tsv"
correspondence <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/regulons/all_correspondence_contigs_regulons.tsv"
annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/merge_and_annotate/gr1234/annotation_merged_contigs_gide_masked_7_5_gr1234.tsv.gz"
n_redundancy <- 80
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/query_cell_lines/input_for_queries/"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
stable_regulons <- read_tsv(stable_regulons, show_col_types = FALSE)

correspondence <- fread(correspondence)

annotation <- fread(annotation)


# ......................................................
#
#   GET MOST STABLE FEATURES ----
#
# ......................................................
pervasive_regulons <- stable_regulons %>% filter(N_redundancy >= n_redundancy & Condition == "non_responder" & Class == "intergenic") %>%
   pull(Regulon)
pervasive_tags <- correspondence[Regulon %in% pervasive_regulons, list(tag, Regulon)]
pervasive_contigs <- merge(pervasive_tags, annotation, by = "tag")

alr_regulons <- stable_regulons %>% filter(N_redundancy >= n_redundancy & Condition == "non_responder" & grepl("ALR", Regulon)) %>%
   pull(Regulon)
alr_tags <- correspondence[Regulon %in% alr_regulons, list(tag, Regulon)]
alr_contigs <- merge(alr_tags, annotation, by = "tag")

# Include responder intergenic contigs
all_pervasive_regulons <- stable_regulons %>% filter(N_redundancy >= n_redundancy & Class == "intergenic")
all_pervasive_tags <- correspondence[Regulon %in% (all_pervasive_regulons %>% pull(Regulon)), list(tag, Regulon)] %>%
   merge(., all_pervasive_regulons, by = "Regulon")
all_pervasive_contigs <- merge(all_pervasive_tags, annotation, by = "tag")


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(pervasive_contigs, paste0(out_dir, "non_responder_top80_stable_pervasive_contig_annotation.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

write.table(alr_contigs, paste0(out_dir, "non_responder_top80_stable_alr_contig_annotation.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)

write.table(all_pervasive_contigs, paste0(out_dir, "all_top80_stable_pervasive_contig_annotation.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
