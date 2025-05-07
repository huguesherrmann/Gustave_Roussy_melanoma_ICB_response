# ......................................................
# Group feature by correlation
# 30/01/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dbscan))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
#source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/prefilter_k_mer_matrix/scripts/functions_for_kmers.R")
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_k_mer_matrix_processing.R")
set.seed(2024)


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/merge_and_annotate/gr1234/merged_contigs_gide_masked_7_5_gr1234.tsv.gz"
design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/merge_and_annotate/gr1234/annotation_merged_contigs_gide_masked_7_5_gr1234.tsv.gz"
intergenic_regions <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/connected.txt"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/unmapped/"
correlation <- 0.6
threads <- 1
verbose <- TRUE


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
design <- read_tsv(design, show_col_types = FALSE)

annotation <- fread(annotation) %>%
                    #select = c("tag", "contig", "chromosome", "start", "end", "mapped_to", "gene_symbol", "gene_biotype", "is_exonic", "is_intronic", "repeats", "nb_hit")) %>%
  filter(nb_hit == 0)
counts <- fread(counts, select = c("tag", design$Sample_ID))
annotated_counts <- merge(annotation, counts, by = "tag")


rm(annotation)
rm(counts)
gc()


# ......................................................
#
#   PREPARE FOR DEFINING REGULONS ----
#
# ......................................................
if (threads >= 2) threads <- threads - 1

samples_id <- design$Sample_ID
col_names_counts <- c("Regulon", samples_id)
n_sample <- length(samples_id)
distance <- sqrt((1 - correlation) * (n_sample * 2))
min_pts <- 2


# Extract all gene symbols, repeats or intergenic regions whose have been mapped
unmapped <- annotated_counts
unmapped_list <- unmapped %>% pull(gene_symbol)
n_unmapped <- nrow(unmapped)


# ......................................................
#
#   DEFINE EXONIC REGULONS IN PARALLEL ----
#
# ......................................................
chimeric <- unmapped %>% filter(mapped_to == "human(chimeric)") %>%
  filter(is.na(repeats))
standardized_chimeric_counts <- chimeric %>% column_to_rownames("tag") %>%
  select(all_of(samples_id)) %>%
  t() %>%
  scale() %>%
  t()
kme <- kmeans(standardized_chimeric_counts, centers = 5, nstart = 100, iter.max = 100)
chimeric_df <- kme$cluster %>% as.data.frame() %>%
  rename(Intra_cluster_feature_id = ".") %>%
  rownames_to_column("tag") %>%
  mutate(Feature = "chimeric") %>%
  inner_join(., annotation, by = "tag")
write.table(chimeric_df, "/mnt/beegfs/scratch/h_herrmann/chimeric_annotations.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


no_mapp <- unmapped %>% filter(is.na(mapped_to)) %>%
  filter(is.na(repeats))
standardized_nomapp_counts <- no_mapp %>% column_to_rownames("tag") %>%
  select(all_of(samples_id)) %>%
  t() %>%
  scale() %>%
  t()
kme2 <- kmeans(standardized_nomapp_counts, centers = 100, nstart = 1, iter.max = 100)
no_mapp_df <- kme2$cluster %>% as.data.frame() %>%
  rename(Intra_cluster_feature_id = ".") %>%
  rownames_to_column("tag") %>%
  mutate(Feature = "unmapped") %>%
  inner_join(., annotation, by = "tag")
write.table(no_mapp_df, "/mnt/beegfs/scratch/h_herrmann/unmapped_annotations.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# dbscan <- dbscan(standardized_chimeric_counts, eps = distance, minPts = min_pts)
# # Associate each tag with a cluster ID and create an unique ID for singletons
# regulons <- data.frame(tag = rownames(standardized_chimeric_counts), Intra_cluster_feature_id = as.character(dbscan$cluster)) %>%
#   mutate(Intra_cluster_feature_id = if_else(Intra_cluster_feature_id == 0, paste0("0_", row_number()), Intra_cluster_feature_id)) %>%
#   mutate(Feature = "chimeric", .after = tag)
# 
# standardized_feature_counts <- unmapped %>% column_to_rownames("tag") %>%
#   select(all_of(samples_id)) %>%
#   t() %>%
#   scale() %>%
#   t()
# dbscan <- dbscan(standardized_feature_counts, eps = distance, minPts = min_pts)
# Associate each tag with a cluster ID and create an unique ID for singletons
# regulons <- data.frame(tag = rownames(standardized_feature_counts), Intra_cluster_feature_id = as.character(dbscan$cluster)) %>%
#   mutate(Intra_cluster_feature_id = if_else(Intra_cluster_feature_id == 0, paste0("0_", row_number()), Intra_cluster_feature_id)) %>%
#   mutate(Feature = "unmapped", .after = tag)
# 
# write.table(regulons, "/mnt/beegfs/scratch/h_herrmann/unmmped_regulons.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
