# ......................................................
# Compute new counts for regulons by summing, averaging or averaging without taking into account zero values contig counts
# of contigs belonging to the same regulon
# 11/04/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_k_mer_matrix_processing.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Compute new counts for regulons.")
parser$add_argument("--counts", type = "character", help = "Path to contig count matrix.")
parser$add_argument("--design", type = "character", help = "Path to design matrix.")
parser$add_argument("--annotation", type = "character", help = "Path to contig annotation table.")
parser$add_argument("--intergenic_regions", type = "character", help = "Path to connected intergenic_regions.")
parser$add_argument("--out_dir", type = "character", help = "Path to regulons directory.")
parser$add_argument("--operation", type = "character", default = "mean", help = "How counts should be aggregated? Either 'mean', 'mean_wo_zeros' or 'sum'.")
args <- parser$parse_args()

counts <- args$counts
design <- args$design
annotation <- args$annotation
intergenic_regions <- args$intergenic_regions
out_dir <- args$out_dir
operation <- args$operation
# counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/merge_and_annotate/gr1234/merged_contigs_gide_masked_7_5_gr1234.tsv.gz"
# design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
# annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/merge_and_annotate/gr1234/annotation_merged_contigs_gide_masked_7_5_gr1234.tsv.gz"
# intergenic_regions <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/subgraphs/connected_intergenics.txt"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/"
# operation <- "sum"


# ......................................................
#
#   LOAD DATA ----
#
# ......................................................
design <- read_tsv(design, show_col_types = FALSE)

annotation <- fread(annotation, 
                    select = c("tag", "contig", "mapped_to", "chromosome", "start", "end", "nb_hit", "gene_symbol", "gene_biotype", "is_exonic", "is_intronic", "repeats"))
counts <- fread(counts, select = c("tag", design$Sample_ID))
annotated_counts <- merge(annotation, counts, by = "tag")

intergenic_regions <- fread(intergenic_regions)
colnames(intergenic_regions) <- c("Chromosome", "tag", "Region") 

rm(annotation)
rm(counts)
gc()

exon_results <- fread(paste0(out_dir, "regulons/exon_regulons.tsv"))
intron_results <- fread(paste0(out_dir, "regulons/intron_regulons.tsv"))
repeat_results <- fread(paste0(out_dir, "regulons/repeat_regulons.tsv"))
intergenic_results <- fread(paste0(out_dir, "regulons/intergenic_regulons.tsv"))
chimeric_results <- fread(paste0(out_dir, "regulons/chimeric_regulons.tsv"))
unmapped_results <- fread(paste0(out_dir, "regulons/unmapped_regulons.tsv"))


# ......................................................
#
#   GET COUNTS TABLES ----
#
# ......................................................
samples_id <- design$Sample_ID

hg38_offsets <- get_chromosome_offset()

# Extract all gene symbols, repeats or intergenic regions whose have been mapped
exons <- annotated_counts[is_exonic == 1 & is_intronic == 0 & !is.na(gene_symbol) & is.na(repeats), ]
introns <- annotated_counts[is_exonic == 0 & is_intronic == 1 & !is.na(gene_symbol) & is.na(repeats), ]
repeats <- annotated_counts[(!is.na(repeats) & is_exonic == 0) | (!is.na(repeats) & nb_hit == 0 & is.na(mapped_to)), ] 
intergenics_per_region <- merge(intergenic_regions, annotated_counts, by = "tag") %>%
  mutate(Region_ID = paste0(Chromosome, "_", Region), .after = tag) %>%
  mutate(Center = start + round((end - start) / 2)) %>%
  mutate(Abs_location = Center + hg38_offsets[chromosome])
chimerics <- annotated_counts[nb_hit == 0 & mapped_to == "human(chimeric)", ]
unmappeds <- annotated_counts[nb_hit == 0 & is.na(mapped_to) & is.na(repeats), ]

rm(annotated_counts)
gc()


# ......................................................
#
#   NAME EACH REGULON AND CALCULATE NEW COUNTS ----
#
# ......................................................
# ---- Exons
class <- "exon"
exon_results <- setDT(exon_results)
full_exon_regulons <- merge(exon_results, exons, by = "tag", sort = FALSE) %>%
  mutate(Class = class, .after = tag) %>%
  mutate(Regulon = paste(Class, Feature, Intra_cluster_feature_id, sep = "_"), .after = tag)

if (operation == "mean") {
  exon_regulon_counts <- full_exon_regulons[, lapply(.SD, mean), by = Regulon, .SDcols = samples_id]  

} else if (grepl("zero", operation)) {
  exon_regulon_counts <- full_exon_regulons[, lapply(.SD, mean_without_zero), by = Regulon, .SDcols = samples_id]

} else if (operation == "sum") {
  exon_regulon_counts <- full_exon_regulons[, lapply(.SD, sum), by = Regulon, .SDcols = samples_id]

} else {
  message("@operation argument should be either 'mean', 'mean_without_zero' or 'sum'.")
}
correspondance_contigs_exon_regulons <- full_exon_regulons[, list(tag, Regulon, Class, Feature, Intra_cluster_feature_id)]

rm(exon_results)
rm(full_exon_regulons)
gc()


# ---- Introns
class <- "intron"
intron_results <- setDT(intron_results)
full_intron_regulons <- merge(intron_results, introns, by = "tag", sort = FALSE) %>%
  mutate(Class = class, .after = tag) %>%
  mutate(Regulon = paste(Class, Feature, Intra_cluster_feature_id, sep = "_"), .after = tag)

if (operation == "mean") {
  intron_regulon_counts <- full_intron_regulons[, lapply(.SD, mean), by = Regulon, .SDcols = samples_id]
  
} else if (grepl("zero", operation)) {
  intron_regulon_counts <- full_intron_regulons[, lapply(.SD, mean_without_zero), by = Regulon, .SDcols = samples_id]
  
} else if (operation == "sum") {
  intron_regulon_counts <- full_intron_regulons[, lapply(.SD, sum), by = Regulon, .SDcols = samples_id]
}
correspondance_contigs_intron_regulons <- full_intron_regulons[, list(tag, Regulon, Class, Feature, Intra_cluster_feature_id)]

rm(full_intron_regulons)
rm(intron_results)
gc()


# ---- Repeats
class <- "repeat"
repeat_results <- setDT(repeat_results)
full_repeat_regulons <- merge(repeat_results, repeats, by = "tag", sort = FALSE) %>%
  mutate(Class = class, .after = tag) %>%
  mutate(Regulon = paste(Class, Feature, Intra_cluster_feature_id, sep = "_"), .after = tag)

if (operation == "mean") {
  repeat_regulon_counts <- full_repeat_regulons[, lapply(.SD, mean), by = Regulon, .SDcols = samples_id]
  
} else if (grepl("zero", operation)) {
  repeat_regulon_counts <- full_repeat_regulons[, lapply(.SD, mean_without_zero), by = Regulon, .SDcols = samples_id]
  
} else if (operation == "sum") {
  repeat_regulon_counts <- full_repeat_regulons[, lapply(.SD, sum), by = Regulon, .SDcols = samples_id]
}
correspondance_contigs_repeat_regulons <- full_repeat_regulons[, list(tag, Regulon, Class, Feature, Intra_cluster_feature_id)]

rm(full_repeat_regulons)
rm(repeat_results)
gc()


# ---- Intergenics
class <- "intergenic"
intergenic_results <- setDT(intergenic_results)
# Give an ID to each regulon
intergenic_results <- intergenic_results %>% mutate(Region = str_split(Feature, "_") %>% map_chr(., 2), .after = Feature) %>% #%>% mutate(Class = class, .after = tag) %>%
  mutate(Regulon = paste(Feature, Intra_cluster_feature_id, sep = "_"), .after = tag) %>%
  setDT()
# Add genomic position to each regulon
tmp_intergenic_results <- merge(intergenics_per_region[, list(tag, chromosome, start, end)], intergenic_results, by = "tag")
tmp_intergenic_results <- tmp_intergenic_results[order(chromosome, Intra_cluster_feature_id, start, end)]
# ----
# Group by $Regulon, if there are more than 2 instances per group -> output first and last instances
# Then add the position $start of the first instance to the last instance
# Then delete the first instance (because there is no new information, the position $start of the first instance was already known)
tmp_intergenic_results <- tmp_intergenic_results[, if (.N == 1) .SD else .SD[c(1, .N)], by = Regulon, .SDcols = c("chromosome", "start", "end", "Region", "Intra_cluster_feature_id")]
tmp_intergenic_results <- tmp_intergenic_results %>% group_by(Regulon) %>% 
  mutate(start = start[match(Regulon, Regulon)]) %>%
  setDT()
final_regulon_names <- tmp_intergenic_results[, if (.N == 1) .SD else .SD[.N], by = Regulon] %>%
  mutate(Genomic_location = paste0(chromosome, ":", start, "-", end)) %>%
  mutate(Regulon_final_name = paste0(Genomic_location, "_", Region, "_", Intra_cluster_feature_id)) %>%
  .[, list(Regulon, Genomic_location, Regulon_final_name)]
# ----
# final_regulon_names <- tmp_intergenic_results[, .SD[1], by = Regulon, .SDcols = c("chromosome", "start", "end", "Region", "Intra_cluster_feature_id")] %>%
#   mutate(Genomic_location = paste0(chromosome, ":", start, "-", end)) %>%
#   mutate(Regulon_final_name = paste0(Genomic_location, "_", Region, "_", Intra_cluster_feature_id)) %>%
#   .[, list(Regulon, Genomic_location, Regulon_final_name)]

# Merge the 2 data.tables
full_intergenic_results <- merge(final_regulon_names, intergenic_results, by = "Regulon") %>%
  merge(., intergenics_per_region, by = "tag") %>%
  select(c("tag", "Regulon_final_name", "Genomic_location", "Intra_cluster_feature_id", all_of(samples_id))) %>%
  rename(Regulon = "Regulon_final_name") %>%
  rename(Feature = "Genomic_location")

if (operation == "mean") {
  intergenic_regulons_counts <- full_intergenic_results[, lapply(.SD, mean), by = Regulon, .SDcols = samples_id]
  
} else if (grepl("zero", operation)) {
  intergenic_regulons_counts <- full_intergenic_results[, lapply(.SD, mean_without_zero), by = Regulon, .SDcols = samples_id]
  
} else if (operation == "sum") {
  intergenic_regulons_counts <- full_intergenic_results[, lapply(.SD, sum), by = Regulon, .SDcols = samples_id]
}
correspondance_contigs_intergenic_regulons <- full_intergenic_results[, list(tag, Regulon, Feature, Intra_cluster_feature_id)] %>%
  mutate(Class = class, .after = Regulon)

rm(intergenics_per_region)
rm(final_regulon_names)
rm(tmp_intergenic_results)
rm(intergenic_results)
gc()


# ---- Chimerics
class <- "chimeric"
chimeric_results <- setDT(chimeric_results)
full_chimeric_regulons <- merge(chimeric_results, chimerics, by = "tag", sort = FALSE) %>%
  mutate(Class = class, .after = tag) %>%
  mutate(Regulon = paste(Class, Feature, Intra_cluster_feature_id, sep = "_"), .after = tag)

if (operation == "mean") {
  chimeric_regulon_counts <- full_chimeric_regulons[, lapply(.SD, mean), by = Regulon, .SDcols = samples_id]  
  
} else if (grepl("zero", operation)) {
  chimeric_regulon_counts <- full_chimeric_regulons[, lapply(.SD, mean_without_zero), by = Regulon, .SDcols = samples_id]
  
} else if (operation == "sum") {
  chimeric_regulon_counts <- full_chimeric_regulons[, lapply(.SD, sum), by = Regulon, .SDcols = samples_id]
  
} else {
  message("@operation argument should be either 'mean', 'mean_without_zero' or 'sum'.")
}
correspondance_contigs_chimeric_regulons <- full_chimeric_regulons[, list(tag, Regulon, Class, Feature, Intra_cluster_feature_id)] %>%
  mutate(Intra_cluster_feature_id = as.character(Intra_cluster_feature_id))


# ---- Unmapped
class <- "unmapped"
unmapped_results <- setDT(unmapped_results)
full_unmapped_regulons <- merge(unmapped_results, unmappeds, by = "tag", sort = FALSE) %>%
  mutate(Class = class, .after = tag) %>%
  mutate(Regulon = paste(Class, Feature, Intra_cluster_feature_id, sep = "_"), .after = tag)

if (operation == "mean") {
  unmapped_regulon_counts <- full_unmapped_regulons[, lapply(.SD, mean), by = Regulon, .SDcols = samples_id]  
  
} else if (grepl("zero", operation)) {
  unmapped_regulon_counts <- full_unmapped_regulons[, lapply(.SD, mean_without_zero), by = Regulon, .SDcols = samples_id]
  
} else if (operation == "sum") {
  unmapped_regulon_counts <- full_unmapped_regulons[, lapply(.SD, sum), by = Regulon, .SDcols = samples_id]
  
} else {
  message("@operation argument should be either 'mean', 'mean_without_zero' or 'sum'.")
}
correspondance_contigs_unmapped_regulons <- full_unmapped_regulons[, list(tag, Regulon, Class, Feature, Intra_cluster_feature_id)] %>%
  mutate(Intra_cluster_feature_id = as.character(Intra_cluster_feature_id))


# ......................................................
#
#   MERGE TABLES ----
#
# ......................................................
all_regulon_counts <- exon_regulon_counts %>% add_row(intron_regulon_counts) %>%
  add_row(repeat_regulon_counts) %>%
  add_row(intergenic_regulons_counts) %>%
  add_row(chimeric_regulon_counts) %>%
  add_row(unmapped_regulon_counts)
all_correspondences <- correspondance_contigs_exon_regulons %>% add_row(correspondance_contigs_intron_regulons) %>%
  add_row(correspondance_contigs_repeat_regulons) %>%
  add_row(correspondance_contigs_intergenic_regulons) %>%
  add_row(correspondance_contigs_chimeric_regulons) %>%
  add_row(correspondance_contigs_unmapped_regulons)


# ......................................................
#
#   EXPORT ----
#
# ......................................................
write.table(all_regulon_counts, paste0(out_dir, "regulons/all_regulon_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(all_correspondences, paste0(out_dir, "regulons/all_correspondence_contigs_regulons.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
