# ......................................................
# Group feature by correlation
# 30/01/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(dbscan))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_k_mer_matrix_processing.R")
set.seed(2024)


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Standardize contig counts and then cluster them with k-means.")
parser$add_argument("--counts", type = "character", help = "Path to contig count matrix.")
parser$add_argument("--design", type = "character", help = "Path to design matrix.")
parser$add_argument("--annotation", type = "character", help = "Path to contig annotation table.")
parser$add_argument("--intergenic_regions", type = "character", help = "Path to connected intergenic_regions.")
parser$add_argument("--out_dir", type = "character", help = "Path to output directory.")
parser$add_argument("--correlation", type = "double", default = 0.8, help = "Pearson correlation threshold to group contigs.")
parser$add_argument("--threads", type = "integer", default = 1, help = "Number of cores to be used.")
parser$add_argument("--k_centers_unmapped", type = "integer", default = 100, help = "Number of centers for kmeans clustering of unmapped contigs.")
parser$add_argument("--k_centers_chimeric", type = "integer", default = 5, help = "Number of centers for kmeans clustering of chimeric contigs.")
parser$add_argument("--verbose", type = "character", default = FALSE, help = "Use to print intermediate messages.")
args <- parser$parse_args()

counts <- args$counts
design <- args$design
annotation <- args$annotation
intergenic_regions <- args$intergenic_regions
out_dir <- args$out_dir
correlation <- args$correlation
threads <- args$threads
k_centers_unmapped <- args$k_centers_unmapped
k_centers_chimeric <- args$k_centers_chimeric
verbose <- args$verbose
# counts <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/merge_and_annotate/gr1234/merged_contigs_gide_masked_7_5_gr1234.tsv.gz"
# design <- "/mnt/beegfs/userdata/h_herrmann/design/gr1234/baseline_curated_design_gr1234.tsv"
# annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/merge_and_annotate/gr1234/annotation_merged_contigs_gide_masked_7_5_gr1234.tsv.gz"
# intergenic_regions <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/subgraphs/connected_intergenics.txt"
# out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/get_stable_features_in_k_mer_matrix/gr1234/regulons/"
# correlation <- 0.5
# threads <- 8
# verbose <- TRUE


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

hg38_offsets <- get_chromosome_offset()

# Extract all gene symbols, repeats or intergenic regions whose have been mapped
exons <- annotated_counts[is_exonic == 1 & is_intronic == 0 & !is.na(gene_symbol) & is.na(repeats), ]
exon_list <- exons %>% distinct(gene_symbol) %>%
  pull(gene_symbol)
n_exon <- length(exon_list)

introns <- annotated_counts[is_exonic == 0 & is_intronic == 1 & !is.na(gene_symbol) & is.na(repeats), ]
intron_list <- introns %>% distinct(gene_symbol) %>%
  pull(gene_symbol)
n_intron <- length(intron_list)

repeats <- annotated_counts[(!is.na(repeats) & is_exonic == 0) | (!is.na(repeats) & nb_hit == 0 & is.na(mapped_to)), ] 
repeat_list <- repeats %>% distinct(repeats) %>%
  pull(repeats)
n_repeat <- length(repeat_list)

intergenics_per_region <- merge(intergenic_regions, annotated_counts, by = "tag") %>%
  mutate(Region_ID = paste0(Chromosome, "_", Region), .after = tag) %>%
  mutate(Center = start + round((end - start) / 2)) %>%
  mutate(Abs_location = Center + hg38_offsets[chromosome])
intergenic_list <- unique(intergenics_per_region, by = "Region_ID") %>% pull(Region_ID)
n_intergenic_region <- length(intergenic_list)

chimerics <- annotated_counts[nb_hit == 0 & mapped_to == "human(chimeric)", ]
chimeric_list <- chimerics %>% pull(tag)

unmappeds <- annotated_counts[nb_hit == 0 & is.na(mapped_to) & is.na(repeats), ]
unmapped_list <- unmappeds %>% pull(tag)

rm(annotated_counts)
gc()


# ......................................................
#
#   DEFINE CHIMERIC AND UNMAPPED REGULONS ----
#
# ......................................................
if (verbose) cat(paste("Group chimeric contigs and define regulons with kmeans...\n"))

class <- "chimeric"
chimeric_results <- define_regulons_per_feature(class = class,
                                                feature_list = chimeric_list, 
                                                index = 1, 
                                                counts = chimerics, 
                                                samples_id = samples_id,
                                                distance = distance,
                                                k_chimeric = k_centers_chimeric)
write.table(chimeric_results, paste0(out_dir, "chimeric_regulons.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

rm(chimeric_results)
gc()
if (verbose) cat("Chimeric regulons done.\n")


if (verbose) cat(paste("Group unmapped contigs and define regulons with kmeans...\n"))

class <- "unmapped"
unmapped_results <- define_regulons_per_feature(class = class,
                                                feature_list = unmapped_list, 
                                                index = 1, 
                                                counts = unmappeds, 
                                                samples_id = samples_id,
                                                distance = distance,
                                                k_unmapped = k_centers_unmapped)
write.table(unmapped_results, paste0(out_dir, "unmapped_regulons.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

rm(unmapped_results)
gc()
if (verbose) cat("Unmapped regulons done.\n")


# ......................................................
#
#   DEFINE EXONIC REGULONS IN PARALLEL ----
#
# ......................................................
if (verbose) cat(paste("Setting up parallel architecture...\n"))

cl <- makePSOCKcluster(threads) # outfile = "/mnt/beegfs/scratch/h_herrmann/r_1_1_group_feature_by_correlation.tmp"
registerDoParallel(cl)
tmp <- clusterEvalQ(cl, library(doParallel))
tmp <- clusterEvalQ(cl, library(dbscan))
tmp <- clusterEvalQ(cl, library(data.table))
tmp <- clusterEvalQ(cl, library(tidyverse))
clusterExport(cl, c("define_regulons_per_feature"))

if (verbose) cat(paste("Group exonic contigs and define regulons with DBSCAN and ", threads, "cores...\n"))

class <- "exon"
exon_results <- foreach(i = 1:n_exon, .combine = "bind_results", .inorder = FALSE, .multicombine = FALSE) %dopar% {
  regulons <- define_regulons_per_feature(class = class,
                                          feature_list = exon_list,
                                          index = i,
                                          counts = exons,
                                          samples_id = samples_id,
                                          distance = distance,
                                          min_pts = min_pts)

  regulons
}
stopCluster(cl)

write.table(exon_results, paste0(out_dir, "exon_regulons.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

rm(exon_results)
gc()
if (verbose) cat("Exonic regulons done.\n")


# ......................................................
#
#   DEFINE INTRONIC REGULONS IN PARALLEL ----
#
# ......................................................
if (verbose) cat(paste("Setting up parallel architecture...\n"))

cl <- makePSOCKcluster(threads) # outfile = "/mnt/beegfs/scratch/h_herrmann/r_1_1_group_feature_by_correlation.tmp"
registerDoParallel(cl)
tmp <- clusterEvalQ(cl, library(doParallel))
tmp <- clusterEvalQ(cl, library(dbscan))
tmp <- clusterEvalQ(cl, library(data.table))
tmp <- clusterEvalQ(cl, library(tidyverse))
clusterExport(cl, c("define_regulons_per_feature"))

if (verbose) cat(paste("Group intronic contigs and define regulons with DBSCAN and ", threads, "cores...\n"))

class <- "intron"
intron_results <- foreach(i = 1:n_intron, .combine = "bind_results", .inorder = FALSE, .multicombine = FALSE) %do% {
  regulons <- define_regulons_per_feature(class = class,
                                          feature_list = intron_list,
                                          index = i,
                                          counts = introns,
                                          samples_id = samples_id,
                                          distance = distance,
                                          min_pts = min_pts)

  regulons
}
stopCluster(cl)

write.table(intron_results, paste0(out_dir, "intron_regulons.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

rm(intron_results)
gc()
if (verbose) cat("Intronic regulons done.\n")


# ......................................................
#
#   DEFINE REPEAT REGULONS IN PARALLEL ----
#
# ......................................................
if (verbose) cat(paste("Setting up parallel architecture...\n"))

cl <- makePSOCKcluster(threads) # outfile = "/mnt/beegfs/scratch/h_herrmann/r_1_1_group_feature_by_correlation.tmp"
registerDoParallel(cl)
tmp <- clusterEvalQ(cl, library(doParallel))
tmp <- clusterEvalQ(cl, library(dbscan))
tmp <- clusterEvalQ(cl, library(data.table))
tmp <- clusterEvalQ(cl, library(tidyverse))
clusterExport(cl, c("define_regulons_per_feature"))

if (verbose) cat(paste("Group repeat contigs and define regulons with DBSCAN and ", threads, "cores...\n"))

class <- "repeat"
repeat_results <- foreach(i = 1:n_repeat, .combine = "bind_results", .inorder = FALSE, .multicombine = FALSE) %dopar% {
  regulons <- define_regulons_per_feature(class = class,
                                          feature_list = repeat_list,
                                          index = i,
                                          counts = repeats,
                                          samples_id = samples_id,
                                          distance = distance,
                                          min_pts = min_pts)

  regulons
}
stopCluster(cl)

write.table(repeat_results, paste0(out_dir, "repeat_regulons.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

rm(repeat_results)
gc()
if (verbose) cat("Repeat regulons done.\n")


# ......................................................
#
#   MERGE AND DEFINE INTERGENIC REGULONS IN PARALLEL ----
#
# ......................................................
if (verbose) cat(paste("Setting up parallel architecture...\n"))

cl <- makePSOCKcluster(threads) # outfile = "/mnt/beegfs/scratch/h_herrmann/r_1_1_group_feature_by_correlation.tmp"
registerDoParallel(cl)
tmp <- clusterEvalQ(cl, library(doParallel))
tmp <- clusterEvalQ(cl, library(dbscan))
tmp <- clusterEvalQ(cl, library(data.table))
tmp <- clusterEvalQ(cl, library(tidyverse))
clusterExport(cl, c("define_regulons_per_feature"))

if (verbose) cat(paste("Group intergenic contigs and define regulons with DBSCAN and ", threads, "cores...\n"))

class <- "intergenic"
intergenic_results <- foreach(i = 1:n_intergenic_region, .combine = "bind_results", .inorder = FALSE, .multicombine = FALSE) %dopar% {
  regulons <- define_regulons_per_feature(class = class,
                                          feature_list = intergenic_list,
                                          index = i,
                                          counts = intergenics_per_region,
                                          samples_id = samples_id,
                                          distance = distance,
                                          min_pts = min_pts)

  regulons
}
stopCluster(cl)

write.table(intergenic_results, paste0(out_dir, "intergenic_regulons.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

rm(intergenic_results)
gc()
if (verbose) cat("Intergenic regulons done.\n")
