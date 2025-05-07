# ......................................................
# Estimate correlations intra and inter feature for DBSCAN clustering
# 08/04/24
# Hugues HERRMANN
# ......................................................
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
source("/mnt/beegfs/userdata/h_herrmann/mela_ici_response/get_stable_features_in_k_mer_matrix/scripts/functions_for_k_mer_matrix_processing.R")


# ......................................................
#
#   PARSE ARGUMENTS ----
#
# ......................................................
parser <- ArgumentParser(description = "Estimate correlations intra and inter feature for DBSCAN clustering.")
parser$add_argument("--contigs", type = "character", help = "Path to contig counts table.")
parser$add_argument("--out_dir", type = "character", help = "Path to directory for outputs.")
args <- parser$parse_args()

contigs <- args$contigs
out_dir <- args$out_dir
contigs <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/mask_k_mer_matrix/gr123/merged_contig_sorted_gide_masked_kmers_6_5_kmers.tsv.gz"
design <- "/mnt/beegfs/userdata/h_herrmann/design/gr123/baseline_curated_design_gr123.tsv"
#annotation <- "/mnt/beegfs/scratch/h_herrmann/gide_mask/test_annotation_merged_contig_sorted_gide_masked_kmers_6_5_kmers.tsv"
annotation <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/mask_k_mer_matrix/gr123/annotation_merged_contig_sorted_gide_masked_kmers_6_5_kmers.tsv.gz"
intergenic_regions <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/fit_k_mer_based_models/gr123/connected_intergenics/connected.txt"
out_dir <- "/mnt/beegfs/userdata/h_herrmann/mela_ici_response_results/process_k_mer_matrix/gr123/correlations_intra_inter_feature/"


# --------------------------------------------
#
#   LOAD DATA
#
# --------------------------------------------
design <- read_tsv(design, show_col_types = FALSE)

annotation <- fread(annotation, 
                    select = c("tag", "contig", "chromosome", "start", "end", "gene_symbol", "gene_biotype", "is_exonic", "is_intronic", "repeats"))
contigs <- fread(contigs, select = c("tag", design$Sample_ID))
annotated_counts <- merge(annotation, contigs, by = "tag")

intergenic_regions <- fread(intergenic_regions)
colnames(intergenic_regions) <- c("Chromosome", "tag", "Region") 

rm(annotation)
rm(contigs)
gc()


# --------------------------------------------
#
#   EXTRACT CLASSES
#
# --------------------------------------------
# Extract all gene symbols, repeats or intergenic regions whose have been mapped
exons <- annotated_counts[is_exonic == 1 & is_intronic == 0 & !is.na(gene_symbol) & is.na(repeats), ]
exon_list <- exons %>% distinct(gene_symbol) %>%
  pull(gene_symbol)

introns <- annotated_counts[is_exonic == 0 & is_intronic == 1 & !is.na(gene_symbol) & is.na(repeats), ]
intron_list <- introns %>% distinct(gene_symbol) %>%
  pull(gene_symbol)

repeats <- annotated_counts[!is.na(repeats) & is_exonic == 0, ] 
repeat_list <- repeats %>% distinct(repeats) %>%
  pull(repeats)

hg38_offsets <- get_chromosome_offset()
intergenics_per_region <- merge(intergenic_regions, annotated_counts, by = "tag") %>%
  mutate(Region_ID = paste0(Chromosome, "_", Region), .after = tag) %>%
  mutate(Center = start + round((end - start) / 2)) %>%
  mutate(Abs_location = Center + hg38_offsets[chromosome])
intergenic_list <- unique(intergenics_per_region, by = "Region_ID") %>% pull(Region_ID)


# --------------------------------------------
#
#   SAMPLE EACH CLASS
#
# --------------------------------------------
# ---------------
# Exons
set.seed(2024)

exon_sample <- sample(exon_list, size = 100, replace = FALSE)
tag_exon_sample <- exons[gene_symbol %in% exon_sample, ] %>% select(tag, gene_symbol)
sampled_exons <- exons[gene_symbol %in% exon_sample, ] %>%
  column_to_rownames("tag") %>%
  select(design$Sample_ID) %>%
  t_df()

correlations <- cor(sampled_exons, method = "pearson")
correlations[upper.tri(correlations, diag = TRUE)] <- NA

correlation_df <- data.frame(Feature = character(),
                             Intra = numeric(),
                             Inter = numeric(),
                             N_contigs = numeric(),
                             stringsAsFactors = FALSE)
# Extract correlations between intra and inter features
indexes <- seq(1, nrow(tag_exon_sample))
for (exon_i in exon_sample) {
  tag_exon_i <- which(tag_exon_sample$gene_symbol == exon_i)
  other_tag_other_exon_i <- setdiff(indexes, tag_exon_i)
  
  intra_cor <- correlations[tag_exon_i, tag_exon_i] %>% mean(., na.rm = TRUE)
  inter_cor <- correlations[tag_exon_i, other_tag_other_exon_i] %>% mean(., na.rm = TRUE)
  
  correlation_df[nrow(correlation_df) + 1, ] <- c(exon_i, intra_cor, inter_cor, length(tag_exon_i))
}
exon_correlation_df <- correlation_df %>% 
  mutate(Intra = as.numeric(Intra), Inter = as.numeric(Inter)) %>% 
  select(-N_contigs) %>%
  pivot_longer(!Feature, names_to = "Statistic", values_to = "Correlation")

# ---------------
# Introns
set.seed(2024)

intron_sample <- sample(intron_list, size = 50, replace = FALSE)
tag_intron_sample <- introns[gene_symbol %in% intron_sample, ] %>% select(tag, gene_symbol)
sampled_introns <- introns[gene_symbol %in% intron_sample, ] %>%
  column_to_rownames("tag") %>%
  select(design$Sample_ID) %>%
  t_df()

correlations <- cor(sampled_introns, method = "pearson")
correlations[upper.tri(correlations, diag = TRUE)] <- NA

correlation_df <- data.frame(Feature = character(),
                             Intra = numeric(),
                             Inter = numeric(),
                             N_contigs = numeric(),
                             stringsAsFactors = FALSE)
# Extract correlations between intra and inter features
indexes <- seq(1, nrow(tag_intron_sample))
for (intron_i in intron_sample) {
  tag_intron_i <- which(tag_intron_sample$gene_symbol == intron_i)
  other_tag_other_intron_i <- setdiff(indexes, tag_intron_i)
  
  intra_cor <- correlations[tag_intron_i, tag_intron_i] %>% mean(., na.rm = TRUE)
  inter_cor <- correlations[tag_intron_i, other_tag_other_intron_i] %>% mean(., na.rm = TRUE)
  
  correlation_df[nrow(correlation_df) + 1, ] <- c(intron_i, intra_cor, inter_cor, length(tag_intron_i))
}
intron_correlation_df <- correlation_df %>% 
  mutate(Intra = as.numeric(Intra), Inter = as.numeric(Inter)) %>%
  select(-N_contigs) %>%
  pivot_longer(!Feature, names_to = "Statistic", values_to = "Correlation")

# ---------------
# Repeats
set.seed(2024)

repeat_sample <- sample(repeat_list, size = 50, replace = FALSE)
tag_repeat_sample <- repeats[repeats %in% repeat_sample, ] %>% select(tag, repeats) %>%
  .[sample(1:nrow(.), 10000), ]
sampled_repeats <- repeats[repeats %in% repeat_sample, ] %>%
  .[sample(1:nrow(.), 10000), ] %>%
  column_to_rownames("tag") %>%
  select(design$Sample_ID) %>%
  t_df()

correlations <- cor(sampled_repeats, method = "pearson")
correlations[upper.tri(correlations, diag = TRUE)] <- NA

correlation_df <- data.frame(Feature = character(),
                             Intra = numeric(),
                             Inter = numeric(),
                             N_contigs = numeric(),
                             stringsAsFactors = FALSE)
# Extract correlations between intra and inter features
indexes <- seq(1, nrow(tag_repeat_sample))
for (repeat_i in repeat_sample) {
  tag_repeat_i <- which(tag_repeat_sample$repeats == repeat_i)
  other_tag_other_repeat_i <- setdiff(indexes, tag_repeat_i)
  
  intra_cor <- correlations[tag_repeat_i, tag_repeat_i] %>% mean(., na.rm = TRUE)
  inter_cor <- correlations[tag_repeat_i, other_tag_other_repeat_i] %>% mean(., na.rm = TRUE)
  
  correlation_df[nrow(correlation_df) + 1, ] <- c(repeat_i, intra_cor, inter_cor, length(tag_repeat_i))
}
repeat_correlation_df <- correlation_df %>% 
  mutate(Intra = as.numeric(Intra), Inter = as.numeric(Inter)) %>%
  select(-N_contigs) %>%
  pivot_longer(!Feature, names_to = "Statistic", values_to = "Correlation")

# ---------------
# Intergenics
set.seed(2024)

intergenic_sample <- sample(intergenic_list, size = 100, replace = FALSE)
tag_intergenic_sample <- intergenics_per_region[Region_ID %in% intergenic_sample, ] %>% select(tag, Region_ID)
sampled_intergenics <- intergenics_per_region[Region_ID %in% intergenic_sample, ] %>%
  column_to_rownames("tag") %>%
  select(design$Sample_ID) %>%
  t_df()

correlations <- cor(sampled_intergenics, method = "pearson")
correlations[upper.tri(correlations, diag = TRUE)] <- NA

correlation_df <- data.frame(Feature = character(),
                             Intra = numeric(),
                             Inter = numeric(),
                             N_contigs = numeric(),
                             stringsAsFactors = FALSE)
# Extract correlations between intra and inter features
indexes <- seq(1, nrow(tag_intergenic_sample))
for (intergenic_i in intergenic_sample) {
  tag_intergenic_i <- which(tag_intergenic_sample$Region_ID == intergenic_i)
  other_tag_other_intergenic_i <- setdiff(indexes, tag_intergenic_i)
  
  intra_cor <- correlations[tag_intergenic_i, tag_intergenic_i] %>% mean(., na.rm = TRUE)
  inter_cor <- correlations[tag_intergenic_i, other_tag_other_intergenic_i] %>% mean(., na.rm = TRUE)
  
  correlation_df[nrow(correlation_df) + 1, ] <- c(intergenic_i, intra_cor, inter_cor, length(tag_intergenic_i))
}
intergenic_correlation_df <- correlation_df %>% 
  mutate(Intra = as.numeric(Intra), Inter = as.numeric(Inter)) %>%
  select(-N_contigs) %>%
  pivot_longer(!Feature, names_to = "Statistic", values_to = "Correlation")


# ......................................................
#
#   PLOT CORRELATIONS ----
#
# ......................................................
# exon_density_plot <- ggplot(exon_correlation_df, aes(x = Correlation, fill = Statistic)) +
#   geom_density(alpha = 0.8) +
#   geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
#   scale_fill_manual(values = c("Inter" = "#F0E442", "Intra" = "#009E73")) +
#   lims(x = c(-1, 1)) +
#   theme_classic() +
#   labs(x = "Pearson correlation", y = "Density") +
#   theme(legend.position = "top", axis.text = element_text(size = 14), text = element_text(size = 14))
# 
# intron_density_plot <- ggplot(intron_correlation_df, aes(x = Correlation, fill = Statistic)) +
#   geom_density(alpha = 0.8) +
#   geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
#   scale_fill_manual(values = c("Inter" = "#F0E442", "Intra" = "#009E73")) +
#   lims(x = c(-1, 1)) +
#   theme_classic() +
#   labs(x = "Pearson correlation", y = "Density") +
#   theme(legend.position = "top", axis.text = element_text(size = 14), text = element_text(size = 14))

all_correlation <- exon_correlation_df %>% mutate(Class = "Exon") %>%
  add_row(intron_correlation_df %>% mutate(Class = "Intron")) %>%
  add_row(repeat_correlation_df %>% mutate(Class = "Repeat")) %>%
  add_row(intergenic_correlation_df %>% mutate(Class = "Intergenic"))

correlation_density_per_class_plot <- ggplot(all_correlation, aes(x = Correlation, fill = Statistic)) +
  geom_density(alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
  scale_fill_manual(values = c("Inter" = "#F0E442", "Intra" = "#009E73")) +
  lims(x = c(-1, 1)) +
  facet_wrap( ~ Class, ncol = 4) +
  theme_classic() +
  labs(x = "Pearson correlation", y = "Density") +
  theme(legend.position = "top", axis.text = element_text(size = 12), text = element_text(size = 12))


# ......................................................
#
#   FOCUS ON 1 FEATURE ----
#
# ......................................................
exon_i <- exon_sample[1]
tag_exon_i <- which(tag_exon_sample$gene_symbol == exon_i)
exon_i_counts <- sampled_exons %>% select(all_of(tag_exon_i))
colnames(exon_i_counts) <- paste0("contig_", 1:ncol(exon_i_counts))

exon_pairs <- ggpairs(exon_i_counts, 
        columns = 1:5, 
        upper = list(continuous = "points"), 
        lower = list(continuous = "cor"), 
        diag  = list(continuous = "barDiag"),
        progress = FALSE) +
  theme_bw() +
  labs(title = exon_i)


# ......................................................
#
#   EXPORT RESULTS ----
#
# ......................................................
ggsave(paste0(out_dir, "correlation_density_per_class.png"), correlation_density_per_class_plot, height = 5, width = 9)
ggsave(paste0(out_dir, "exon_example_pair_plot.png"), exon_pairs)